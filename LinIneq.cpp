#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

LinCon **LinIneq::li;

int LinIneq::M;
int LinIneq::N;
int LinIneq::nr;

/**
 * initialize the static lists and variables
 * @param M_in input dimension of sp space
 * @param N_in input nr of particles
 * @param nr_in input nr of constraints
 */
void LinIneq::init(int M_in,int N_in,int nr_in){

   M = M_in;
   N = N_in;
   nr = nr_in;

   li = new LinCon * [nr];

   for(int i = 0;i < nr;++i)
      li[i] = new LinCon();

}

/**
 * function that deallocates the static lists
 */
void LinIneq::clear(){

   for(int i = 0;i < nr;++i)
      delete li[i];

   delete [] li;

}

/**
 * constructor of a LinIneq object
 */
LinIneq::LinIneq(){

   proj = new double [nr];

}

/**
 * copy constructor
 * @param li_copy The LinIneq object to be copied
 */
LinIneq::LinIneq(const LinIneq &li_copy){

   this->tr = li_copy.gtr();

   proj = new double [nr];

   for(int i = 0;i < nr;++i)
      proj[i] = li_copy.gproj(i);

}

/**
 * destructor
 */
LinIneq::~LinIneq(){

   delete [] proj;

}

/**
 * @return the number of currently applied linear constraints
 */
int LinIneq::gnr(){

   return nr;

}

/**
 * read and write access to your LinCon object
 * @param i row number
 * @return the entry on index i
 */
LinCon &LinIneq::operator[](int i){

   return *li[i];

}

/**
 * read only access to your LinCon object
 * @param i row number
 * @return the entry on index i
 */
const LinCon &LinIneq::operator[](int i) const{

   return *li[i];

}

/**
 * @return nr of particles
 */
int LinIneq::gN() const{

   return N;

}

/**
 * @return nr of sp orbitals
 */
int LinIneq::gM() const{

   return M;

}

/**
 * @param li_epsi the LinIneq object filled with the Newton direction found by the cg loop.
 * @return the first singularity in the potential along the "epsilon" direction.
 */
double LinIneq::min_end(const LinIneq &li_epsi) const{

   double min_end = -constraint(0)/li_epsi.constraint(0);

   if(min_end < 0.0)
      min_end = 1.0e+15;

   for(int i = 1;i < nr;++i){

      double tmp = -constraint(i)/li_epsi.constraint(i);

      if(tmp < 0.0)
         tmp = 1.0e+15;

      if(tmp < min_end)
         min_end = tmp;

   }

   return min_end;

}

/**
 * Calculate the projection of the input TPM object on the constraint matrices of the elements of LinIneq
 * @param tpm The input TPM.
 */
void LinIneq::fill(const TPM &tpm){

   for(int i = 0;i < nr;++i)
      proj[i] = (li[i]->gI()).ddot(tpm);

   tr = 2.0*tpm.trace()/(N*(N - 1.0));

}

/**
 * @param i the index of the constraint we are interested in.
 * @return the projection on the constaint with index i
 */
double LinIneq::gproj(int i) const {

   return proj[i];

}

/**
 * @return the array containing the projection on the constraints.
 */
double *LinIneq::gproj() {

   return proj;

}

/**
 * @param index the index...
 * @return the value of the constraint labeled with index index
 */
double LinIneq::constraint(int index) const{

   return proj[index] - li[index]->gi()*tr;

}

/**
 * @return the trace of the input TPM scaled with N(N-1)/2
 */
double LinIneq::gtr() const{

   return tr;

}

/**
 * @param c The stepsize along the newtonstep epsilon
 * @param li_epsi The LinIneq object containing the linear constraint info about the step epsilon.
 * @return the value of the linear inequality part of the line search function with abcis c along the Newton step epsilon.
 */
double LinIneq::lsfunc(double c,const LinIneq &li_epsi) const {

   double ward = 0.0;

   for(int i = 0;i < nr;++i)
      ward += li_epsi.constraint(i)/(constraint(i) + c*li_epsi.constraint(i));

   return ward;

}
