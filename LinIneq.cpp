#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

int LinIneq::nr;

double LinIneq::a;
double LinIneq::c;

int LinIneq::M;
int LinIneq::N;

LinCon **LinIneq::li;

double *LinIneq::coef;

/**
 * Static function that initialializes the static members: it allocates the LinCon objects and 
 * initializes them.
 * @param M_in nr of sp orbitals
 * @param N_in nr of particles
 * @param si input SphInt object for construction of subsystem constraints
 */
void LinIneq::init(int M_in,int N_in,const SphInt &si){

   M = M_in;
   N = N_in;

   nr = 4;

   //allocate
   li = new LinCon * [nr];

   for(int i = 0;i < nr;++i)
      li[i] = new LinCon();

   for(int i = 0;i < nr;++i)
      li[i]->fill_Random();
/*
   //make the subsystem objects
   SubSys ss_Be(0,si);
   ss_Be.setBe();

   SubSys ss_B(1,si);
   ss_B.setB();

   li[0]->subcon(ss_Be,0);
   li[1]->subcon(ss_Be,1);
   li[2]->subcon(ss_B,0);
   li[3]->subcon(ss_B,1);
*/
   //what are the coef's of the overlap matrix without the linear constraints:
   init_overlap();

   int n = 2*nr;

   //make the linear system for the inverse overlapmatrix coefficient calculations
   coef = new double [n*n];

   //now make some variables needed for the making of the system

   //overlap
   Matrix I_overlap(nr);

   for(int i = 0;i < nr;++i)
      for(int j = i;j < nr;++j)
         I_overlap(i,j) = li[i]->gI().ddot(li[j]->gI());

   I_overlap.symmetrize();

   //tenslotte
   Matrix I_bar_overlap(nr);

   for(int i = 0;i < nr;++i)
      for(int j = i;j < nr;++j)
         I_bar_overlap(i,j) = li[i]->gI_bar().ddot(li[j]->gI_bar());

   I_bar_overlap.symmetrize();

   for(int k = 0;k < nr;++k){//columns

      for(int i = 0;i < nr;++i){//rows

         coef[k*n + i] = I_overlap(i,k);

         if(i == k)
            coef[k*n + i] += a;

      }

      for(int i = nr;i < 2*nr;++i)//rows
         coef[k*n + i] = 0.25*I_bar_overlap(i - nr,k);

   }

   for(int k = nr;k < 2*nr;++k){//columns

      for(int i = 0;i < nr;++i)//rows
         coef[k*n + i] = 0.0;

      coef[k*n + k - nr] = -4.0*c;

      for(int i = nr;i < 2*nr;++i)
         coef[k*n + i] = 0.0;

      coef[k*n + k] = a - c*(M - 2.0);

   }

   int *ipiv = new int [n];

   int lwork = 3*n;

   double *work = new double [lwork];

   int info;

   //done, now invert the mofo
   dgetrf_(&n,&n,coef,&n,ipiv,&info);

   dgetri_(&n,coef,&n,ipiv,work,&lwork,&info);

   delete [] ipiv;
   delete [] work;

}

/**
 * delete the statics
 */
void LinIneq::clear(){

   for(int i = 0;i < nr;++i)
      delete li[i];

   delete [] li;

   delete [] coef;

}

/**
 * Will calculate the parameters needed for the overlapmatrix-map: a and c.
 */
void LinIneq::init_overlap(){

   a = 1.0;
   c = 0.0;

#ifdef __Q_CON

   a += 1.0;
   c += (2.0*N - M)/((N - 1.0)*(N - 1.0));

#endif

#ifdef __G_CON

   a += 4.0;
   c += (2.0*N - M - 2.0)/((N - 1.0)*(N - 1.0));

#endif

#ifdef __T1_CON

   a += M - 4.0;
   c -= (M*M + 2.0*N*N - 4.0*M*N - M + 8.0*N - 4.0)/( 2.0*(N - 1.0)*(N - 1.0) );

#endif

#ifdef __T2_CON

   a += 5.0*M - 8.0;
   c += (2.0*N*N + (M - 2.0)*(4.0*N - 3.0) - M*M)/(2.0*(N - 1.0)*(N - 1.0));

#endif

}

/**
 * constructor of a LinIneq object
 */
LinIneq::LinIneq(){

   proj = new double [nr];
   proj_bar = new double [nr];

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

   proj_bar = new double [nr];

   for(int i = 0;i < nr;++i)
      proj_bar[i] = li_copy.gproj_bar(i);

}

/**
 * destructor
 */
LinIneq::~LinIneq(){

   delete [] proj;
   delete [] proj_bar;

}

/**
 * overload the equality operator
 * @param li_copy The LinIneq you want to be copied into this
 */
LinIneq &LinIneq::operator=(const LinIneq &li_copy){

   int inc = 1;

   dcopy_(&nr,li_copy.proj,&inc,proj,&inc);

   return *this;

}

/**
 * Make all the numbers in your LinIneq projection object array equal to the number a.
 * @param a the number
 */
LinIneq &LinIneq::operator=(double a){

   for(int i = 0;i < nr;++i)
      proj[i] = a;

   return *this;

}

/**
 * overload the += operator
 * @param li_pl The LinIneq object you want to add to this
 */
LinIneq &LinIneq::operator+=(const LinIneq &li_pl){

   int inc = 1;
   double alpha = 1.0;

   daxpy_(&nr,&alpha,li_pl.proj,&inc,proj,&inc);

   return *this;

}

/**
 * overload the -= operator
 * @param li_pl The LinIneq object you want to deduct from this
 */
LinIneq &LinIneq::operator-=(const LinIneq &li_pl){

   int inc = 1;
   double alpha = -1.0;

   daxpy_(&nr,&alpha,li_pl.proj,&inc,proj,&inc);

   return *this;

}

/**
 * the daxpy function, add another object scaled with alpha.
 * @param alpha the scaling factor
 * @param li_pl The LinIneq object you want to add to this
 */
LinIneq &LinIneq::daxpy(double alpha,const LinIneq &li_pl){

   int inc = 1;

   daxpy_(&nr,&alpha,li_pl.proj,&inc,proj,&inc);

   return *this;

}

/**
 * static function returning the number of constraints
 */
int LinIneq::gnr() {

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
 * Calculate the projection of the input TPM object on the constraint matrices of the elements of LinIneq
 * @param tpm The input TPM.
 */
void LinIneq::fill(const TPM &tpm){

   tr = tpm.trace();

   for(int i = 0;i < nr;++i)
      proj[i] = (li[i]->gI()).ddot(tpm) + li[i]->gI_tr() * tr;

   SPM spm;
   spm.bar(1.0,tpm);

   for(int i = 0;i < nr;++i)
      proj_bar[i] = (li[i]->gI_bar()).ddot(spm) + li[i]->gI_tr() * (M - 1.0) * 2.0 * tr;

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
 * @param i the index of the constraint we are interested in.
 * @return the barred projection on the barred constaint with index i
 */
double LinIneq::gproj_bar(int i) const {

   return proj_bar[i];

}

/**
 * @return the array containing the barred projection on the barred constraints.
 */
double *LinIneq::gproj_bar(){

   return proj_bar;

}

/**
 * @return the trace of the input TPM
 */
double LinIneq::gtr() const{

   return tr;

}

/**
 * @return the value of the "a" coefficient of the overlapmatrix-map
 */ 
double LinIneq::ga() const {

   return a;

}

/**
 * @return the value of the "c" coefficient of the overlapmatrix-map
 */ 
double LinIneq::gc() const {

   return c;

}

/**
 * The "alpha" function, see notes, projects a LinIneq (actually a TPM) onto a scalar.
 * @param index the index of the alpha function
 * @return the alpha function value for this LinIneq object
 */
double LinIneq::alpha(int index) const {

   double tmp = 0.0;

   int n = 2*nr;

   //alpha_1^i
   for(int k = 0;k < nr;++k)
      tmp += 4.0*coef[k*n + index] * proj[k];

   //alpha_2^i
   for(int k = nr;k < 2*nr;++k)
      tmp += coef[k*n + index] * proj_bar[k - nr];

   return tmp;

}

/**
 * @return inproduct of (*this) LinIneq with li_i, defined as (A^T B)
 * @param li_i input matrix
 */
double LinIneq::ddot(const LinIneq &li_i) const {

   int inc = 1;

   return ddot_(&nr,proj,&inc,li_i.proj,&inc);

}

ostream &operator<<(ostream &output,const LinIneq &li_p){

   for(int i = 0;i < LinIneq::gnr();++i)
      output << i << "\t" << li_p.gproj(i) << endl;

   return output;

}

/**
 * Multiply LinIneq object left and right with another LinIneq object
 * @param map LinIneq object that will be multiplied to the left en to the right of the matrix object
 * @param object central LinIneq object
 */
void LinIneq::L_map(const LinIneq &map,const LinIneq &object){

   for(int i = 0;i < nr;++i)
      proj[i] = map.gproj(i)*object.gproj(i)*map.gproj(i);

}

/**
 * Take the square root out of the projections onto the Linear Constraints
 * @param option = 1, positive square root, = -1, inverse square root.
 */
void LinIneq::sqrt(int option){

   if(option == 1){

      for(int i = 0;i < nr;++i)
         proj[i] = std::sqrt(proj[i]);

   }
   else{

      for(int i = 0;i < nr;++i)
         proj[i] = 1.0/std::sqrt(proj[i]);

   }

}

/**
 * Multiply the projections with a constant factor  alpha
 * @param alpha the scaling factor
 */
void LinIneq::dscal(double alpha){

   int inc = 1;

   dscal_(&nr,&alpha,proj,&inc);

}

/**
 * Invert the projections
 */
void LinIneq::invert(){

   for(int i = 0;i < nr;++i)
      proj[i] = 1.0/proj[i];

}

/**
 * Fill the projection randomly
 */
void LinIneq::fill_Random(){

   for(int i = 0;i < nr;++i)
      proj[i] = (double) rand()/RAND_MAX;

}

/**
 * Test function that prints the coefficients of the inverse linear system of equations for the overlapmatrix.
 */
void LinIneq::print_coef() {

   for(int i = 0;i < 4*nr*nr;++i)
      cout << i << "\t" << coef[i] << endl;

}

/**
 * seperate into positive and negative part.
 * @param p positive part
 * @param m negative part
 */
void LinIneq::sep_pm(LinIneq &p,LinIneq &m){

   for(int i = 0;i < nr;++i){

      if(proj[i] > 0.0){

         p.proj[i] = proj[i];
         m.proj[i] = 0.0;

      }
      else{

         m.proj[i] = proj[i];
         p.proj[i] = 0.0;

      }

   }

}