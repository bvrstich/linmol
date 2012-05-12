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

int LinCon::M;
int LinCon::N;

/**
 * initialize the static lists and variables
 * @param M_in input dimension of sp space
 * @param N_in input nr of particles
 */
void LinCon::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * Constructor of a LinCon object
 * @param M The constraint matrix
 * @param N the minimal projection
 */
LinCon::LinCon(){

   I_c = new TPM();

   I_c_bar = new SPM();

}

/**
 * copy constructor
 * @param lc_copy The LinCon object to be copied
 */
LinCon::LinCon(const LinCon &lc_copy){

   I_c = new TPM(lc_copy.gI());

   I_c_bar = new SPM(lc_copy.gI_bar());

   i_c = lc_copy.gi();

   this->I_c_tr = lc_copy.gI_tr();

}

/**
 * destructor
 */
LinCon::~LinCon(){

   delete I_c;

   delete I_c_bar;

}

/**
 * @return the constraint TPM object
 */
const TPM &LinCon::gI() const{

   return *I_c;

}

/**
 * @return the partially trace constraint, the SPM object I_c_bar.
 */
const SPM &LinCon::gI_bar() const{

   return *I_c_bar;

}

/**
 * @return the scaled trace of the constraint: 2*Tr(I_c)/M(M-1)
 */
double LinCon::gI_tr() const{

   return I_c_tr;

}

/**
 * @return The minimal projection
 */
double LinCon::gi() const{

   return i_c;

}

/**
 * set the constraint value
 * @param i the value that the minimal projection will be set to.
 */
void LinCon::si(double i){

   i_c = i;

}

/**
 * set the constraint Matrix, warning, first set the value i_c, because otherwise the shift will be wrong
 * @param I the input constraint Matrix
 */
void LinCon::sI(const TPM &I){

   *I_c = I;

   for(int B = 0;B < I.gnr();++B)
      for(int i = 0;i < I.gdim(B);++i)//shift it for convenience
         (*I_c)(B,i,i) -= 2.0*i_c/(N*(N - 1.0));

   I_c_tr = 2.0 * I_c->trace()/(double) (M*(M-1.0) );

   //project onto traceless matrix space
   I_c->proj_Tr();

   //make the bar
   I_c_bar->bar(1.0,I);

}

ostream &operator<<(ostream &output,const LinCon &lc_p){

   output << "The minimal projection:\t" << lc_p.gi() << endl;
   output << endl;

   output << "The shifted, traceless, Constraint matrix:" << endl;
   output << endl;

   output << lc_p.gI() << endl;

   return output;

}

/**
 * @return nr of sp orbitals
 */
int LinCon::gM() const{

   return M;

}

/**
 * @return nr of particles
 */
int LinCon::gN() const{

   return N;

}

/**
 * fills the constraints object randomly
 */
void LinCon::fill_Random(){

   i_c = (double) rand()/RAND_MAX;

   I_c->fill_Random();

   for(int B = 0;B < I_c->gnr();++B)
      for(int i = 0;i < I_c->gdim(B);++i)
         (*I_c)(B,i,i) -= 2.0*i_c/(N*(N - 1.0));

   I_c_tr = 2.0 * I_c->trace()/(double) (M*(M-1.0));

   I_c->proj_Tr();

   I_c_bar->bar(1.0,*I_c);

}

/**
 * construct the spin matrix as the spin matrix
 * @param spin the value of the spinconstraint
 */
void LinCon::spincon(double spin){

   I_c->set_S_2();

   for(int B = 0;B < I_c->gnr();++B)
      for(int i = 0;i < I_c->gdim(B);++i)
         (*I_c)(B,i,i) -= 2.0*spin/(N*(N - 1.0));

   I_c_tr = 2.0 * I_c->trace()/ (double)(M*(M - 1.0));

   I_c->proj_Tr();

   I_c_bar->bar(1.0,*I_c);

   i_c = spin;

}
