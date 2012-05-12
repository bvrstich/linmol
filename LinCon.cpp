#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

int LinCon::N;
int LinCon::M;

/**
 * initialize the static lists and variables
 * @param M_in input dimension of sp space
 * @param N_in input nr of particles
 */
void DPM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * Constructor of a LinCon object
 */
LinCon::LinCon(){

   I_c = new TPM();

}

/**
 * copy constructor
 * @param lc_copy The LinCon object to be copied
 */
LinCon::LinCon(const LinCon &lc_copy){

   I_c = new TPM(lc_copy.gI());

   i_c = lc_copy.gi();

}

/**
 * destructor
 */
LinCon::~LinCon(){

   delete I_c;

}

/**
 * @return the constraint TPM object
 */
TPM &LinCon::gI() const{

   return *I_c;

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
 * set the constraint Matrix
 * @param I the input constraint Matrix
 */
void LinCon::sI(const TPM &I){

   *I_c = I;

}

ostream &operator<<(ostream &output,const LinCon &lc_p){

   cout << "The minimal projection:\t" << lc_p.gi() << endl;
   cout << endl;

   cout << "The Constraint matrix:" << endl;
   cout << endl;

   cout << lc_p.gI() << endl;

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
 * construct the spin matrix as the spin matrix
 */
void LinCon::spincon(){

   I_c->set_S_2();

}
