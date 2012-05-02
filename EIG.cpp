#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

int EIG::M;
int EIG::N;
int EIG::dim;

/**
 * initialize the statics
 * @param M_in the nr of sites
 * @param N_in the nr of particles
 */
void EIG::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

   dim = M*(M - 1)/2;

#ifdef __Q_CON
   dim += M*(M - 1)/2;
#endif

#ifdef __G_CON
   dim += M*M;
#endif

#ifdef __T1_CON
   dim += M*(M - 1)*(M - 2)/6;
#endif

#ifdef __T2_CON
   dim += M*M*(M - 1)/2;
#endif

}

/**
 * standard constructor with initialization on the eigenvalues of a SUP object.
 * @param sup input SUP object that will be destroyed after this function is called. The eigenvectors
 * of the matrix will be stored in the columns of the original SUP matrix.
 */
EIG::EIG(SUP &sup){

   v_I = new BlockVector<TPM>(sup.gI());

#ifdef __Q_CON
   v_Q = new BlockVector<TPM>(sup.gQ());
#endif

#ifdef __G_CON
   v_G = new BlockVector<PHM>(sup.gG());
#endif

#ifdef __T1_CON
   v_T1 = new BlockVector<DPM>(sup.gT1());
#endif

#ifdef __T2_CON
   v_T2 = new BlockVector<PPHM>(sup.gT2());
#endif

}

/**
 * Copy constructor\n
 * allocates the memory for the eigenvalues of a SUP object and copies the content of eig_c into it.
 * @param eig_c The input EIG that will be copied into this.
 */
EIG::EIG(const EIG &eig_c){

   v_I = new BlockVector<TPM>(eig_c.gv_I());

#ifdef __Q_CON
   v_Q = new BlockVector<TPM>(eig_c.gv_Q());
#endif

#ifdef __G_CON
   v_G = new BlockVector<PHM>(eig_c.gv_G());
#endif

#ifdef __T1_CON
   v_T1 = new BlockVector<DPM>(eig_c.gv_T1());
#endif

#ifdef __T2_CON
   v_T2 = new BlockVector<PPHM>(eig_c.gv_T2());
#endif

}

/**
 * overload equality operator
 * @param eig_c object that will be copied into this.
 */
EIG &EIG::operator=(const EIG &eig_c){

   *v_I = eig_c.gv_I();

#ifdef __Q_CON
   *v_Q = eig_c.gv_Q();
#endif

#ifdef __G_CON
   *v_G = eig_c.gv_G();
#endif

#ifdef __T1_CON
   *v_T1 = eig_c.gv_T1();
#endif

#ifdef __T2_CON
   *v_T2 = eig_c.gv_T2();
#endif

   return *this;

}

/**
 * Destructor, deallocation of the memory
 */
EIG::~EIG(){

   delete v_I;

#ifdef __Q_CON
   delete v_Q;
#endif

#ifdef __G_CON
   delete v_G;
#endif

#ifdef __T1_CON
   delete v_T1;
#endif

#ifdef __T2_CON
   delete v_T2;
#endif

}

/**
 * Diagonalize a SUP matrix when the memory has allready been allocated before
 * @param sup matrix to be diagonalized
 */
void EIG::diagonalize(SUP &sup){

   v_I->diagonalize(sup.gI());

#ifdef __Q_CON
   v_Q->diagonalize(sup.gQ());
#endif

#ifdef __G_CON
   v_G->diagonalize(sup.gG());
#endif

#ifdef __T1_CON
   v_T1->diagonalize(sup.gT1());
#endif

#ifdef __T2_CON
   v_T2->diagonalize(sup.gT2());
#endif

}

ostream &operator<<(ostream &output,const EIG &eig_p){

   std::cout << eig_p.gv_I() << std::endl;

#ifdef __Q_CON
   std::cout << eig_p.gv_Q() << std::endl;
#endif

#ifdef __G_CON
   std::cout << eig_p.gv_G() << std::endl;
#endif

#ifdef __T1_CON
   std::cout << eig_p.gv_T1() << std::endl;
#endif

#ifdef __T2_CON
   std::cout << eig_p.gv_T2() << std::endl;
#endif

   return output;

}

/**
 * @return nr of particles
 */
int EIG::gN() const{

   return N;

}

/**
 * @return dimension of sp space
 */
int EIG::gM() const{

   return M;

}

/** 
 * get the BlockVector<TPM> object containing the eigenvalues of the TPM block I
 * @return a BlockVector<TPM> object containing the desired eigenvalues
 */
BlockVector<TPM> &EIG::gv_I(){

   return *v_I;

}

/** 
 * const version\n\n
 * get the BlockVector<TPM> object containing the eigenvalues of the TPM block I
 * @return a BlockVector<TPM> object containing the desired eigenvalues
 */
const BlockVector<TPM> &EIG::gv_I() const{

   return *v_I;

}

#ifdef __Q_CON

/** 
 * get the BlockVector<TPM> object containing the eigenvalues of the TPM block Q
 * @return a BlockVector<TPM> object containing the desired eigenvalues
 */
BlockVector<TPM> &EIG::gv_Q(){

   return *v_Q;

}

/** 
 * const version\n\n
 * get the BlockVector<TPM> object containing the eigenvalues of the TPM block Q
 * @return a BlockVector<TPM> object containing the desired eigenvalues
 */
const BlockVector<TPM> &EIG::gv_Q() const{

   return *v_Q;

}

#endif

#ifdef __G_CON

/** 
 * get the BlockVector<PHM> object containing the eigenvalues of the PHM block G
 * @return a BlockVector<PHM> object containing the desired eigenvalues
 */
BlockVector<PHM> &EIG::gv_G(){

   return *v_G;

}

/** 
 * get the BlockVector<PHM> object containing the eigenvalues of the PHM block G, const version
 * @return a BlockVector<PHM> object containing the desired eigenvalues
 */
const BlockVector<PHM> &EIG::gv_G() const{

   return *v_G;

}

#endif

#ifdef __T1_CON

/** 
 * get the BlockVector<DPM> object containing the eigenvalues of the DPM block T1 of the SUP matrix
 * @return a BlockVector<DPM> object containing the desired eigenvalues
 */
BlockVector<DPM> &EIG::gv_T1(){

   return *v_T1;

}

/** 
 * get the BlockVector<DPM> object containing the eigenvalues of the DPM block T1 of the SUP matrix, const version
 * @return a BlockVector<DPM> object containing the desired eigenvalues
 */
const BlockVector<DPM> &EIG::gv_T1() const{

   return *v_T1;

}

#endif

#ifdef __T2_CON

/** 
 * get the BlockVector<PPHM> object containing the eigenvalues of the PPHM block T2 of the SUP matrix
 * @return a BlockVector<PPHM> object containing the desired eigenvalues
 */
BlockVector<PPHM> &EIG::gv_T2(){

   return *v_T2;

}

/** 
 * get the BlockVector<PPHM> object containing the eigenvalues of the PPHM block T2 of the SUP matrix, const version
 * @return a BlockVector<PPHM> object containing the desired eigenvalues
 */
const BlockVector<PPHM> &EIG::gv_T2() const{

   return *v_T2;

}

#endif

/**
 * @return total dimension of the EIG object
 */
int EIG::gdim() const{

   return dim;

}


/**
 * @return the minimal element present in this EIG object.
 * watch out, only works when EIG is filled with the eigenvalues of a diagonalized SUP matrix
 */
double EIG::min() const{

   //lowest eigenvalue of P block
   double ward = v_I->min();

#ifdef __Q_CON
   //lowest eigenvalue of Q block
   if(ward > v_Q->min())
      ward = v_Q->min();
#endif

#ifdef __G_CON
   //lowest eigenvalue of G block
   if(ward > v_G->min())
      ward = v_G->min();
#endif

#ifdef __T1_CON
   //lowest eigenvalue of the T1 block
   if(ward > v_T1->min())
      ward = v_T1->min();
#endif

#ifdef __T2_CON
   //lowest eigenvalue of the T2 block
   if(ward > v_T2->min())
      ward = v_T2->min();
#endif

   return ward;

}

/**
 * @return the maximum element present in this EIG object.
 * watch out, only works when EIG is filled with the eigenvalues of a diagonalized SUP matrix
 */
double EIG::max() const{

   //highest eigenvalue of P block
   double ward = v_I->max();

#ifdef __Q_CON
   //highest eigenvalue of Q block
   if(ward < v_Q->max())
      ward = v_Q->max();
#endif

#ifdef __G_CON
   //highest eigenvalue of G block
   if(ward < v_G->max())
      ward = v_G->max();
#endif

#ifdef __T1_CON
   //highest eigenvalue of the T1 block
   if(ward < v_T1->max())
      ward = v_T1->max();
#endif

#ifdef __T2_CON
   //highest eigenvalue of the T2 block
   if(ward < v_T2->max())
      ward = v_T2->max();
#endif

   return ward;

}
