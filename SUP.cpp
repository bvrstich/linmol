#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

int SUP::M;
int SUP::N;
int SUP::dim;

/**
 * initialize the statics
 * @param M_in the dimension of sp Hilbert space
 * @param N_in the nr of particles
 */
void SUP::init(int M_in,int N_in){

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
 * standard constructor\n
 * Allocates two TPM matrices and optionally a PHM, DPM or PPHM matrix.
 */
SUP::SUP(){

   I = new TPM();

#ifdef __Q_CON
   Q = new TPM();
#endif

#ifdef __G_CON
   G = new PHM();
#endif

#ifdef __T1_CON
   T1 = new DPM();
#endif

#ifdef __T2_CON
   T2 = new PPHM();
#endif

}

/**
 * standard constructor\n
 * Allocates two TPM matrices and optionally a PHM, DPM or PPHM matrix, then copies the content of
 * input SUP sup_c into it.
 * @param sup_c input SUP
 */
SUP::SUP(const SUP &sup_c){

   I = new TPM(sup_c.gI());

#ifdef __Q_CON
   Q = new TPM(sup_c.gQ());
#endif

#ifdef __G_CON
   G = new PHM(sup_c.gG());
#endif

#ifdef __T1_CON
   T1 = new DPM(sup_c.gT1());
#endif

#ifdef __T2_CON
   T2 = new PPHM(sup_c.gT2());
#endif

}

/**
 * Destructor
 */
SUP::~SUP(){

   delete I;

#ifdef __Q_CON
   delete Q;
#endif

#ifdef __G_CON
   delete G;
#endif

#ifdef __T1_CON
   delete T1;
#endif

#ifdef __T2_CON
   delete T2;
#endif

}

/**
 * Overload += operator
 * @param sup_p The SUP matrix that has to be added to this
 */
SUP &SUP::operator+=(const SUP &sup_p){

   *I += sup_p.gI();

#ifdef __Q_CON
   *Q += sup_p.gQ();
#endif

#ifdef __G_CON
   *G += sup_p.gG();
#endif

#ifdef __T1_CON
   *T1 += sup_p.gT1();
#endif

#ifdef __T2_CON
   *T2 += sup_p.gT2();
#endif

   return *this;

}

/**
 * Overload -= operator
 * @param sup_p The SUP that will be deducted from this
 */
SUP &SUP::operator-=(const SUP &sup_p){

   *I -= sup_p.gI();

#ifdef __Q_CON
   *Q -= sup_p.gQ();
#endif

#ifdef __G_CON
   *G -= sup_p.gG();
#endif

#ifdef __T1_CON
   *T1 -= sup_p.gT1();
#endif

#ifdef __T2_CON
   *T2 -= sup_p.gT2();
#endif

   return *this;

}

/**
 * Overload equality operator, copy sup_c into this
 * @param sup_c SUP object to be copied into this
 */
SUP &SUP::operator=(const SUP &sup_c){

   *I = sup_c.gI();

#ifdef __Q_CON
   *Q = sup_c.gQ();
#endif

#ifdef __G_CON
   (*G) = sup_c.gG();
#endif

#ifdef __T1_CON
   (*T1) = sup_c.gT1();
#endif

#ifdef __T2_CON
   (*T2) = sup_c.gT2();
#endif

   return *this;

}

/**
 * overload operator = number, all the blockmatrices in SUP are put equal to the number a.
 * e.g. SZ = 0 makes all the Matrix elements zero.
 * @param a the number
 */
SUP &SUP::operator=(double a){

   *I = a;

#ifdef __Q_CON
   *Q = a;
#endif

#ifdef __G_CON
   (*G) = a;
#endif

#ifdef __T1_CON
   (*T1) = a;
#endif

#ifdef __T2_CON
   (*T2) = a;
#endif

   return *this;

}

/**
 * @return pointer to the I object
 */
TPM &SUP::gI(){

   return *I;

}

/**
 * const version
 * @return pointer to the I object
 */
const TPM &SUP::gI() const{

   return *I;

}

#ifdef __Q_CON

/**
 * @return pointer to the individual Q object
 */
TPM &SUP::gQ(){

   return *Q;

}

/**
 * const version
 * @return pointer to the Q object, const version
 */
const TPM &SUP::gQ() const{

   return *Q;

}

#endif

#ifdef __G_CON

/**
 * @return pointer to the PHM block: G
 */
PHM &SUP::gG(){

   return *G;

}

/**
 * const version
 * @return pointer to the PHM block: G
 */
const PHM &SUP::gG() const{

   return *G;

}

#endif

#ifdef __T1_CON

/**
 * @return pointer to the DPM block: T1
 */
DPM &SUP::gT1(){

   return *T1;

}

/**
 * const version
 * @return pointer to the DPM block: T1
 */
const DPM &SUP::gT1() const{

   return *T1;

}

#endif

#ifdef __T2_CON

/**
 * const version
 * @return pointer to the PPHM block: T2
 */
const PPHM &SUP::gT2() const{

   return *T2;

}

/**
 * @return pointer to the PPHM block: T2
 */
PPHM &SUP::gT2(){

   return *T2;

}

#endif

/**
 * Fill the SUP matrix with random elements, watch out, not necessarily positive definite
 */
void SUP::fill_Random(){

   I->fill_Random();

#ifdef __Q_CON
   Q->fill_Random();
#endif

#ifdef __G_CON
   G->fill_Random();
#endif

#ifdef __T1_CON
   T1->fill_Random();
#endif

#ifdef __T2_CON
   T2->fill_Random();
#endif

}

/**
 * @return number of particles
 */
int SUP::gN() const {

   return N;

}

/**
 * @return dimension of sp space
 */
int SUP::gM() const{

   return M;

}

/**
 * @return total dimension of SUP (carrier) space
 */
int SUP::gdim() const{

   return dim;

}

/**
 * @param sup_i input SUP object
 * @return inproduct between this and input matrix sup_i, defined as Tr(this sup_i)
 */
double SUP::ddot(const SUP &sup_i) const{

   double ward = I->ddot(sup_i.gI());

#ifdef __Q_CON
   ward += Q->ddot(sup_i.gQ());
#endif

#ifdef __G_CON
   ward += G->ddot(sup_i.gG());
#endif

#ifdef __T1_CON
   ward += T1->ddot(sup_i.gT1());
#endif

#ifdef __T2_CON
   ward += T2->ddot(sup_i.gT2());
#endif

   return ward;

}

/**
 * Invert all the Matrices in SUP and put it in this, watch out, destroys original matrices.
 * Makes use of cholesky decomposition, so only positive definite matrices can be used as input!
 */
void SUP::invert(){

   I->invert();

#ifdef __Q_CON
   Q->invert();
#endif

#ifdef __G_CON
   G->invert();
#endif

#ifdef __T1_CON
   T1->invert();
#endif

#ifdef __T2_CON
   T2->invert();
#endif

}

/**
 * Scale all the elements in the blockmatrix with parameter alpha
 * @param alpha the scalefactor
 */
void SUP::dscal(double alpha){

   I->dscal(alpha);

#ifdef __Q_CON
   Q->dscal(alpha);
#endif

#ifdef __G_CON
   G->dscal(alpha);
#endif

#ifdef __T1_CON
   T1->dscal(alpha);
#endif

#ifdef __T2_CON
   T2->dscal(alpha);
#endif

}

/**
 * Take the square root out of the positive semidefinite SUP matrix this and put it in this, watch out, original matrix is destroyed.
 * @param option = 1 positive square root, = -1 negative square root.
 */
void SUP::sqrt(int option){

   I->sqrt(option);

#ifdef __Q_CON
   Q->sqrt(option);
#endif

#ifdef __G_CON
   G->sqrt(option);
#endif

#ifdef __T1_CON
   T1->sqrt(option);
#endif

#ifdef __T2_CON
   T2->sqrt(option);
#endif

}

/**
 * Multiply symmetric SUP blockmatrix object left en right with symmetric SUP blockmatrix map to 
 * form another symmetric matrix and put it in (*this): this = map*object*map
 * @param map SUP that will be multiplied to the left en to the right of matrix object
 * @param object central SUP
 */
void SUP::L_map(const SUP &map,const SUP &object){

   I->L_map(map.gI(),object.gI());

#ifdef __Q_CON
   Q->L_map(map.gQ(),object.gQ());
#endif

#ifdef __G_CON
   G->L_map(map.gG(),object.gG());
#endif

#ifdef __T1_CON
   T1->L_map(map.gT1(),object.gT1());
#endif

#ifdef __T2_CON
   T2->L_map(map.gT2(),object.gT2());
#endif

}

/**
 * add the SUP sup_p times the constant alpha to this
 * @param alpha the constant to multiply the sup_p with
 * @param sup_p the SUP to be multiplied by alpha and added to (*this)
 */
void SUP::daxpy(double alpha,const SUP &sup_p){

   I->daxpy(alpha,sup_p.gI());

#ifdef __Q_CON
   Q->daxpy(alpha,sup_p.gQ());
#endif

#ifdef __G_CON
   G->daxpy(alpha,sup_p.gG());
#endif

#ifdef __T1_CON
   T1->daxpy(alpha,sup_p.gT1());
#endif

#ifdef __T2_CON
   T2->daxpy(alpha,sup_p.gT2());
#endif

}

/**
 * General matrixproduct between two SUP matrices, act with Matrix::mprod on every block
 * 
 * @param A left hand matrix
 * @param B right hand matrix
 * @return The product AB
 */
SUP &SUP::mprod(const SUP &A,const SUP &B){

   I->mprod(A.gI(),B.gI());

#ifdef __Q_CON
   Q->mprod(A.gQ(),B.gQ());
#endif

#ifdef __G_CON
   G->mprod(A.gG(),B.gG());
#endif

#ifdef __T1_CON
   T1->mprod(A.gT1(),B.gT1());
#endif

#ifdef __T2_CON
   T2->mprod(A.gT2(),B.gT2());
#endif

   return *this;

}

/**
 * Fill the SUP matrix (*this) with a TPM matrix like: this = diag[tpm  Q(tpm)  ( G(tpm) T1(tpm) T2(tpm) ) ]
 * @param tpm input TPM
 */
void SUP::fill(const TPM &tpm){

   *I = tpm;

#ifdef __Q_CON
   Q->Q(1,tpm);
#endif

#ifdef __G_CON
   G->G(tpm);
#endif

#ifdef __T1_CON
   T1->T(tpm);
#endif

#ifdef __T2_CON
   T2->T(tpm);
#endif

}

/**
 * fill the SUP matrix with the TPM matrix stored in the first block:\n\n
 */
void SUP::fill(){

#ifdef __Q_CON
   Q->Q(1,*I);
#endif 

#ifdef __G_CON
   G->G(*I);
#endif 

#ifdef __T1_CON
   T1->T(*I);
#endif 

#ifdef __T2_CON
   T2->T(*I);
#endif 

}

ostream &operator<<(ostream &output,const SUP &sup_p){

   output << sup_p.gI() << std::endl;

#ifdef __Q_CON
   output << std::endl;
   output << sup_p.gQ();
#endif

#ifdef __G_CON
   output << std::endl;
   output << sup_p.gG();
#endif

#ifdef __T1_CON
   output << std::endl;
   output << sup_p.gT1();
#endif

#ifdef __T2_CON
   output << std::endl;
   output << sup_p.gT2();
#endif

   return output;

}

/**
 * Seperate SUP into two SUP's, a positive and negative semidefinite part.
 * @param p positive (plus) output part
 * @param m negative (minus) output part
 */
void SUP::sep_pm(SUP &p,SUP &m){

   I->sep_pm(p.gI(),m.gI());

#ifdef __Q_CON
   Q->sep_pm(p.gQ(),m.gQ());
#endif

#ifdef __G_CON
   G->sep_pm(p.gG(),m.gG());
#endif

#ifdef __T1_CON
   T1->sep_pm(p.gT1(),m.gT1());
#endif

#ifdef __T2_CON
   T2->sep_pm(p.gT2(),m.gT2());
#endif

}

/* vim: set ts=3 sw=3 expandtab :*/
