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
   *Q -= sup_c.gQ();
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
SUP &SUP::operator=(double &a){

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
 * Initialization of the SUP matrix S, is just u^0: see primal_dual.pdf for more information
 */
void SUP::init_S(){

   I->unit();

   this->fill();

}

/**
 * Initialisation for dual SUP matrix Z, see primal_dual.pdf for info.
 */
void SUP::init_Z(double alpha,const TPM &ham,const SUP &u_0) {

   this->fill_Random();

   this->proj_C(ham);

   //nog een eenheidsmatrix maal constante bijtellen zodat Z positief definiet is:
   this->daxpy(alpha,u_0); 

}

/**
 * Project the general SUP matrix (*this) orthogonally onto the linear space for which\n\n
 * Tr(Z u^i) = h^i      with h^i = Tr(tpm f^i)\n\n
 * is valid.
 * @param tpm input TPM (mostly the hamiltonian of the problem)
 */
void SUP::proj_C(const TPM &tpm) {

   TPM hulp;

   hulp.collaps(0,*this);

   hulp -= tpm;

   hulp.proj_Tr();

   //Z_res is the orthogonal piece of this that will be deducted,
   //so the piece of this in the U-space - ham
   SUP Z_res;

   //apply iverse S to it and put it in Z_res.tpm(0)
   Z_res.gI().S(-1,hulp);

   //and fill it up Johnny
   Z_res.fill();

   *this -= Z_res;

}


/**
 * @return Deviation from the central path measured through the logarithmic potential, it's a measure for
 * the deviation of the product of the primal with the dual matrix (SZ) from the unit matrix.\n
 * Usage of the function: S.center_dev(Z) gives returns the deviation.\n\n
 * (*this) = S = primal matrix of the problem
 * @param Z = dual matrix of the problem
 */
double SUP::center_dev(const SUP &Z) const {

   SUP sqrt_S(*this);

   sqrt_S.sqrt(1);

   SUP SZ;
   SZ.L_map(sqrt_S,Z);

   EIG eig(SZ);

   return eig.center_dev();

}

/**
 * Construct the D matrix and put it in this, D is the matrix matrix of the hessian, see primal_dual.pdf for more information
 * @param S The primal SUP matrix S
 * @param Z The dual SUP matrix Z
 */
void SUP::D(const SUP &S,const SUP &Z) {

   //positieve vierkantswortel uit Z
   SUP Z_copy(Z);

   Z_copy.sqrt(1);

   //links en rechts vermenigvuldigen met wortel Z
   SUP hulp;

   hulp.L_map(Z_copy,S);

   hulp.sqrt(1);

   //negatieve vierkantswortel uit Z
   Z_copy = Z;

   Z_copy.sqrt(-1);

   //en links en rechts hulp vermenigvuldigen met die wortel, en in this steken:
   this->L_map(Z_copy,hulp);

}

/**
 * Orthogonal projection of a general SUP matrix [ M M_Q ( M_G M_T1 M_T2 ) ] onto the orthogonal complement of the U space (C space)
 * See primal_dual.pdf for more information
 */
void SUP::proj_C(){

   SUP Z_copy(*this);

   //projecteer op de U ruimte
   Z_copy.proj_U();

   //en het orthogonaal complement nemen:
   *this -= Z_copy;

}

/**
 * Orthogonal projection of a general SUP matrix diag[ M M_Q ( M_G M_T1 M_T2 ) ] onto U space: diag[ M_u Q(M_u) ( G(M_u) T1(M_u) T2(M_u) ) ]
 * for more information see primal_dual.pdf
 */
void SUP::proj_U(){
  
   //eerst M_Gamma + Q(M_Q) + ( G(M_G) + T1(M_T1) + T2(M_T2) ) in O stoppen
   TPM O;

   O.collaps(1,*this);

   //dan de inverse overlapmatrix hierop laten inwerken en in this[0] stoppen
   I->S(-1,O);

   //fill up the rest with the right maps
   this->fill();

}

/**
 * Implementation of the linear conjugate gradient algorithm for the solution of the dual Newton system\n\n
 * H(*this) = B in which H is the dual hessian map
 * @param B right hand side of the equation
 * @param D SUP matrix that defines the structure of the hessian map (the metric) (inverse of the primal Newton equation hessian)
 * @return return the number of iteration required to converge
 */
int SUP::solve(SUP &B,const SUP &D) {

   SUP HB;
   HB.H(*this,D);

   B -= HB;

   //de r initialiseren op B - H DZ
   SUP r(B);

   double rr = r.ddot(r);
   double rr_old,ward;

   int cg_iter = 0;

   while(rr > 1.0e-5){

      ++cg_iter;

      HB.H(B,D);

      ward = rr/B.ddot(HB);

      //delta += ward*b
      this->daxpy(ward,B);

      //r -= ward*HB
      r.daxpy(-ward,HB);

      //nieuwe r_norm maken
      rr_old = rr;
      rr = r.ddot(r);

      //eerst herschalen van b
      B.dscal(rr/rr_old);

      //dan r er bijtellen
      B += r;

   }
   
   return cg_iter;

}

/**
 * The dual hessian map:\n\n
 * HB = DBD (dus SUP::L_map), projected onto C-space (SUP::proj_C)
 * @param B SUP matrix onto which the hessian works.
 * @param D SUP matrix that defines the structure of the map (metric)
 */
void SUP::H(const SUP &B,const SUP &D) {

   this->L_map(D,B);

   this->proj_C();

}

/**
 * Line search function that checks how large a step you can take in a given primal dual predictor direction (DS,DZ), starting from 
 * the current primal dual point (S,Z), before deviating beyond max_dev from the central path.\n\n
 * (*this) = DS --> primal search direction
 * @param DZ dual search direction
 * @param S Current primal point
 * @param Z Current dual point
 * @param max_dev number (double) input by which you can tell the function how far you want to deviate from the central path after the step.
 */
double SUP::line_search(const SUP &DZ,const SUP &S,const SUP &Z,double max_dev) const
{
   //eerst de huidige deviatie van het centraal pad nemen:
   double center_dev = S.center_dev(Z);

   //eigenwaarden zoeken van S^{-1/2} DS S^{-1/2} en Z^{-1/2} DZ Z^{-1/2}

   //kopieer S in de zogeheten wortel:
   SUP wortel(S);

   //maak negatieve vierkantswortel uit S
   wortel.sqrt(-1);

   //de L_map
   SUP hulp;
   hulp.L_map(wortel,*this);

   //eigenwaarden in eigen_S stoppen
   EIG eigen_S(hulp);

   //nu idem voor Z
   wortel = Z;

   wortel.sqrt(-1);

   hulp.L_map(wortel,DZ);

   EIG eigen_Z(hulp);

   //nog c_S en c_Z uitrekenen:
   double pd_gap = S.ddot(Z);

   //c_S = Tr (DS Z)/Tr (SZ)
   double c_S = this->ddot(Z)/pd_gap;

   //c_Z = Tr (S DZ)/Tr (SZ)
   double c_Z = S.ddot(DZ)/pd_gap;

   //waar zitten de singulariteiten: tot waar mag ik zoeken?
   double a_max = -1.0/eigen_S.min();
   double b_max = -1.0/eigen_Z.min();

   //a_max is de waarde tot waar ik zal zoeken:
   if(b_max < a_max)
      a_max = b_max;

   double a = 0.0;
   double b = a_max;

   double c = (a + b)/2.0;

   //bissectiemethode om stapgrootte te bepalen:
   while(b - a > 1.0e-5){

      c = (a + b)/2.0;

      if( (center_dev + eigen_S.centerpot(c,eigen_Z,c_S,c_Z) - max_dev) < 0.0 )
         a = c;
      else
         b = c;

   }

   return c;

}
