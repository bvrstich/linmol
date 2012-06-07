#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::complex;

#include "include.h"

vector< vector<int> > SphInt::s2inlm;
int ****SphInt::inlm2s;

vector< vector<int> > SphInt::t2s;
int **SphInt::s2t;

int SphInt::dim;
int SphInt::N;
int SphInt::N_Z;
int SphInt::n_max;
int SphInt::l_max;

double SphInt::NucRepEn;

/** 
 * static function that allocates the static lists and calculates the dimensions and such
 */
void SphInt::init(){

   N_Z = CartInt::gN_Z();
   n_max = CartInt::gn_max();
   l_max = CartInt::gl_max();
   N = CartInt::gN();
   NucRepEn = CartInt::gNucRepEn();

   //allocate
   inlm2s = new int *** [N_Z];

   for(int i = 0;i < N_Z;++i){

      inlm2s[i] = new int ** [n_max];

      for(int n = 0;n < n_max;++n){

         inlm2s[i][n] = new int * [l_max + 1];

         for(int l = 0;l <= l_max;++l)
            inlm2s[i][n][l] = new int [2*l + 1];

      }
   }

   //construct
   vector<int> v(4);

   for(int s = 0;s < CartInt::gdim();++s){

      v[0] = CartInt::gs2inlxyz(s,0);//i
      v[1] = CartInt::gs2inlxyz(s,1);//n
      v[2] = CartInt::gs2inlxyz(s,2);//l

      for(int m = -v[2];m <= v[2];++m){

         v[3] = m;//m

         s2inlm.push_back(v);

      }

      s += (v[2] + 2)*(v[2] + 1)/2 - 1;

   }

   dim = s2inlm.size();

   for(int s = 0;s < dim;++s){

      v = s2inlm[s];

      inlm2s[v[0]][v[1] - v[2] - 1][v[2]][v[3] + v[2]] = s;

   }

   s2t = new int * [dim];

   for(int i = 0;i < dim;++i)
      s2t[i] = new int [dim];

   vector<int> vst(2);

   int iter = 0;

   for(int i = 0;i < dim;++i)
      for(int j = 0;j < dim;++j){

         vst[0] = i;
         vst[1] = j;

         t2s.push_back(vst);

         s2t[i][j] = iter;

         ++iter;

      }

}

/** 
 * function that deallocates the static members
 */
void SphInt::clear(){

   for(int i = 0;i < N_Z;++i){

      for(int n = 0;n < n_max;++n){

         for(int l = 0;l <= l_max;++l)
            delete [] inlm2s[i][n][l];

         delete [] inlm2s[i][n];

      }

      delete [] inlm2s[i];

   }

   delete [] inlm2s;

   for(int i = 0;i < dim;++i)
      delete [] s2t[i];

   delete [] s2t;

}

/** 
 * Standard constructor: allocates the matrices and fills them with the correct elements by transforming a CartInt object
 * @param ci input CartInt object
 */
SphInt::SphInt(const CartInt &ci){ 

   S = new Matrix(dim);
   T = new Matrix(dim);
   U = new Matrix(dim);

   V = new Matrix(dim*dim);

   int i,n_i,l_i,m_i;
   int j,n_j,l_j,m_j;

   //start with overlap
   for(int s_i = 0;s_i < dim;++s_i){

      i = s2inlm[s_i][0];
      n_i = s2inlm[s_i][1];
      l_i = s2inlm[s_i][2];
      m_i = s2inlm[s_i][3];

      Transform tf_i(i,n_i,l_i,m_i);

      for(int s_j = s_i;s_j < dim;++s_j){

         j = s2inlm[s_j][0];
         n_j = s2inlm[s_j][1];
         l_j = s2inlm[s_j][2];
         m_j = s2inlm[s_j][3];

         Transform tf_j(j,n_j,l_j,m_j);

         complex<double> c_s(0.0,0.0);

         //S
         for(int d_i = 0;d_i < tf_i.gdim();++d_i)
            for(int d_j = 0;d_j < tf_j.gdim();++d_j)
               c_s += conj(tf_i.gcoef(d_i)) * tf_j.gcoef(d_j) * (ci.gS())(tf_i.gind(d_i),tf_j.gind(d_j));

         (*S)(s_i,s_j) = real(c_s);

         complex<double> c_t(0.0,0.0);

         //T
         for(int d_i = 0;d_i < tf_i.gdim();++d_i)
            for(int d_j = 0;d_j < tf_j.gdim();++d_j)
               c_t += conj(tf_i.gcoef(d_i)) * tf_j.gcoef(d_j) * (ci.gT())(tf_i.gind(d_i),tf_j.gind(d_j));

         (*T)(s_i,s_j) = real(c_t);

         complex<double> c_u(0.0,0.0);

         //U
         for(int d_i = 0;d_i < tf_i.gdim();++d_i)
            for(int d_j = 0;d_j < tf_j.gdim();++d_j)
               c_u += conj(tf_i.gcoef(d_i)) * tf_j.gcoef(d_j) * (ci.gU())(tf_i.gind(d_i),tf_j.gind(d_j));

         (*U)(s_i,s_j) = real(c_u);

      }

   }

   S->symmetrize();
   T->symmetrize();
   U->symmetrize();

   int s_i,s_j,s_k,s_l;
   int k,n_k,l_k,m_k;
   int l,n_l,l_l,m_l;

   for(int t_i = 0;t_i < dim*dim;++t_i){

      s_i = t2s[t_i][0];
      s_j = t2s[t_i][1];

      i = s2inlm[s_i][0];
      n_i = s2inlm[s_i][1];
      l_i = s2inlm[s_i][2];
      m_i = s2inlm[s_i][3];

      Transform tf_i(i,n_i,l_i,m_i);
      
      j = s2inlm[s_j][0];
      n_j = s2inlm[s_j][1];
      l_j = s2inlm[s_j][2];
      m_j = s2inlm[s_j][3];

      Transform tf_j(j,n_j,l_j,m_j);

      for(int t_j = t_i;t_j < dim*dim;++t_j){

         s_k = t2s[t_j][0];
         s_l = t2s[t_j][1];

         k = s2inlm[s_k][0];
         n_k = s2inlm[s_k][1];
         l_k = s2inlm[s_k][2];
         m_k = s2inlm[s_k][3];

         Transform tf_k(k,n_k,l_k,m_k);

         l = s2inlm[s_l][0];
         n_l = s2inlm[s_l][1];
         l_l = s2inlm[s_l][2];
         m_l = s2inlm[s_l][3];

         Transform tf_l(l,n_l,l_l,m_l);

         complex<double> c_v(0.0,0.0);

         for(int d_i = 0;d_i < tf_i.gdim();++d_i)
            for(int d_j = 0;d_j < tf_j.gdim();++d_j)
               for(int d_k = 0;d_k < tf_k.gdim();++d_k)
                  for(int d_l = 0;d_l < tf_l.gdim();++d_l){

                     c_v += conj(tf_i.gcoef(d_i)) * conj(tf_j.gcoef(d_j)) * tf_k.gcoef(d_k) * tf_l.gcoef(d_l)
                     
                        * ci.gV(tf_i.gind(d_i),tf_j.gind(d_j),tf_k.gind(d_k),tf_l.gind(d_l));

                  }

         if(fabs(c_v.imag()) > 1.0e-15)
               cout << t_i << "\t" << t_j << "\t|\t" << c_v << endl;

         (*V)(t_i,t_j) = real(c_v);

      }

   }

   V->symmetrize();

}

/** 
 * copy constructor
 * @param si_c SphInt object to be copied in the newly constructed object
 */
SphInt::SphInt(const SphInt &si_c){ 

   S = new Matrix(si_c.gS());
   T = new Matrix(si_c.gT());
   U = new Matrix(si_c.gU());

   V = new Matrix(si_c.gV());

}

/**
 * standard destructor
 */
SphInt::~SphInt(){ 

   delete S;
   delete T;
   delete U;

   delete V;

}

/** 
 * @return the overlapmatrix, const version
 */
const Matrix &SphInt::gS() const { 

   return *S;

}

/** 
 * @return the overlapmatrix
 */
Matrix &SphInt::gS() { 

   return *S;

}

/** 
 * @return the kinetic energy matrix, const version
 */
const Matrix &SphInt::gT() const { 

   return *T; 
}

/** 
 * @return the kinetic energy matrix
 */
Matrix &SphInt::gT() { 

   return *T;

}

/** 
 * @return the nuclear attraction matrix, const version
 */
const Matrix &SphInt::gU() const { 

   return *U; 
}

/** 
 * @return the nuclear attraction matrix
 */
Matrix &SphInt::gU() { 

   return *U;

}

/** 
 * @return the electronic repulsion matrix
 */
const Matrix &SphInt::gV() const { 

   return *V; 
}

/** 
 * @return the electronic repulsion matrix
 */
Matrix &SphInt::gV() { 

   return *V;

}

/**
 * @return the dimension of spatial sp space
 */
int SphInt::gdim() {

   return dim;

}

/**
 * static function
 * @return nr of electrons
 */
int SphInt::gN(){

   return N;

}

ostream &operator<<(ostream &output,SphInt &si_p){

   output << endl;
   output << "Overlap Matrix:" << endl;
   output << endl;

   for(int s_i = 0;s_i < si_p.gdim();++s_i)
      for(int s_j = s_i;s_j < si_p.gdim();++s_j){

         output << si_p.s2inlm[s_i][0] << "\t" << si_p.s2inlm[s_i][1] << "\t" << si_p.s2inlm[s_i][2] << "\t" << si_p.s2inlm[s_i][3]

            << "\t|\t" << si_p.s2inlm[s_j][0] << "\t" << si_p.s2inlm[s_j][1] << "\t" << si_p.s2inlm[s_j][2]

            << "\t" << si_p.s2inlm[s_j][3] << "\t|\t" << (si_p.gS())(s_i,s_j) << endl;

      }

   output << endl;
   output << "Kinetic energy:" << endl;
   output << endl;

   for(int s_i = 0;s_i < si_p.gdim();++s_i)
      for(int s_j = s_i;s_j < si_p.gdim();++s_j){

         output << si_p.s2inlm[s_i][0] << "\t" << si_p.s2inlm[s_i][1] << "\t" << si_p.s2inlm[s_i][2] << "\t" << si_p.s2inlm[s_i][3]

            << "\t|\t" << si_p.s2inlm[s_j][0] << "\t" << si_p.s2inlm[s_j][1] << "\t" << si_p.s2inlm[s_j][2]

            << "\t" << si_p.s2inlm[s_j][3] << "\t|\t" << (si_p.gT())(s_i,s_j) << endl;

      }

   output << endl;
   output << "Nuclear attraction:" << endl;
   output << endl;

   for(int s_i = 0;s_i < si_p.gdim();++s_i)
      for(int s_j = s_i;s_j < si_p.gdim();++s_j){

         output << si_p.s2inlm[s_i][0] << "\t" << si_p.s2inlm[s_i][1] << "\t" << si_p.s2inlm[s_i][2] << "\t" << si_p.s2inlm[s_i][3]

            << "\t|\t" << si_p.s2inlm[s_j][0] << "\t" << si_p.s2inlm[s_j][1] << "\t" << si_p.s2inlm[s_j][2]

            << "\t" << si_p.s2inlm[s_j][3] << "\t|\t" << (si_p.gU())(s_i,s_j) << endl;

      }

   output << endl;
   output << "Electronic repulsion energy:" << endl;
   output << endl;

   int s_i,s_j,s_k,s_l;

   for(int t_i = 0;t_i < si_p.dim*si_p.dim;++t_i){

      s_i = si_p.t2s[t_i][0];
      s_j = si_p.t2s[t_i][1];

      for(int t_j = t_i;t_j < si_p.dim*si_p.dim;++t_j){

         s_k = si_p.t2s[t_j][0];
         s_l = si_p.t2s[t_j][1];

         output << "[\t" << si_p.s2inlm[s_i][0] << "\t" << si_p.s2inlm[s_i][1] << "\t" << si_p.s2inlm[s_i][2] << "\t" << si_p.s2inlm[s_i][3] << "\t|\t"

            << si_p.s2inlm[s_j][0] << "\t" << si_p.s2inlm[s_j][1] << "\t" << si_p.s2inlm[s_j][2] << "\t" << si_p.s2inlm[s_j][3] << "\t]\t||\t[" 

            << si_p.s2inlm[s_k][0] << "\t" << si_p.s2inlm[s_k][1] << "\t" << si_p.s2inlm[s_k][2] << "\t" << si_p.s2inlm[s_k][3] << "\t|\t"

            << si_p.s2inlm[s_l][0] << "\t" << si_p.s2inlm[s_l][1] << "\t" << si_p.s2inlm[s_l][2] << "\t" << si_p.s2inlm[s_l][3] << "\t]\t"
            
            << (si_p.gV())(t_i,t_j) << endl;

      }
   }

   return output;

}

/**
 * orthogonalizes the basis: inverse sqrt of S
 */
void SphInt::orthogonalize() {

   //first inverse sqrt of S
   S->sqrt(-1);

   Matrix T_copy(dim);
   Matrix U_copy(dim);

   T_copy = 0.0;
   U_copy = 0.0;

   //transform T
   for(int i = 0;i < dim;++i)
      for(int j = 0;j < dim;++j){

         for(int k = 0;k < dim;++k){

            T_copy(i,j) += (*S)(i,k) * (*T)(k,j);
            U_copy(i,j) += (*S)(i,k) * (*U)(k,j);

         }

      }

   *T = 0.0;
   *U = 0.0;

   for(int i = 0;i < dim;++i)
      for(int j = i;j < dim;++j){

         for(int k = 0;k < dim;++k){

            (*T)(i,j) += T_copy(i,k) * (*S)(k,j);
            (*U)(i,j) += U_copy(i,k) * (*S)(k,j);

         }

      }

   Matrix V_copy(dim*dim);

   V_copy = 0.0;

   int a,b,c,d;

   //contract a
   for(int i = 0;i < dim*dim;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = 0;j < dim*dim;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         for(int a_ = 0;a_ < dim;++a_)
            V_copy(i,j) += (*S)(a,a_) * (*V)(s2t[a_][b],j); 

      }
   }

   *V = 0.0;

   //contract b
   for(int i = 0;i < dim*dim;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = 0;j < dim*dim;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         for(int b_ = 0;b_ < dim;++b_)
            (*V)(i,j) += (*S)(b,b_) * V_copy(s2t[a][b_],j); 

      }
   }

   V_copy = 0.0;

   //contract c
   for(int i = 0;i < dim*dim;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = 0;j < dim*dim;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         for(int c_ = 0;c_ < dim;++c_)
            V_copy(i,j) += (*V)(i,s2t[c_][d]) * (*S)(c_,c); 

      }
   }

   *V = 0.0;

   //contract d
   for(int i = 0;i < dim*dim;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = 0;j < dim*dim;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         for(int d_ = 0;d_ < dim;++d_)
            (*V)(i,j) += V_copy(i,s2t[c][d_]) * (*S)(d_,d); 

      }
   }

   T->symmetrize();
   U->symmetrize();

}

/**
 * access to the individual elements of the matrices
 */
double SphInt::gS(int i,int j) const {

   return (*S)(i,j);

}

/**
 * access to the individual elements of the matrices
 */
double SphInt::gT(int i,int j) const {

   return (*T)(i,j);

}

/**
 * access to the individual elements of the matrices
 */
double SphInt::gU(int i,int j) const {

   return (*U)(i,j);

}

/**
 * access to the individual elements of the matrices
 */
double SphInt::gV(int a,int b,int c,int d) const {

   return (*V)(s2t[a][b],s2t[c][d]);

}

/**
 * @return the highest orbital angular momentum
 */
int SphInt::gl_max() {

   return l_max;

}

/**
 * static function that allows for access to the private lists.
 * @param s the single-particle index
 * @param option determines what quantumnumber is returned:
 * if option == 0 : return i
 * if option == 1 : return n
 * if option == 2 : return l
 * if option == 3 : return m
 */
int SphInt::gs2inlm(int s,int option) {

   return s2inlm[s][option];

}

/**
 * static function that allows for access to the private lists.
 */
int SphInt::ginlm2s(int i,int n,int l,int m){

   return inlm2s[i][n - l - 1][l][m + l];

}

/**
 * @return the nr of cores
 */
int SphInt::gN_Z() {

   return N_Z;

}

/**
 * @return the highest value of the main quantumnumber n
 */
int SphInt::gn_max() {

   return n_max;

}

/**
 * @param g input SPM single-particle index
 * @return the matrixelement single particle index s when given a SPM single-particle index g
 */
int SphInt::gg2s(int g){

   int m = SPM::gg2ms(g,0);
   int sp_in = SPM::gg2ms(g,1);

   int i = SPM::gs2inl(m,sp_in,0);
   int n = SPM::gs2inl(m,sp_in,1);
   int l = SPM::gs2inl(m,sp_in,2);

   return ginlm2s(i,n,l,m);

}

/* vim: set ts=3 sw=3 expandtab :*/
