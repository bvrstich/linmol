#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <hdf5.h>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::vector;

#include "include.h"

int SI_TPM::tp_dim;

vector< vector<int> > SI_TPM::t2s;
int **SI_TPM::s2t;

/** 
 * static function that reads in the input data and makes the corresponding lists and dimensions
 */
void SI_TPM::init(){

   tp_dim = SI_SPM::gdim()*SI_SPM::gdim();

   s2t = new int * [SI_SPM::gdim()];

   for(int i = 0;i < SI_SPM::gdim();++i)
      s2t[i] = new int [SI_SPM::gdim()];

   vector<int> vst(2);

   int iter = 0;

   for(int i = 0;i < SI_SPM::gdim();++i)
      for(int j = 0;j < SI_SPM::gdim();++j){

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
void SI_TPM::clear(){

   for(int i = 0;i < SI_SPM::gdim();++i)
      delete [] s2t[i];

   delete [] s2t;

}

/** 
 * Standard constructor: allocates the matrices and fills them with the correct elements
 */
SI_TPM::SI_TPM() : Matrix(tp_dim) { }

/** 
 * copy constructor
 * @param si_c SI_TPM object to be copied in the newly constructed object
 */
SI_TPM::SI_TPM(const SI_TPM &si_c) : Matrix(si_c){ }

/**
 * standard destructor
 */
SI_TPM::~SI_TPM(){ }


/**
 * access the elements of the matrix using the SI_SPM indices
 * @return the number on place SI_TPM(i,j)
 */
double SI_TPM::operator()(int a,int b,int c,int d) const{

   return (*this)(gs2t(a,b),gs2t(c,d));

}

/**
 * access the elements of the matrix using the physical quantum number i,n,l,x,y,z.
 * @return the number on place SI_TPM(i,j)
 */
double SI_TPM::operator()(int i_a,int n_a,int l_a,int m_a,int i_b,int n_b,int l_b,int m_b,

      int i_c,int n_c,int l_c,int m_c,int i_d,int n_d,int l_d,int m_d) const {

   return (*this)(SI_SPM::ginlm2s(i_a,n_a,l_a,m_a),SI_SPM::ginlm2s(i_b,n_b,l_b,m_b),SI_SPM::ginlm2s(i_c,n_c,l_c,m_c),SI_SPM::ginlm2s(i_d,n_d,l_d,m_d));

}

ostream &operator<<(ostream &output,SI_TPM &si_p){

   int s_i,s_j,s_k,s_l;

   for(int t_i = 0;t_i < si_p.tp_dim;++t_i){

      s_i = si_p.t2s[t_i][0];
      s_j = si_p.t2s[t_i][1];

      for(int t_j = t_i;t_j < si_p.tp_dim;++t_j){

         s_k = si_p.t2s[t_j][0];
         s_l = si_p.t2s[t_j][1];

         output << "[\t" << SI_SPM::gs2inlm(s_i,0) << "\t" << SI_SPM::gs2inlm(s_i,1) << "\t" << SI_SPM::gs2inlm(s_i,2) << "\t" << SI_SPM::gs2inlm(s_i,3)
         
            << "\t|\t" << SI_SPM::gs2inlm(s_j,0) << "\t" << SI_SPM::gs2inlm(s_j,1) << "\t" << SI_SPM::gs2inlm(s_j,2) << "\t" << SI_SPM::gs2inlm(s_j,3)
            
            << "\t]\t||\t[" << SI_SPM::gs2inlm(s_k,0) << "\t" << SI_SPM::gs2inlm(s_k,1) << "\t" << SI_SPM::gs2inlm(s_k,2) << "\t"
            
            << SI_SPM::gs2inlm(s_k,3) << "\t|\t" << SI_SPM::gs2inlm(s_l,0) << "\t" << SI_SPM::gs2inlm(s_l,1) << "\t" << SI_SPM::gs2inlm(s_l,2)
            
            << "\t" << SI_SPM::gs2inlm(s_l,3) << "\t]\t" << si_p(t_i,t_j) << endl;

      }
   }

   return output;

}

/**
 * public static function to access the lists safely from other classes
 * @param a the first single particle index
 * @param b the second single particle index
 * @return the two-particle index corresponding to a,b
 */
int SI_TPM::gs2t(int a,int b) {

   return s2t[a][b];

}

/**
 * public static function to access the lists safely from other classes
 * @param t the two-particle index
 * @param option determines which single-particle index to return
 * @return if option == 0, return a, if option == 1 return b
 */
int SI_TPM::gt2s(int t,int option) {

   return t2s[t][option];

}

/**
 * Fill *this object by transforming a cartesian CI_TPM object to a spherical one
 * @param ci input CI_TPM object
 */
void SI_TPM::transform(const CI_TPM &ci){

   int a,b,c,d;

   int i_a,n_a,l_a,m_a;
   int i_b,n_b,l_b,m_b;
   int i_c,n_c,l_c,m_c;
   int i_d,n_d,l_d,m_d;

   for(int t_i = 0;t_i < tp_dim;++t_i){

      a = t2s[t_i][0];
      b = t2s[t_i][1];

      i_a = SI_SPM::gs2inlm(a,0);
      n_a = SI_SPM::gs2inlm(a,1);
      l_a = SI_SPM::gs2inlm(a,2);
      m_a = SI_SPM::gs2inlm(a,3);

      Transform tf_a(i_a,n_a,l_a,m_a);
      
      i_b = SI_SPM::gs2inlm(b,0);
      n_b = SI_SPM::gs2inlm(b,1);
      l_b = SI_SPM::gs2inlm(b,2);
      m_b = SI_SPM::gs2inlm(b,3);

      Transform tf_b(i_b,n_b,l_b,m_b);

      for(int t_j = t_i;t_j < tp_dim;++t_j){

         c = t2s[t_j][0];
         d = t2s[t_j][1];

         i_c = SI_SPM::gs2inlm(c,0);
         n_c = SI_SPM::gs2inlm(c,1);
         l_c = SI_SPM::gs2inlm(c,2);
         m_c = SI_SPM::gs2inlm(c,3);

         Transform tf_c(i_c,n_c,l_c,m_c);

         i_d = SI_SPM::gs2inlm(d,0);
         n_d = SI_SPM::gs2inlm(d,1);
         l_d = SI_SPM::gs2inlm(d,2);
         m_d = SI_SPM::gs2inlm(d,3);

         Transform tf_d(i_d,n_d,l_d,m_d);

         complex<double> c(0.0,0.0);

         for(int d_a = 0;d_a < tf_a.gdim();++d_a)
            for(int d_b = 0;d_b < tf_b.gdim();++d_b)
               for(int d_c = 0;d_c < tf_c.gdim();++d_c)
                  for(int d_d = 0;d_d < tf_d.gdim();++d_d){

                     c += conj(tf_a.gcoef(d_a)) * conj(tf_b.gcoef(d_b)) * tf_c.gcoef(d_c) * tf_d.gcoef(d_d)
                     
                        * ci(tf_a.gind(d_a),tf_b.gind(d_b),tf_c.gind(d_c),tf_d.gind(d_d));

                  }

         (*this)(t_i,t_j) = real(c);

      }

   }

   this->symmetrize();

}

/**
 * @return the dimension of the SI_TPM object
 */
int SI_TPM::gtp_dim(){

   return tp_dim;

}
