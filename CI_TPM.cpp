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

int CI_TPM::tp_dim;

vector< vector<int> > CI_TPM::t2s;
int **CI_TPM::s2t;

/** 
 * static function that reads in the input data and makes the corresponding lists and dimensions
 */
void CI_TPM::init(){

   tp_dim = CI_SPM::gdim()*CI_SPM::gdim();

   s2t = new int * [CI_SPM::gdim()];

   for(int i = 0;i < CI_SPM::gdim();++i)
      s2t[i] = new int [CI_SPM::gdim()];

   vector<int> vst(2);

   int iter = 0;

   for(int i = 0;i < CI_SPM::gdim();++i)
      for(int j = 0;j < CI_SPM::gdim();++j){

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
void CI_TPM::clear(){

   for(int i = 0;i < CI_SPM::gdim();++i)
      delete [] s2t[i];

   delete [] s2t;

}

/** 
 * Standard constructor: allocates the matrices and fills them with the correct elements
 */
CI_TPM::CI_TPM() : Matrix(tp_dim) { }

/** 
 * copy constructor
 * @param ci_c CI_TPM object to be copied in the newly constructed object
 */
CI_TPM::CI_TPM(const CI_TPM &ci_c) : Matrix(ci_c){ }

/**
 * standard destructor
 */
CI_TPM::~CI_TPM(){ }


/**
 * access the elements of the matrix using the CI_SPM indices
 * @return the number on place CI_TPM(i,j)
 */
double CI_TPM::operator()(int a,int b,int c,int d) const{

   return (*this)(gs2t(a,b),gs2t(c,d));

}

/**
 * access the elements of the matrix using the physical quantum number i,n,l,x,y,z.
 * @return the number on place CI_TPM(i,j)
 */
double CI_TPM::operator()(int i_a,int n_a,int l_a,int x_a,int y_a,int z_a,int i_b,int n_b,int l_b,int x_b,int y_b,int z_b,

      int i_c,int n_c,int l_c,int x_c,int y_c,int z_c,int i_d,int n_d,int l_d,int x_d,int y_d,int z_d) const {

   return (*this)(CI_SPM::ginlxyz2s(i_a,n_a,l_a,x_a,y_a,z_a),CI_SPM::ginlxyz2s(i_b,n_b,l_b,x_b,y_b,z_b),
   
         CI_SPM::ginlxyz2s(i_c,n_c,l_c,x_c,y_c,z_c),CI_SPM::ginlxyz2s(i_d,n_d,l_d,x_d,y_d,z_d));

}

ostream &operator<<(ostream &output,CI_TPM &ci_p){

   int s_i,s_j,s_k,s_l;

   for(int t_i = 0;t_i < CI_TPM::tp_dim;++t_i){

      s_i = ci_p.gt2s(t_i,0);
      s_j = ci_p.gt2s(t_i,1);

      for(int t_j = 0;t_j < CI_TPM::tp_dim;++t_j){

         s_k = ci_p.gt2s(t_j,0);
         s_l = ci_p.gt2s(t_j,1);

         output << "[\t" << CI_SPM::gs2inlxyz(s_i,0) << "\t" << CI_SPM::gs2inlxyz(s_i,1) << "\t" << CI_SPM::gs2inlxyz(s_i,2)

            << "\t(" << CI_SPM::gs2inlxyz(s_i,3) << "," << CI_SPM::gs2inlxyz(s_i,4) << "," << CI_SPM::gs2inlxyz(s_i,5) << ")\t|\t"

            << CI_SPM::gs2inlxyz(s_j,0) << "\t" << CI_SPM::gs2inlxyz(s_j,1) << "\t" << CI_SPM::gs2inlxyz(s_j,2)

            << "\t(" << CI_SPM::gs2inlxyz(s_j,3) << "," << CI_SPM::gs2inlxyz(s_j,4) << "," << CI_SPM::gs2inlxyz(s_j,5) << ")\t]\t||\t[" 

            << CI_SPM::gs2inlxyz(s_k,0) << "\t" << CI_SPM::gs2inlxyz(s_k,1) << "\t" << CI_SPM::gs2inlxyz(s_k,2)

            << "\t(" << CI_SPM::gs2inlxyz(s_k,3) << "," << CI_SPM::gs2inlxyz(s_k,4) << "," << CI_SPM::gs2inlxyz(s_k,5) << ")\t|\t"

            << CI_SPM::gs2inlxyz(s_l,0) << "\t" << CI_SPM::gs2inlxyz(s_l,1) << "\t" << CI_SPM::gs2inlxyz(s_l,2)

            << "\t(" << CI_SPM::gs2inlxyz(s_l,3) << "," << CI_SPM::gs2inlxyz(s_l,4) << "," << CI_SPM::gs2inlxyz(s_l,5) << ")\t]\t" << ci_p(t_i,t_j) << endl;

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
int CI_TPM::gs2t(int a,int b) {

   return s2t[a][b];

}

/**
 * public static function to access the lists safely from other classes
 * @param t the two-particle index
 * @param option determines which single-particle index to return
 * @return if option == 0, return a, if option == 1 return b
 */
int CI_TPM::gt2s(int t,int option) {

   return t2s[t][option];

}

/**
 * Fill the CI_TPM object with the electronic repulsion integrals
 */
void CI_TPM::calc_eri(){

   int a,b,c,d;

   int i_a,j_a,x_a,y_a,z_a;
   int i_b,j_b,x_b,y_b,z_b;
   int i_c,j_c,x_c,y_c,z_c;
   int i_d,j_d,x_d,y_d,z_d;

   for(int ti = 0;ti < tp_dim;++ti){

      a = t2s[ti][0];
      b = t2s[ti][1];

      //a
      i_a = CI_SPM::gs2inlxyz(a,0);
      j_a = CI_SPM::gs2gauss(a,1);

      x_a = CI_SPM::gs2inlxyz(a,3);
      y_a = CI_SPM::gs2inlxyz(a,4);
      z_a = CI_SPM::gs2inlxyz(a,5);

      //b
      i_b = CI_SPM::gs2inlxyz(b,0);
      j_b = CI_SPM::gs2gauss(b,1);

      x_b = CI_SPM::gs2inlxyz(b,3);
      y_b = CI_SPM::gs2inlxyz(b,4);
      z_b = CI_SPM::gs2inlxyz(b,5);

      for(int tj = ti;tj < tp_dim;++tj){

         c = t2s[tj][0];
         d = t2s[tj][1];

         //c
         i_c = CI_SPM::gs2inlxyz(c,0);
         j_c = CI_SPM::gs2gauss(c,1);

         x_c = CI_SPM::gs2inlxyz(c,3);
         y_c = CI_SPM::gs2inlxyz(c,4);
         z_c = CI_SPM::gs2inlxyz(c,5);

         //d
         i_d = CI_SPM::gs2inlxyz(d,0);
         j_d = CI_SPM::gs2gauss(d,1);

         x_d = CI_SPM::gs2inlxyz(d,3);
         y_d = CI_SPM::gs2inlxyz(d,4);
         z_d = CI_SPM::gs2inlxyz(d,5);

         (*this)(ti,tj) = 0.0;

         for(int k_a = 0;k_a < input::gGaussInfo(i_a).gNcontr(j_a);++k_a)
            for(int k_b = 0;k_b < input::gGaussInfo(i_b).gNcontr(j_b);++k_b)
               for(int k_c = 0;k_c < input::gGaussInfo(i_c).gNcontr(j_c);++k_c)
                  for(int k_d = 0;k_d < input::gGaussInfo(i_d).gNcontr(j_d);++k_d){

                     (*this)(ti,tj) += input::gGaussInfo(i_a).gprefactors(j_a,k_a) * input::gGaussInfo(i_b).gprefactors(j_b,k_b) 

                        * input::gGaussInfo(i_c).gprefactors(j_c,k_c) * input::gGaussInfo(i_d).gprefactors(j_d,k_d)

                        * Tools::gnorm(a,k_a) * Tools::gnorm(b,k_b) * Tools::gnorm(c,k_c) * Tools::gnorm(d,k_d)

                        * Tools::gV()(i_a,j_a,k_a,x_a,y_a,z_a,i_c,j_c,k_c,x_c,y_c,z_c,i_b,j_b,k_b,x_b,y_b,z_b,i_d,j_d,k_d,x_d,y_d,z_d);

                  }

      }
   }

   this->symmetrize();

}

/**
 * @return the dimension of the CI_TPM object
 */
int CI_TPM::gtp_dim(){

   return tp_dim;

}
