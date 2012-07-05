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

vector< vector<int> > SI_SPM::s2inlm;
int ****SI_SPM::inlm2s;

int SI_SPM::l_max;
int SI_SPM::n_max;

int SI_SPM::dim;

/** 
 * static function that reads in the input data and makes the corresponding lists and dimensions
 */
void SI_SPM::init(){

   n_max = CI_SPM::gn_max();
   l_max = CI_SPM::gl_max();

   //allocate
   inlm2s = new int *** [input::gN_Z()];

   for(int i = 0;i < input::gN_Z();++i){

      inlm2s[i] = new int ** [n_max];

      for(int n = 0;n < n_max;++n){

         inlm2s[i][n] = new int * [l_max + 1];

         for(int l = 0;l <= l_max;++l)
            inlm2s[i][n][l] = new int [2*l + 1];

      }
   }

   //construct
   vector<int> v(4);

   for(int s = 0;s < CI_SPM::gdim();++s){

      v[0] = CI_SPM::gs2inlxyz(s,0);//i
      v[1] = CI_SPM::gs2inlxyz(s,1);//n
      v[2] = CI_SPM::gs2inlxyz(s,2);//l

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

}

/** 
 * function that deallocates the static members
 */
void SI_SPM::clear(){

   for(int i = 0;i < input::gN_Z();++i){

      for(int n = 0;n < n_max;++n){

         for(int l = 0;l <= l_max;++l)
            delete [] inlm2s[i][n][l];

         delete [] inlm2s[i][n];

      }

      delete [] inlm2s[i];

   }

   delete [] inlm2s;

}

/** 
 * Standard constructor: allocates the matrices and fills them with the correct elements
 */
SI_SPM::SI_SPM() : Matrix(dim) { }

/** 
 * copy constructor
 * @param si_c SI_SPM object to be copied in the newly constructed object
 */
SI_SPM::SI_SPM(const SI_SPM &si_c) : Matrix(si_c){ }

/**
 * standard destructor
 */
SI_SPM::~SI_SPM(){ }

/**
 * access the elements of the matrix using the physical quantum numbers i,n,l,m
 * @return the number on place SI_SPM(i,j)
 */
double SI_SPM::operator()(int i_a,int n_a,int l_a,int m_a,int i_b,int n_b,int l_b,int m_b) const {

   return (*this)(ginlm2s(i_a,n_a,l_a,m_a),ginlm2s(i_b,n_b,l_b,m_b));

}

/**
 * access to the lists from outside of the class
 * @param s single-particle index
 * @param option what qm to return
 */
int SI_SPM::gs2inlm(int s,int option){

   return s2inlm[s][option];

}
/**
 * static function that allows for access to the private lists.
 */
int SI_SPM::ginlm2s(int i,int n,int l,int m){

   return inlm2s[i][n - l - 1][l][m + l];

}

/**
 * @return the highest value of the main quantumnumber n
 */
int SI_SPM::gn_max() {

   return n_max;

}

/**
 * @return the highest value of the angular momentum quantumnumber l
 */
int SI_SPM::gl_max() {

   return l_max;

}

/**
 * @return the dimension of the spherical single-particle space
 */
int SI_SPM::gdim() {

   return dim;

}

/**
 * Fill *this object by transforming a cartesian CI_SPM object to a spherical one
 * @param ci input CI_SPM object
 */
void SI_SPM::transform(const CI_SPM &ci){

   int i_a,n_a,l_a,m_a;
   int i_b,n_b,l_b,m_b;

   //start with overlap
   for(int a = 0;a < dim;++a){

      i_a = s2inlm[a][0];
      n_a = s2inlm[a][1];
      l_a = s2inlm[a][2];
      m_a = s2inlm[a][3];

      Transform tf_a(i_a,n_a,l_a,m_a);

      for(int b = a;b < dim;++b){

         i_b = s2inlm[b][0];
         n_b = s2inlm[b][1];
         l_b = s2inlm[b][2];
         m_b = s2inlm[b][3];

         Transform tf_b(i_b,n_b,l_b,m_b);

         complex<double> c(0.0,0.0);

         for(int d_a = 0;d_a < tf_a.gdim();++d_a)
            for(int d_b = 0;d_b < tf_b.gdim();++d_b)
               c += conj(tf_a.gcoef(d_a)) * tf_b.gcoef(d_b) * ci(tf_a.gind(d_a),tf_b.gind(d_b));

         (*this)(a,b) = real(c);

      }
   }

   this->symmetrize();

}

ostream &operator<<(ostream &output,SI_SPM &si_p){

   for(int s_i = 0;s_i < si_p.gdim();++s_i)
      for(int s_j = s_i;s_j < si_p.gdim();++s_j){

         output << si_p.s2inlm[s_i][0] << "\t" << si_p.s2inlm[s_i][1] << "\t" << si_p.s2inlm[s_i][2] << "\t" << si_p.s2inlm[s_i][3]

            << "\t|\t" << si_p.s2inlm[s_j][0] << "\t" << si_p.s2inlm[s_j][1] << "\t" << si_p.s2inlm[s_j][2]

            << "\t" << si_p.s2inlm[s_j][3] << "\t|\t" << si_p(s_i,s_j) << endl;

      }

   return output;

}
