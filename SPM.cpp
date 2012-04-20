#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

vector< vector<int> > *SPM::s2inl;
int ****SPM::inl2s;

int SPM::M;
int SPM::N;

int SPM::l_max;

/**
 * initialize the static lists and variables
 * @param M_in input dimension of sp space
 * @param N_in input nr of particles
 */
void SPM::init(int M_in,int N_in) {

   M = M_in;
   N = N_in;

   l_max = SphInt::gl_max();

   //allocate
   s2inl = new vector< vector<int> > [2*l_max + 1];

   inl2s = new int *** [2*l_max + 1];
   
   for(int m = 0;m < 2*l_max + 1;++m){

      inl2s[m] = new int ** [SphInt::gN_Z()];

      for(int i = 0;i < SphInt::gN_Z();++i){

         inl2s[m][i] = new int * [SphInt::gn_max()];

         for(int n = 0;n < SphInt::gn_max();++n)
            inl2s[m][i][n] = new int [SphInt::gl_max() + 1];

      }

   }

   //construct
   vector<int> v(3);

   for(int m = -l_max;m <= l_max;++m){

      int iter = 0;

      for(int s = 0;s < SphInt::gdim();++s){

         if(SphInt::gs2inlm(s,3) == m){

            v[0] = SphInt::gs2inlm(s,0);//i
            v[1] = SphInt::gs2inlm(s,1);//n
            v[2] = SphInt::gs2inlm(s,2);//l

            s2inl[m + l_max].push_back(v);

            inl2s[m + l_max][v[0]][v[1] - v[2] - 1][v[2]] = iter;

            ++iter;

         }

      }

   }

}

/**
 * deallocate the static lists
 */
void SPM::clear(){

   delete [] s2inl;

   for(int m = 0;m < 2*l_max + 1;++m){

      for(int i = 0;i < SphInt::gN_Z();++i){

         for(int n = 0;n < SphInt::gn_max();++n)
            delete [] inl2s[m][i][n];

         delete [] inl2s[m][i];

      }

      delete [] inl2s[m];

   }

   delete [] inl2s;

}

/**
 * standard constructor, initializes a BlockMatrix with 2*l_max + 1 blocks corresponding to the possible m projections
 */
SPM::SPM() : BlockMatrix(2*l_max + 1) {
   
   for(int m = 0;m < 2*l_max + 1;++m)
      this->setMatrixDim(m,s2inl[m].size(),2);
   
}

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(const SPM &spm_copy) : BlockMatrix(spm_copy) { }

/**
 * destructor
 */
SPM::~SPM(){ }

/**
 * @return nr of particles
 */
int SPM::gN() const{

   return N;

}

/**
 * @return dimension of sp space
 */
int SPM::gM() const{

   return M;

}

ostream &operator<<(ostream &output,const SPM &spm_p){

   int i,n_i,l_i;
   int j,n_j,l_j;

   for(int m = -spm_p.gl_max();m <= spm_p.gl_max();++m){

      output << endl;
      output << "Block with L_z = " << m << endl;
      output << endl;

      for(int s_i = 0;s_i < spm_p.gdim(m + spm_p.gl_max());++s_i){

         i = spm_p.gs2inl(m,s_i,0);
         n_i = spm_p.gs2inl(m,s_i,1);
         l_i = spm_p.gs2inl(m,s_i,2);

         for(int s_j = s_i;s_j < spm_p.gdim(m + spm_p.gl_max());++s_j){

            j = spm_p.gs2inl(m,s_j,0);
            n_j = spm_p.gs2inl(m,s_j,1);
            l_j = spm_p.gs2inl(m,s_j,2);

            output << "( " << i << " " << n_i << " " << l_i << ")\t|\t(" << j << " " << n_j << " " << l_j
            
               << ")\t||\t" << spm_p(m+spm_p.gl_max(),s_i,s_j) << endl;

         }
      }

   }

   return output;

}

/**
 * @return the highest orbital angular momentum in the basisset
 */
int SPM::gl_max(){

   return l_max;

}

/**
 * static function that allows for access to the private lists.
 * @param m block index corresponding to the z projection of angular mommentum
 * @param s the single-particle index
 * @param option determines what quantumnumber is returned:
 * if option == 0 : return i
 * if option == 1 : return n
 * if option == 2 : return l
 */
int SPM::gs2inl(int m,int s,int option) {

   return s2inl[m + l_max][s][option];

}
