#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

vector< vector<int> > *SPM::s2inl;
int ****SPM::inl2s;

vector< vector<int> > SPM::g2ms;
int **SPM::ms2g;

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

   ms2g = new int * [2*l_max + 1];

   for(int m = 0;m < 2*l_max + 1;++m)
      ms2g[m] = new int [s2inl[m].size()];

   vector<int> vg(2);

   int iter = 0;

   for(int m = -l_max;m <= l_max;++m)
      for(unsigned int s = 0;s < s2inl[m + l_max].size();++s){

         vg[0] = m;
         vg[1] = s;

         g2ms.push_back(vg);
         ms2g[m + l_max][s] = iter;

         ++iter;

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

      delete [] ms2g[m];

   }

   delete [] inl2s;

   delete [] ms2g;

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

/**
 * static function that allows for access to the private lists.
 * @param m block index corresponding to the z projection of angular momentum
 * @param i index of core location
 * @param n main quantumnumber of orbital
 * @param l orbital angular momentum quantumnumber
 */
int SPM::ginl2s(int m,int i,int n,int l){

   return inl2s[m + l_max][i][n - l - 1][l];

}

/**
 * function the gives acces to the private ms2g list
 * @param m the l_z projection
 * @param s the single-particle index inside the "m"-block
 */
int SPM::gms2g(int m,int s) {

   return ms2g[m + l_max][s];

}

/**
 * function the gives acces to the private g2ms list
 * @param g the global single-particle index
 * @param option determines what quantumnumber is returned
 * if option == 0 : return m
 * if option == 1 : return s
 */
int SPM::gg2ms(int g,int option) {

   return g2ms[g][option];

}

/**
 * map a TPM object on an SPM object by tracing one pair of indices
 * @param scale number with which to scale the TPM with
 * @param tpm the input TPM object
 */
void SPM::bar(double scale,const TPM &tpm) {

   double ward,hard;

   for(int m = -l_max;m <= l_max;++m){

      for(int a = 0;a < gdim(m + l_max);++a)
         for(int c = a;c < gdim(m + l_max);++c){

            (*this)(l_max + m,a,c) = 0.0;

            for(int S = 0;S < 2;++S){

               hard = 0.0;

               for(int m_b = -l_max;m_b <= l_max;++m_b)
                  for(int b = 0;b < gdim(m_b + l_max);++b){

                     ward = tpm(S,m + m_b,m,a,m_b,b,m,c,m_b,b);

                     if(m == m_b){

                        if(a == b)
                           ward *= std::sqrt(2.0);

                        if(c == b)
                           ward *= std::sqrt(2.0);

                     }

                     hard += ward;

                  }

               (*this)(l_max + m,a,c) += (2.0*S + 1.0) * hard;

            }

            (*this)(l_max + m,a,c) *= scale * 0.5;

         }

   }

   this->symmetrize();

}

/**
 * Trace out a set of indices to create the "bar" matrix of a PHM, slight difference from the bar(TPM) function (normalization of the tp basisset).
 * @param scale the factor u want the SPM to be scaled with
 * @param phm the PHM out of which the SPM will be filled
 */
void SPM::bar(double scale,const PHM &phm){

   double ward;

   for(int m = -l_max;m <= l_max;++m){

      for(int a = 0;a < gdim(m + l_max);++a)
         for(int c = a;c < gdim(m + l_max);++c){

            (*this)(l_max + m,a,c) = 0.0;

            for(int S = 0;S < 2;++S){

               ward = 0.0;

               for(int m_b = -l_max;m_b <= l_max;++m_b)
                  for(int b = 0;b < gdim(m_b + l_max);++b)
                     ward += phm(S,m + m_b,m,a,m_b,b,m,c,m_b,b);

               (*this)(l_max + m,a,c) += (2.0*S + 1.0) * ward;

            }

            (*this)(l_max + m,a,c) *= scale * 0.5;

         }

   }

   this->symmetrize();

}

/** 
 * This bar function maps a PPHM object directly onto a SPM object, scaling it with a factor scale
 * @param scale the scalefactor
 * @param pphm Input PPHM object
 */
void SPM::bar(double scale,const PPHM &pphm){

   int l,k;
   int m_l,m_k;

   for(int m = -l_max;m <= l_max;++m){

      for(int a = 0;a < gdim(m + l_max);++a)
         for(int c = a;c < gdim(m + l_max);++c){

            (*this)(m + l_max,a,c) = 0.0;

            //first S = 1/2 part
            for(int S_lk = 0;S_lk < 2;++S_lk){

               for(int gl = 0;gl < M/2;++gl){

                  m_l = g2ms[gl][0];
                  l = g2ms[gl][1];

                  for(int gk = 0;gk < M/2;++gk){

                     m_k = g2ms[gk][0];
                     k = g2ms[gk][1];

                     if(k == l && m_k == m_l)
                        (*this)(m + l_max,a,c) += 2.0 * pphm(0,m_l+m_k-m,S_lk,m_l,l,m_k,k,-m,a,S_lk,m_l,l,m_k,k,-m,c);
                     else
                        (*this)(m + l_max,a,c) += pphm(0,m_l+m_k-m,S_lk,m_l,l,m_k,k,-m,a,S_lk,m_l,l,m_k,k,-m,c);

                  }
               }

            }

            //then S = 3/2 part:
            for(int gl = 0;gl < M/2;++gl){

               m_l = g2ms[gl][0];
               l = g2ms[gl][1];

               for(int gk = 0;gk < M/2;++gk){

                  m_k = g2ms[gk][0];
                  k = g2ms[gk][1];

                  (*this)(m+l_max,a,c) += 2.0 * pphm(1,m_l+m_k-m,1,m_l,l,m_k,k,-m,a,1,m_l,l,m_k,k,-m,c);

               }
            }

            //scaling
            (*this)(m+l_max,a,c) *= scale;

         }

   }

   this->symmetrize();

}

/**
 * construct the SPM object which, when the dotproduct is taken with a 1DM, gives back the subsystem occupation
 * @param S the overlapmatrix !!positive sqrt!!
 */
void SPM::subocc_op(int core,const Matrix &S){

   *this = 0.0;

   int i_a,n_a,l_a;
   int i_b,n_b,l_b;

   for(int m = -l_max;m <= l_max;++m){

      for(int a = 0;a < gdim(m + l_max);++a){

         i_a = s2inl[m + l_max][a][0];
         n_a = s2inl[m + l_max][a][1];
         l_a = s2inl[m + l_max][a][2];

         for(int b = a;b < gdim(m + l_max);++b){

            i_b = s2inl[m + l_max][b][0];
            n_b = s2inl[m + l_max][b][1];
            l_b = s2inl[m + l_max][b][2];
            
            for(int s = 0;s < SphInt::gdim();++s)
               if(core == SphInt::gs2inlm(s,0))
                  (*this)[m + l_max](a,b) += S(SphInt::ginlm2s(i_a,n_a,l_a,m),s) * S(SphInt::ginlm2s(i_b,n_b,l_b,m),s);

         }
      }

   }

   this->symmetrize();

}

/**
 * construct the SPM object which, when the dotproduct is taken with a 1DM, gives back the subsystem occupation
 * @param S the overlapmatrix !!positive sqrt!!
 */
void SPM::si21dm(const Matrix &S){

   int i_a,n_a,l_a;
   int i_b,n_b,l_b;

   for(int m = -l_max;m <= l_max;++m){

      for(int a = 0;a < gdim(m + l_max);++a){

         i_a = s2inl[m + l_max][a][0];
         n_a = s2inl[m + l_max][a][1];
         l_a = s2inl[m + l_max][a][2];

         for(int b = a;b < gdim(m + l_max);++b){

            i_b = s2inl[m + l_max][b][0];
            n_b = s2inl[m + l_max][b][1];
            l_b = s2inl[m + l_max][b][2];

            int ga = ms2g[m + l_max][a];
            int gb = ms2g[m + l_max][b];

            (*this)[m + l_max](a,b) = S(ga,gb);
            
         }
      }

   }

   this->symmetrize();

}
