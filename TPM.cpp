#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

vector< vector< vector<int> > > TPM::t2s;
int ***TPM::s2t;

vector< vector<int> > TPM::B2SM;
int **TPM::SM2B;

int TPM::M;
int TPM::N;

int TPM::l_max;

/**
 * initialize the static lists and variables
 * @param M_in input dimension of sp space
 * @param N_in input nr of particles
 */
void TPM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

   l_max = SphInt::gl_max();

   //allocate the block list
   SM2B = new int * [2];

   for(int S = 0;S < 2;++S)
      SM2B[S] = new int [4*l_max + 1];

   vector<int> v(2);

   int B = 0;

   for(int S = 0;S < 2;++S)
      for(int L_z = -2*l_max;L_z <= 2*l_max;++L_z){

         vector< vector<int> > bv;

         for(int i = 0;i < M/2;++i)
            for(int j = i + S;j < M/2;++j){

               if(L_z == SPM::gg2ms(i,0) + SPM::gg2ms(j,0)){//correct L_z projection

                  v[0] = i;
                  v[1] = j;

                  bv.push_back(v);

               }

            }

         if(bv.size() != 0){

            t2s.push_back(bv);

            v[0] = S;
            v[1] = L_z;

            B2SM.push_back(v);

            SM2B[S][L_z + 2*l_max] = B;

         }

         ++B;

      }

   s2t = new int ** [B2SM.size()];

   for(unsigned int B = 0;B < B2SM.size();++B){

      s2t[B] = new int * [M/2];

      for(int i = 0;i < M/2;++i)
         s2t[B][i] = new int [M/2];

   }

   for(unsigned int B = 0;B < B2SM.size();++B){

      for(unsigned int t = 0;t < t2s[B].size();++t){

         s2t[B][t2s[B][t][0]][t2s[B][t][1]] = t;
         s2t[B][t2s[B][t][1]][t2s[B][t][0]] = t;

      }

   }

}

/**
 * deallocate the static lists
 */
void TPM::clear(){

   for(unsigned int B = 0;B < B2SM.size();++B){

      for(int i = 0;i < M/2;++i)
         delete [] s2t[B][i];

      delete [] s2t[B];

   }

   delete [] s2t;

   for(int S = 0;S < 2;++S)
      delete [] SM2B[S];

   delete [] SM2B;

}

/**
 * standard constructor for a spinsymmetrical tp matrix: constructs BlockMatrix object with 2 blocks, for S = 0 or 1,
 */
TPM::TPM() : BlockMatrix(B2SM.size()) {

   for(int B = 0;B < gnr();++B)
      this->setMatrixDim(B,t2s[B].size(),2*B2SM[B][0] + 1);

}

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpm_c
 * if counter == 0, the lists containing the relationship between sp and tp basis.
 * @param tpm_c object that will be copied into this.
 */
TPM::TPM(const TPM &tpm_c) : BlockMatrix(tpm_c){ }

/**
 * destructor: if counter == 1 the memory for the static lists t2s en s2t will be deleted.
 * 
 */
TPM::~TPM(){ }

/**
 * @return number of particles
 */
int TPM::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int TPM::gM() const{

   return M;

}

ostream &operator<<(ostream &output,const TPM &tpm_p){

   for(int B = 0;B < tpm_p.gnr();++B){

      output << endl;
      output << "Block " << B << " with S = " << tpm_p.B2SM[B][0] << " and M = " << tpm_p.B2SM[B][1] << endl;
      output << endl;

      int s_i,s_j,s_k,s_l;

      for(int t_i = 0;t_i < tpm_p.gdim(B);++t_i){

         s_i = tpm_p.t2s[B][t_i][0];
         s_j = tpm_p.t2s[B][t_i][1];

         for(int t_j = t_i;t_j < tpm_p.gdim(B);++t_j){

            s_k = tpm_p.t2s[B][t_j][0];
            s_l = tpm_p.t2s[B][t_j][1];

            output << "[ (" << SPM::gg2ms(s_i,0) << "," << SPM::gg2ms(s_i,1) << ")\t" << "(" << SPM::gg2ms(s_j,0) << "," << SPM::gg2ms(s_j,1) << ") ]"

               << "\t|\t[ (" << SPM::gg2ms(s_k,0) << "," << SPM::gg2ms(s_k,1) << ")\t" << "(" << SPM::gg2ms(s_l,0) << "," << SPM::gg2ms(s_l,1) << ") ]"

               << "\t||\t" << tpm_p(B,t_i,t_j) << endl;

         }
      }

   }

   return output;

}
