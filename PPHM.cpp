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

vector< vector< vector<int> > > PPHM::pph2s;
int *****PPHM::s2pph;

vector< vector<int> > PPHM::B2SM;
int **PPHM::SM2B;

int PPHM::M;
int PPHM::N;

int PPHM::l_max;

/**
 * initialize the static lists and variables
 * @param M_in input dimension of sp space
 * @param N_in input nr of particles
 */
void PPHM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

   l_max = SphInt::gl_max();

   //allocate the block list
   SM2B = new int * [2];

   for(int S = 0;S < 2;++S)
      SM2B[S] = new int [6*l_max + 1];

   vector<int> v(4);

   vector<int> bq(2);

   int B = 0;

   //first S == 1/2
   for(int L_z = -3*l_max;L_z <= 3*l_max;++L_z){

      vector< vector<int> > bv;

      //S_ab = 0/1 ga <(=) gb < gc
      for(int S_ab = 0;S_ab < 2;++S_ab)
         for(int ga = 0;ga < M/2;++ga)
            for(int gb = ga + S_ab;gb < M/2;++gb)
               for(int gc = 0;gc < M/2;++gc){

                  if(L_z == SPM::gg2ms(ga,0) + SPM::gg2ms(gb,0) + SPM::gg2ms(gc,0)){//correct L_z projection

                     v[0] = S_ab;
                     v[1] = ga;
                     v[2] = gb;
                     v[3] = gc;

                     bv.push_back(v);

                  }

               }

      if(bv.size() != 0){

         pph2s.push_back(bv);

         bq[0] = 0;
         bq[1] = L_z;

         B2SM.push_back(bq);

         SM2B[0][L_z + 3*l_max] = B;

         ++B;

      }

   }

   //then S == 3/2
   for(int L_z = -3*l_max;L_z <= 3*l_max;++L_z){

      vector< vector<int> > bv;

      //only S_ab = 1 ; ga < gb < gc
      for(int ga = 0;ga < M/2;++ga)
         for(int gb = ga + 1;gb < M/2;++gb)
            for(int gc = 0;gc < M/2;++gc){

               if(L_z == SPM::gg2ms(ga,0) + SPM::gg2ms(gb,0) + SPM::gg2ms(gc,0)){//correct L_z projection

                  v[0] = 1;
                  v[1] = ga;
                  v[2] = gb;
                  v[3] = gc;

                  bv.push_back(v);

               }

            }

      if(bv.size() != 0){

         pph2s.push_back(bv);

         bq[0] = 1;
         bq[1] = L_z;

         B2SM.push_back(bq);

         SM2B[1][L_z + 3*l_max] = B;

         ++B;

      }

   }

   s2pph = new int **** [B2SM.size()];

   for(unsigned int B = 0;B < B2SM.size();++B){

      s2pph[B] = new int *** [2];

      for(int S_ab = 0;S_ab < 2;++S_ab){

         s2pph[B][S_ab] = new int ** [M/2];

         for(int ga = 0;ga < M/2;++ga){

            s2pph[B][S_ab][ga] = new int * [M/2];

            for(int gb = 0;gb < M/2;++gb)
               s2pph[B][S_ab][ga][gb] = new int [M/2];

         }
      }
   }

   for(unsigned int B = 0;B < B2SM.size();++B)
      for(unsigned int d = 0;d < pph2s[B].size();++d)
         s2pph[B][pph2s[B][d][0]][pph2s[B][d][1]][pph2s[B][d][2]][pph2s[B][d][3]] = d;

}

/**
 * deallocate the static lists
 */
void PPHM::clear(){

   for(unsigned int B = 0;B < B2SM.size();++B){

      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int ga = 0;ga < M/2;++ga){

            for(int gb = 0;gb < M/2;++gb)
               delete [] s2pph[B][S_ab][ga][gb];

            delete [] s2pph[B][S_ab][ga];

         }

         delete [] s2pph[B][S_ab];

      }

      delete [] s2pph[B];

   }

   delete [] s2pph;

   for(int S = 0;S < 2;++S)
      delete [] SM2B[S];

   delete [] SM2B;

}

/**
 * standard constructor for a spin and axial symmetrical tp matrix.
 */
PPHM::PPHM() : BlockMatrix(B2SM.size()) {

   for(int B = 0;B < gnr();++B)
      this->setMatrixDim(B,pph2s[B].size(),2*B2SM[B][0] + 2);

}

/**
 * copy constructor:
 * @param pphm_c object that will be copied into this.
 */
PPHM::PPHM(const PPHM &pphm_c) : BlockMatrix(pphm_c){ }

/**
 * destructor
 */
PPHM::~PPHM(){ }

/**
 * @return number of particles
 */
int PPHM::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int PPHM::gM() const{

   return M;

}

ostream &operator<<(ostream &output,const PPHM &pphm_p){

   for(int B = 0;B < pphm_p.gnr();++B){

      output << endl;
      output << "Block " << B << " with S = " << pphm_p.B2SM[B][0] << " and M = " << pphm_p.B2SM[B][1] << endl;
      output << endl;

      int ga,gb,gc,gd,ge,gz;
      int S_ab,S_de;

      for(int pph_i = 0;pph_i < pphm_p.gdim(B);++pph_i){

         S_ab = pphm_p.pph2s[B][pph_i][0];

         ga = pphm_p.pph2s[B][pph_i][1];
         gb = pphm_p.pph2s[B][pph_i][2];
         gc = pphm_p.pph2s[B][pph_i][3];

         for(int pph_j = pph_i;pph_j < pphm_p.gdim(B);++pph_j){

            S_de = pphm_p.pph2s[B][pph_j][0];

            gd = pphm_p.pph2s[B][pph_j][1];
            ge = pphm_p.pph2s[B][pph_j][2];
            gz = pphm_p.pph2s[B][pph_j][3];

            output << "[ " << S_ab << "\t(" << SPM::gg2ms(ga,0) << "," << SPM::gg2ms(ga,1) << ")\t"

               << "(" << SPM::gg2ms(gb,0) << "," << SPM::gg2ms(gb,1) << ")\t(" <<  SPM::gg2ms(gc,0) << "," << SPM::gg2ms(gc,1) << ") ]"

               << "\t|\t[ " << S_de << "\t(" << SPM::gg2ms(gd,0) << "," << SPM::gg2ms(gd,1) << ")\t"

               << "(" << SPM::gg2ms(ge,0) << "," << SPM::gg2ms(ge,1) << ")\t(" << SPM::gg2ms(gz,0) << "," << SPM::gg2ms(gz,1) << ") ]"

               << "\t||\t" << pphm_p(B,pph_i,pph_j) << endl;

         }
      }

   }

   return output;

}
