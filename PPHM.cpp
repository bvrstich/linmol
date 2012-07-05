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

   l_max = SI_SPM::gl_max();

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

/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * @param S the pph-spin index
 * @param Lz the pph angular momentum projection index
 * @param S_ab The intermediate spinquantumnumber of a and b.
 * @param m_a angular momentum projection of the first particle of the row
 * @param a spatial index of the first particle of the row 
 * @param m_b angular momentum projection of the second particle of the row 
 * @param b spatial index of the second particle of the row 
 * @param m_c angular momentum projection of the hole of the row 
 * @param c spatial index of the hole of the row 
 * @param S_de The intermediate spinquantumnumber of d and e.
 * @param m_d angular momentum projection of the first particle of the column
 * @param d spatial index of the first particle of the column
 * @param m_e angular momentum projection of the second particle of the column
 * @param e spatial index of the second particle of the column
 * @param m_z angular momentum projection of the hole of the column
 * @param z spatial index of the hole of the column
 * @return the number on place PPHM(S,i,j) with the right phase.
 */
double PPHM::operator()(int S,int Lz,int S_ab,int m_a,int a,int m_b,int b,int m_c,int c,int S_de,int m_d,int d,int m_e,int e,int m_z,int z) const {

   if(m_a + m_b + m_c != m_d + m_e + m_z)
      return 0.0;

   if(m_a + m_b + m_c != Lz)
      return 0.0;

   int i,j;

   int ga = SPM::gms2g(m_a,a);
   int gb = SPM::gms2g(m_b,b);
   int gc = SPM::gms2g(m_c,c);

   int phase_i = get_inco(S,Lz,S_ab,ga,gb,gc,i);

   if(phase_i == 0)
      return 0;

   int gd = SPM::gms2g(m_d,d);
   int ge = SPM::gms2g(m_e,e);
   int gz = SPM::gms2g(m_z,z);

   int phase_j = get_inco(S,Lz,S_de,gd,ge,gz,j);

   if(phase_j == 0)
      return 0;

   return phase_i*phase_j* (*this)(SM2B[S][Lz + 3*l_max],i,j);

}

/** 
 * Static member function that gets the pph-index and phase corresponding to the sp indices S,Lz,S_ab,ga,gb,gc.
 * @param S pph-spin index of the state, 0 -> S = 1/2, 1 -> S = 3/2
 * @param Lz pph angular momentum index of the state, 0 -> S = 1/2, 1 -> S = 3/2
 * @param S_ab intermediate spincoupling of a and b. = 0 or 1
 * @param ga first sp orbital
 * @param gb second sp orbital
 * @param gc third sp orbital
 * @param i the corresponding pph index will be stored in this int after calling the function
 * @return the phase needed to get to a normal ordering of indices that corresponds to a pph index i
 */
int PPHM::get_inco(int S,int Lz,int S_ab,int ga,int gb,int gc,int &i){

   if(S == 0){//S = 1/2

      if(S_ab == 0){//symmetric in spatial sp's

         if(ga <= gb)
            i = s2pph[SM2B[S][Lz + 3*l_max]][0][ga][gb][gc];
         else
            i = s2pph[SM2B[S][Lz + 3*l_max]][0][gb][ga][gc];

         return 1;

      }
      else{//antisymmetric in spatial sp's

         if(ga == gb)
            return 0;

         if(ga < gb){

            i = s2pph[SM2B[S][Lz + 3*l_max]][1][ga][gb][gc];

            return 1;

         }
         else{

            i = s2pph[SM2B[S][Lz + 3*l_max]][1][gb][ga][gc];

            return -1;

         }

      }

   }
   else{//S = 3/2

      if(S_ab == 0)
         return 0;

      if(ga == gb)
         return 0;

      if(ga < gb){

         i = s2pph[SM2B[S][Lz + 3*l_max]][1][ga][gb][gc];

         return 1;

      }
      else{

         i = s2pph[SM2B[S][Lz + 3*l_max]][1][gb][ga][gc];

         return -1;

      }

   }

}

/**
 * The spincoupled T2 map, maps a TPM onto a PPHM object. See notes for more info
 * @param tpm Input TPM matrix
 */
void PPHM::T(const TPM &tpm){

   SPM spm;
   spm.bar(1.0/(N - 1.0),tpm);

   int ga,gb,gc,gd,ge,gz;

   int a,b,c,d,e,z;
   int m_a,m_b,m_c,m_d,m_e,m_z;

   int S_ab,S_de;

   double norm_ab,norm_de;
   int sign_ab,sign_de;

   int S;

   for(int B = 0;B < gnr();++B){

      S = B2SM[B][0];

      for(int i = 0;i < gdim(B);++i){

         S_ab = pph2s[B][i][0];

         ga = pph2s[B][i][1];
         gb = pph2s[B][i][2];
         gc = pph2s[B][i][3];

         m_a = SPM::gg2ms(ga,0);
         a = SPM::gg2ms(ga,1);

         m_b = SPM::gg2ms(gb,0);
         b = SPM::gg2ms(gb,1);

         m_c = SPM::gg2ms(gc,0);
         c = SPM::gg2ms(gc,1);

         sign_ab = 1 - 2*S_ab;

         norm_ab = 1.0;

         if(ga == gb)
            norm_ab /= std::sqrt(2.0);

         for(int j = i;j < gdim(B);++j){

            S_de = pph2s[B][j][0];

            gd = pph2s[B][j][1];
            ge = pph2s[B][j][2];
            gz = pph2s[B][j][3];

            m_d = SPM::gg2ms(gd,0);
            d = SPM::gg2ms(gd,1);

            m_e = SPM::gg2ms(ge,0);
            e = SPM::gg2ms(ge,1);

            m_z = SPM::gg2ms(gz,0);
            z = SPM::gg2ms(gz,1);

            sign_de = 1 - 2*S_de;

            norm_de = 1.0;

            if(gd == ge)
               norm_de /= std::sqrt(2.0);

            //start the map:
            (*this)(B,i,j) = 0.0;

            //tp(1)
            if(gc == gz)
               if(S_ab == S_de)
                  (*this)(B,i,j) += tpm(S_ab,m_a+m_b,m_a,a,m_b,b,m_d,d,m_e,e);

            if(ga == gd){

               //sp(1) first term
               if(gb == ge)
                  if(S_ab == S_de)
                     (*this)(B,i,j) += norm_ab * norm_de * spm(-m_c + l_max,c,z);

               //tp(2)
               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g9j(S,Z,S_ab,S_de) * tpm(Z,-m_c+m_e,-m_c,c,m_e,e,-m_z,z,m_b,b);

               ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

               if(c == e && -m_c == m_e)
                  ward *= std::sqrt(2.0);

               if(z == b && -m_z == m_b)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= ward;

            }

            if(gb == gd){

               //sp(1) second term
               if(ga == ge)
                  if(S_ab == S_de)
                     (*this)(B,i,j) += sign_ab * norm_ab * norm_de * spm(-m_c +l_max,c,z);

               //tp(3)
               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g9j(S,Z,S_ab,S_de) * tpm(Z,-m_c+m_e,-m_c,c,m_e,e,-m_z,z,m_a,a);

               ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

               if(c == e && -m_c == m_e)
                  ward *= std::sqrt(2.0);

               if(z == a && -m_z == m_a)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= sign_ab * ward;

            }

            //tp(4)
            if(ga == ge){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g9j(S,Z,S_ab,S_de) * tpm(Z,-m_c+m_d,-m_c,c,m_d,d,-m_z,z,m_b,b);

               ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

               if(c == d && -m_c == m_d)
                  ward *= std::sqrt(2.0);

               if(z == b && -m_z == m_b)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= sign_de * ward;

            }

            //tp(5)
            if(gb == ge){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2.0*Z + 1.0) * Tools::g9j(S,Z,S_ab,S_de) * tpm(Z,-m_c+m_d,-m_c,c,m_d,d,-m_z,z,m_a,a);

               ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

               if(c == d && -m_c == m_d)
                  ward *= std::sqrt(2.0);

               if(z == a && -m_z == m_a)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= sign_ab * sign_de * ward;

            }

         }

      }

   }

   this->symmetrize();

}

/**
 * print the TPM object in a non-axially symmetric form
 */
void PPHM::printnax(const char *filename) const {

   ofstream out(filename);
   out.precision(15);

   int S;
   int S_ab,S_de;

   int ga,gb,gc,gd,ge,gz;

   for(int B = 0;B < gnr();++B){

      S = B2SM[B][0];

      for(int i = 0;i < gdim(B);++i){

         S_ab = pph2s[B][i][0];
         ga = pph2s[B][i][1];
         gb = pph2s[B][i][2];
         gc = pph2s[B][i][3];

         for(int j = i;j < gdim(B);++j){

            S_de = pph2s[B][j][0];
            gd = pph2s[B][j][1];
            ge = pph2s[B][j][2];
            gz = pph2s[B][j][3];

            out << S << "\t" << S_ab << "\t" << ga << "\t" << gb << "\t" << gc << "\t"
            
               << S_de << "\t" << gd << "\t" << ge << "\t" << gz << "\t" << (*this)(B,i,j) << endl;

         }
      }

   }

}

/* vim: set ts=3 sw=3 expandtab :*/
