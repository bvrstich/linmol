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

vector< vector< vector<int> > > PHM::ph2s;
int ***PHM::s2ph;

vector< vector<int> > PHM::B2SM;
int **PHM::SM2B;

int PHM::M;
int PHM::N;

int PHM::l_max;

/**
 * initialize the static lists and variables
 * @param M_in input dimension of sp space
 * @param N_in input nr of particles
 */
void PHM::init(int M_in,int N_in){

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
            for(int j = 0;j < M/2;++j){

               if(L_z == SPM::gg2ms(i,0) + SPM::gg2ms(j,0)){//correct L_z projection

                  v[0] = i;
                  v[1] = j;

                  bv.push_back(v);

               }

            }

         if(bv.size() != 0){

            ph2s.push_back(bv);

            v[0] = S;
            v[1] = L_z;

            B2SM.push_back(v);

            SM2B[S][L_z + 2*l_max] = B;

            ++B;

         }

      }

   s2ph = new int ** [B2SM.size()];

   for(unsigned int B = 0;B < B2SM.size();++B){

      s2ph[B] = new int * [M/2];

      for(int i = 0;i < M/2;++i)
         s2ph[B][i] = new int [M/2];

   }

   for(unsigned int B = 0;B < B2SM.size();++B){

      for(unsigned int ph = 0;ph < ph2s[B].size();++ph)
         s2ph[B][ph2s[B][ph][0]][ph2s[B][ph][1]] = ph;

   }

}

/**
 * deallocate the static lists
 */
void PHM::clear(){

   for(unsigned int B = 0;B < B2SM.size();++B){

      for(int i = 0;i < M/2;++i)
         delete [] s2ph[B][i];

      delete [] s2ph[B];

   }

   delete [] s2ph;

   for(int S = 0;S < 2;++S)
      delete [] SM2B[S];

   delete [] SM2B;

}

/**
 * standard constructor for a spin and axial symmetrical tp matrix.
 */
PHM::PHM() : BlockMatrix(B2SM.size()) {

   for(int B = 0;B < gnr();++B)
      this->setMatrixDim(B,ph2s[B].size(),2*B2SM[B][0] + 1);

}

/**
 * copy constructor:
 * @param phm_c object that will be copied into this.
 */
PHM::PHM(const PHM &phm_c) : BlockMatrix(phm_c){ }

/**
 * destructor
 */
PHM::~PHM(){ }

/**
 * @return number of particles
 */
int PHM::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int PHM::gM() const{

   return M;

}

ostream &operator<<(ostream &output,const PHM &phm_p){

   for(int B = 0;B < phm_p.gnr();++B){

      output << endl;
      output << "Block " << B << " with S = " << phm_p.B2SM[B][0] << " and M = " << phm_p.B2SM[B][1] << endl;
      output << endl;

      int s_i,s_j,s_k,s_l;

      for(int ph_i = 0;ph_i < phm_p.gdim(B);++ph_i){

         s_i = phm_p.ph2s[B][ph_i][0];
         s_j = phm_p.ph2s[B][ph_i][1];

         for(int ph_j = ph_i;ph_j < phm_p.gdim(B);++ph_j){

            s_k = phm_p.ph2s[B][ph_j][0];
            s_l = phm_p.ph2s[B][ph_j][1];

            output << "[ (" << SPM::gg2ms(s_i,0) << "," << SPM::gg2ms(s_i,1) << ")\t" << "(" << SPM::gg2ms(s_j,0) << "," << SPM::gg2ms(s_j,1) << ") ]"

               << "\t|\t[ (" << SPM::gg2ms(s_k,0) << "," << SPM::gg2ms(s_k,1) << ")\t" << "(" << SPM::gg2ms(s_l,0) << "," << SPM::gg2ms(s_l,1) << ") ]"

               << "\t||\t" << phm_p(B,ph_i,ph_j) << endl;

         }
      }

   }

   return output;

}

/**
 * access the elements of the the blocks in sp mode, the symmetry or antisymmetry of the blocks is automatically accounted for:\n\n
 * Antisymmetrical for S = 1, symmetrical in the sp orbitals for S = 0\n\n
 * @param B The block index
 * @param a first sp index that forms the ph row index i of spin S, together with b
 * @param b second sp index that forms the ph row index i of spin S, together with a
 * @param c first sp index that forms the ph column index j of spin S, together with d
 * @param d second sp index that forms the ph column index j of spin S, together with c
 * @return the number on place PHM(B,i,j) with the right phase.
 */
double PHM::operator()(int B,int a,int b,int c,int d) const{

   return (*this)(B2SM[B][0],B2SM[B][1],SPM::gg2ms(a,0),SPM::gg2ms(a,1),SPM::gg2ms(b,0),SPM::gg2ms(b,1),
   
         SPM::gg2ms(c,0),SPM::gg2ms(c,1),SPM::gg2ms(d,0),SPM::gg2ms(d,1));

}

/**
 * access the elements of the the blocks in sp mode, the symmetry or antisymmetry of the blocks is automatically accounted for:\n\n
 * Antisymmetrical for S = 1, symmetrical in the sp orbitals for S = 0\n\n
 * @param S the two-particle spin
 * @param L_z the two-particle angularmomentum projection
 * @param m_a angular momentum projection of first index
 * @param a first index
 * @param m_b angular momentum projection of second index
 * @param b second index
 * @param m_c angular momentum projection of third index
 * @param c third index
 * @param m_d angular momentum projection of fourth index
 * @param d fourth index
 * @return the number on place PHM(B,i,j) with the right phase.
 */
double PHM::operator()(int S,int L_z,int m_a,int a,int m_b,int b,int m_c,int c,int m_d,int d) const{

   if(m_a + m_b != m_c + m_d)
      return 0.0;

   if(L_z != m_a + m_b)
      return 0.0;

   int ga = SPM::gms2g(m_a,a);
   int gb = SPM::gms2g(m_b,b);
   int gc = SPM::gms2g(m_c,c);
   int gd = SPM::gms2g(m_d,d);

   int B = SM2B[S][L_z + 2*l_max];

   int i = s2ph[B][ga][gb];
   int j = s2ph[B][gc][gd];

   return (*this)(B,i,j);

}

/**
 * The G map, maps a TPM object on a PHM object.
 * @param tpm input TPM
 */
void PHM::G(const TPM &tpm){

   //construct the SPM corresponding to the TPM
   SPM spm;
   spm.bar(1.0/(N - 1.0),tpm);

   int ga,gb,gc,gd;

   int S;

   int a,b,c,d;
   int m_a,m_b,m_c,m_d;

   for(int B = 0;B < gnr();++B){

      S = B2SM[B][0];

      for(int i = 0;i < gdim(B);++i){

         ga = ph2s[B][i][0];
         gb = ph2s[B][i][1];

         m_a = SPM::gg2ms(ga,0);
         a = SPM::gg2ms(ga,1);

         m_b = SPM::gg2ms(gb,0);
         b = SPM::gg2ms(gb,1);

         for(int j = i;j < gdim(B);++j){

            gc = ph2s[B][j][0];
            gd = ph2s[B][j][1];

            m_c = SPM::gg2ms(gc,0);
            c = SPM::gg2ms(gc,1);

            m_d = SPM::gg2ms(gd,0);
            d = SPM::gg2ms(gd,1);

            //tp part
            (*this)(B,i,j) = -Tools::g6j(0,0,S,0) * tpm(0,m_a-m_d,m_a,a,-m_d,d,m_c,c,-m_b,b) 
            
               - 3.0 * Tools::g6j(0,0,S,1) * tpm(1,m_a-m_d,m_a,a,-m_d,d,m_c,c,-m_b,b);

            //norm
            if(a == d && m_a == -m_d)
               (*this)(B,i,j) *= std::sqrt(2.0);

            if(c == b && m_c == -m_b)
               (*this)(B,i,j) *= std::sqrt(2.0);

            //sp part
            if(gb == gd)
               (*this)(B,i,j) += spm(m_a + l_max,a,c);

         }
      }

   }

   this->symmetrize();

}

/**
 * The bar function that maps a PPHM object onto a PHM object by tracing away the first pair of incdices of the PPHM
 * @param pphm Input PPHM object
 */
void PHM::bar(const PPHM &pphm){

   int ga,gb,gc,gd;
   int a,b,c,d,l;

   int m_a,m_b,m_c,m_d,m_l;

   int S;

   double ward,hard;

   for(int B = 0;B < gnr();++B){//loop over blocks of PHM

      S = B2SM[B][0];

      for(int i = 0;i < gdim(B);++i){

         ga = ph2s[B][i][0];
         gb = ph2s[B][i][1];

         m_a = SPM::gg2ms(ga,0);
         a = SPM::gg2ms(ga,1);

         m_b = SPM::gg2ms(gb,0);
         b = SPM::gg2ms(gb,1);

         for(int j = i;j < gdim(B);++j){

            gc = ph2s[B][j][0];
            gd = ph2s[B][j][1];

            m_c = SPM::gg2ms(gc,0);
            c = SPM::gg2ms(gc,1);

            m_d = SPM::gg2ms(gd,0);
            d = SPM::gg2ms(gd,1);

            (*this)(B,i,j) = 0.0;

            //first the S = 1/2 block of the PPHM matrix
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_de = 0;S_de < 2;++S_de){

                  ward = 2.0 * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) ) * Tools::g6j(0,0,S,S_ab) * Tools::g6j(0,0,S,S_de);

                  for(int gl = 0;gl < M/2;++gl){

                     m_l = SPM::gg2ms(gl,0);
                     l = SPM::gg2ms(gl,1);

                     hard = ward * pphm(0,m_l+m_a+m_b,S_ab,m_l,l,m_a,a,m_b,b,S_de,m_l,l,m_c,c,m_d,d);

                     //norms
                     if(gl == ga)
                        hard *= std::sqrt(2.0);

                     if(gl == gc)
                        hard *= std::sqrt(2.0);

                     (*this)(B,i,j) += hard;

                  }

               }

            //then the S = 3/2 block
            if(S == 1){

               for(int gl = 0;gl < M/2;++gl){

                  m_l = SPM::gg2ms(gl,0);
                  l = SPM::gg2ms(gl,1);

                  (*this)(B,i,j) += 4.0/3.0 * pphm(1,m_a+m_b+m_l,1,m_l,l,m_a,a,m_b,b,1,m_l,l,m_c,c,m_d,d);

               }

            }

         }
      }

   }

   this->symmetrize();

}
