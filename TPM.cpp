#include <iostream>
#include <cmath>
#include <fstream>
#include <hdf5.h>

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

double TPM::Sa = 1;
double TPM::Sc = 0;

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

            ++B;

         }

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

   init_overlap();

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
 * initialize the overlapmatrix parameters
 */
void TPM::init_overlap(){

#ifdef __Q_CON
   Sa += 1.0;
   Sc += (2.0*N - M)/((N - 1.0)*(N - 1.0));
#endif

#ifdef __G_CON
   Sa += 4.0;
   Sc += (2.0*N - M - 2.0)/((N - 1.0)*(N - 1.0));
#endif

#ifdef __T1_CON
   Sa += M - 4.0;
   Sc -= (M*M + 2.0*N*N - 4.0*M*N - M + 8.0*N - 4.0)/( 2.0*(N - 1.0)*(N - 1.0) );
#endif

#ifdef __T2_CON
   Sa += 5.0*M - 8.0;
   Sc += (2.0*N*N + (M - 2.0)*(4.0*N - 3.0) - M*M)/(2.0*(N - 1.0)*(N - 1.0));
#endif

}

/**
 * standard constructor for a spin and axial symmetrical tp matrix.
 */
TPM::TPM() : BlockMatrix(B2SM.size()) {

   for(int B = 0;B < gnr();++B)
      this->setMatrixDim(B,t2s[B].size(),2*B2SM[B][0] + 1);

}

/**
 * copy constructor:
 * @param tpm_c object that will be copied into this.
 */
TPM::TPM(const TPM &tpm_c) : BlockMatrix(tpm_c){ }

/**
 * destructor
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

/**
 * access the elements of the the blocks in sp mode, the symmetry or antisymmetry of the blocks is automatically accounted for:\n\n
 * Antisymmetrical for S = 1, symmetrical in the sp orbitals for S = 0\n\n
 * @param B The block index
 * @param a first sp index that forms the tp row index i of spin S, together with b
 * @param b second sp index that forms the tp row index i of spin S, together with a
 * @param c first sp index that forms the tp column index j of spin S, together with d
 * @param d second sp index that forms the tp column index j of spin S, together with c
 * @return the number on place TPM(B,i,j) with the right phase.
 */
double TPM::operator()(int B,int a,int b,int c,int d) const{

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
 * @return the number on place TPM(B,i,j) with the right phase.
 */
double TPM::operator()(int S,int L_z,int m_a,int a,int m_b,int b,int m_c,int c,int m_d,int d) const{

   if(m_a + m_b != m_c + m_d)
      return 0.0;

   if(L_z != m_a + m_b)
      return 0.0;

   int ga = SPM::gms2g(m_a,a);
   int gb = SPM::gms2g(m_b,b);
   int gc = SPM::gms2g(m_c,c);
   int gd = SPM::gms2g(m_d,d);

   if(S == 0){

      int B = SM2B[S][L_z + 2*l_max];

      int i = s2t[B][ga][gb];
      int j = s2t[B][gc][gd];

      return (*this)(B,i,j);

   }
   else{

      if( ( ga == gb ) || ( gc == gd ) )
         return 0;
      else{

         int B = SM2B[S][L_z + 2*l_max];

         int i = s2t[B][ga][gb];
         int j = s2t[B][gc][gd];

         int phase = 1;

         if(ga > gb)
            phase *= -1;
         if(gc > gd)
            phase *= -1;

         return phase*(*this)(B,i,j);

      }

   }

}

/**
 * fill the (*this) object with the Hamiltonian that defines the molecular system.
 * @param si the input SphInt object containing the spherical matrixelements
 */
void TPM::molecule(const SphInt &si){

   int a,b,c,d;
   int sign;

   int a_me,b_me,c_me,d_me;

   int S;

   double norm;

   for(int B = 0;B < gnr();++B){

      S = B2SM[B][0];

      sign = 1 - 2*S;

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         a_me = SphInt::gg2s(a);
         b_me = SphInt::gg2s(b);

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            c_me = SphInt::gg2s(c);
            d_me = SphInt::gg2s(d);

            //determine the norm for the basisset
            norm = 1.0;

            if(S == 0){

               if(a == b)
                  norm /= std::sqrt(2.0);

               if(c == d)
                  norm /= std::sqrt(2.0);

            }

            (*this)(B,i,j) = 0.0;

            if(b == d)
               (*this)(B,i,j) +=  1.0/(N - 1.0) * (si.gT(a_me,c_me) + si.gU(a_me,c_me));

            if(a == d)
               (*this)(B,i,j) +=  sign /(N - 1.0) * (si.gT(b_me,c_me) + si.gU(b_me,c_me));

            if(b == c)
               (*this)(B,i,j) +=  sign /(N - 1.0) * (si.gT(a_me,d_me) + si.gU(a_me,d_me));

            if(a == c)
               (*this)(B,i,j) +=  1.0/(N - 1.0) * (si.gT(b_me,d_me) + si.gU(b_me,d_me));

            //(*this)(B,i,j) += si.gV(a_me,b_me,c_me,d_me) + sign * si.gV(a_me,b_me,d_me,c_me);

            (*this)(B,i,j) *= norm;

         }
      }

   }

   this->symmetrize();

}

/**
 * initialize this onto the unitmatrix with trace N*(N - 1)/2
 */
void TPM::unit(){

   double ward = N*(N - 1.0)/(M*(M - 1.0));

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i){

         (*this)(B,i,i) = ward;

         for(int j = i + 1;j < gdim(B);++j)
            (*this)(B,i,j) = (*this)(B,j,i) = 0.0;

      }
   }

}

/**
 * orthogonal projection onto the space of traceless matrices
 */
void TPM::proj_Tr(){

   double ward = (2.0 * this->trace())/(M*(M - 1));

   this->min_unit(ward);

}

/**
 * Deduct the unitmatrix times a constant (scale) from this.\n\n
 * this -= scale* 1
 * @param scale the constant
 */

void TPM::min_unit(double scale){

   for(int B = 0;B < gnr();++B)
      for(int i = 0;i < gdim(B);++i)
         (*this)(B,i,i) -= scale;

}

/**
 * The spincoupled and axially symmetric Q map
 * @param option = 1, regular Q map , = -1 inverse Q map
 * @param tpm_d the TPM of which the Q map is taken and saved in this.
 */
void TPM::Q(int option,const TPM &tpm_d){

   double a = 1;
   double b = 1.0/(N*(N - 1.0));
   double c = 1.0/(N - 1.0);

   this->Q(option,a,b,c,tpm_d);

}

/**
 * The spincoupled and axially symmetric Q-like map: see primal-dual.pdf for more info (form: Q^S(A,B,C)(TPM) )
 * @param option = 1, regular Q-like map , = -1 inverse Q-like map
 * @param A factor in front of the two particle piece of the map
 * @param B factor in front of the no particle piece of the map
 * @param C factor in front of the single particle piece of the map
 * @param tpm_d the TPM of which the Q-like map is taken and saved in this.
 */
void TPM::Q(int option,double A,double B,double C,const TPM &tpm_d){

   //for inverse
   if(option == -1){

      B = (B*A + B*C*M - 2.0*C*C)/( A * (C*(M - 2.0) -  A) * ( A + B*M*(M - 1.0) - 2.0*C*(M - 1.0) ) );
      C = C/(A*(C*(M - 2.0) - A));
      A = 1.0/A;

   }

   SPM spm;
   spm.bar(C,tpm_d);

   //de trace*2 omdat mijn definitie van trace in berekeningen over alle (alpha,beta) loopt
   double ward = B*tpm_d.trace()*2.0;

   int ga,gb,gc,gd;

   int a,b,c,d;
   int m_a,m_b;

   int sign;

   double norm;

   for(int B = 0;B < gnr();++B){

      sign = 1 - 2*B2SM[B][0];

      for(int i = 0;i < gdim(B);++i){

         ga = t2s[B][i][0];
         gb = t2s[B][i][1];

         m_a = SPM::gg2ms(ga,0);
         a = SPM::gg2ms(ga,1);

         m_b = SPM::gg2ms(gb,0);
         b = SPM::gg2ms(gb,1);

         for(int j = i;j < gdim(B);++j){

            gc = t2s[B][j][0];
            gd = t2s[B][j][1];

            c = SPM::gg2ms(gc,1);
            d = SPM::gg2ms(gd,1);

            norm = 1.0;

            if(ga == gb)
               norm /= std::sqrt(2.0);

            if(gc == gd)
               norm /= std::sqrt(2.0);

            (*this)(B,i,j) = A * tpm_d(B,i,j);

            if(i == j)
               (*this)(B,i,j) += ward;

            if(ga == gc)
               (*this)(B,i,j) -= norm * spm(m_b + l_max,b,d);

            if(gb == gc)
               (*this)(B,i,j) -= norm * sign * spm(m_a + l_max,a,d);

            if(ga == gd)
               (*this)(B,i,j) -= norm * sign * spm(m_b + l_max,b,c);

            if(gb == gd)
               (*this)(B,i,j) -= norm * spm(m_a + l_max,a,c);

         }
      }

   }

   this->symmetrize();

}

/**
 * The G down map, maps a PHM object onto a TPM object using the G map.
 * @param phm input PHM
 */
void TPM::G(const PHM &phm){

   SPM spm;
   spm.bar(1.0/(N - 1.0),phm);

   int sign;

   int ga,gb,gc,gd;

   int a,b,c,d;
   int m_a,m_b,m_c,m_d;

   int S;

   for(int B = 0;B < gnr();++B){

      S = B2SM[B][0];

      sign = 1 - 2*S;

      for(int i = 0;i < gdim(B);++i){

         ga = t2s[B][i][0];
         gb = t2s[B][i][1];

         m_a = SPM::gg2ms(ga,0);
         a = SPM::gg2ms(ga,1);

         m_b = SPM::gg2ms(gb,0);
         b = SPM::gg2ms(gb,1);

         for(int j = i;j < gdim(B);++j){

            gc = t2s[B][j][0];
            gd = t2s[B][j][1];

            m_c = SPM::gg2ms(gc,0);
            c = SPM::gg2ms(gc,1);

            m_d = SPM::gg2ms(gd,0);
            d = SPM::gg2ms(gd,1);

            //init
            (*this)(B,i,j) = 0.0;

            //ph part
            for(int Z = 0;Z < 2;++Z){

               (*this)(B,i,j) -= (2*Z + 1.0) * Tools::g6j(0,0,S,Z) * ( phm(Z,m_a-m_d,m_a,a,-m_d,d,m_c,c,-m_b,b)

                     + phm(Z,m_b-m_c,m_b,b,-m_c,c,m_d,d,-m_a,a) + sign * phm(Z,m_b-m_d,m_b,b,-m_d,d,m_c,c,-m_a,a)

                     + sign * phm(Z,m_a-m_c,m_a,a,-m_c,c,m_d,d,-m_b,b) );

            }

            //4 sp parts
            if(gb == gd)
               (*this)(B,i,j) += spm(m_a+l_max,a,c);

            if(ga == gc)
               (*this)(B,i,j) += spm(m_b+l_max,b,d);

            if(ga == gd)
               (*this)(B,i,j) += sign * spm(m_b+l_max,b,c);

            if(gb == gc)
               (*this)(B,i,j) += sign * spm(m_a+l_max,a,d);

            //norm of the basisset:
            if(ga == gb)
               (*this)(B,i,j) /= std::sqrt(2.0);

            if(gc == gd)
               (*this)(B,i,j) /= std::sqrt(2.0);

         }

      }

   }

   this->symmetrize();

}

/**
 * Construct a spincoupled TPM matrix out of a spincoupled DPM matrix, for the definition and derivation see symmetry.pdf
 * @param dpm input DPM
 */
void TPM::bar(const DPM &dpm){

   int ga,gb,gc,gd;

   int a,b,c,d,l;
   int m_a,m_b,m_c,m_d,m_l;

   double ward;

   int S;

   for(int B = 0;B < gnr();++B){

      S = B2SM[B][0];

      if(S == 0){

         //first the S = 0 part, easiest:
         for(int i = 0;i < gdim(B);++i){

            ga = t2s[B][i][0];
            gb = t2s[B][i][1];

            m_a = SPM::gg2ms(ga,0);
            a = SPM::gg2ms(ga,1);

            m_b = SPM::gg2ms(gb,0);
            b = SPM::gg2ms(gb,1);

            for(int j = i;j < gdim(B);++j){

               gc = t2s[B][j][0];
               gd = t2s[B][j][1];

               m_c = SPM::gg2ms(gc,0);
               c = SPM::gg2ms(gc,1);

               m_d = SPM::gg2ms(gd,0);
               d = SPM::gg2ms(gd,1);

               (*this)(B,i,j) = 0.0;

               //only total S = 1/2 can remain because cannot couple to S = 3/2 with intermediate S = 0
               for(int gl = 0;gl < M/2;++gl){

                  m_l = SPM::gg2ms(gl,0);
                  l = SPM::gg2ms(gl,1);

                  (*this)(B,i,j) += 2.0 * dpm(0,m_a+m_b+m_l,0,m_a,a,m_b,b,m_l,l,0,m_c,c,m_d,d,m_l,l);

               }

            }
         }

      }
      else{

         //then the S = 1 part:
         for(int i = 0;i < gdim(B);++i){

            ga = t2s[B][i][0];
            gb = t2s[B][i][1];

            m_a = SPM::gg2ms(ga,0);
            a = SPM::gg2ms(ga,1);

            m_b = SPM::gg2ms(gb,0);
            b = SPM::gg2ms(gb,1);

            for(int j = i;j < gdim(B);++j){

               gc = t2s[B][j][0];
               gd = t2s[B][j][1];

               m_c = SPM::gg2ms(gc,0);
               c = SPM::gg2ms(gc,1);

               m_d = SPM::gg2ms(gd,0);
               d = SPM::gg2ms(gd,1);

               (*this)(B,i,j) = 0.0;

               for(int Z = 0;Z < 2;++Z){//loop over the dpm blocks: S = 1/2 and 3/2 = Z + 1/2

                  ward = 0.0;

                  for(int gl = 0;gl < M/2;++gl){

                     m_l = SPM::gg2ms(gl,0);
                     l = SPM::gg2ms(gl,1);

                     ward += dpm(Z,m_a+m_b+m_l,1,m_a,a,m_b,b,m_l,l,1,m_c,c,m_d,d,m_l,l);

                  }

                  ward *= (2 * (Z + 0.5) + 1.0)/3.0;

                  (*this)(B,i,j) += ward;

               }

            }
         }

      }

   }

   this->symmetrize();

}

/** 
 * The T1-down map that maps a DPM on TPM. This is just a Q-like map using the TPM::bar (dpm) as input.
 * @param dpm the input DPM matrix
 */
void TPM::T(const DPM &dpm){

   TPM tpm;
   tpm.bar(dpm);

   double a = 1;
   double b = 1.0/(3.0*N*(N - 1.0));
   double c = 0.5/(N - 1.0);

   this->Q(1,a,b,c,tpm);

}

/**
 * print the TPM object in a non-axially symmetric form
 */
void TPM::printnax(const char *filename) const {

   ofstream out(filename);
   out.precision(15);

   int S;

   int ga,gb,gc,gd;

   for(int B = 0;B < gnr();++B){

      S = B2SM[B][0];

      for(int i = 0;i < gdim(B);++i){

         ga = t2s[B][i][0];
         gb = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            gc = t2s[B][j][0];
            gd = t2s[B][j][1];

            out << S << "\t" << ga << "\t" << gb << "\t" << gc << "\t" << gd << "\t" << (*this)(B,i,j) << endl;

         }
      }

   }
}

/**
 * The bar function that maps a PPHM object onto a TPM object by tracing away the last pair of incdices of the PPHM
 * @param pphm Input PPHM object
 */
void TPM::bar(const PPHM &pphm){

   int ga,gb,gc,gd;

   int a,b,c,d,l;
   int m_a,m_b,m_c,m_d,m_l;

   double ward;

   int Z;

   for(int B = 0;B < gnr();++B){

      Z = B2SM[B][0];

      for(int i = 0;i < gdim(B);++i){

         ga = t2s[B][i][0];
         gb = t2s[B][i][1];

         m_a = SPM::gg2ms(ga,0);
         a = SPM::gg2ms(ga,1);

         m_b = SPM::gg2ms(gb,0);
         b = SPM::gg2ms(gb,1);

         for(int j = i;j < gdim(B);++j){

            gc = t2s[B][j][0];
            gd = t2s[B][j][1];

            m_c = SPM::gg2ms(gc,0);
            c = SPM::gg2ms(gc,1);

            m_d = SPM::gg2ms(gd,0);
            d = SPM::gg2ms(gd,1);

            (*this)(B,i,j) = 0.0;

            for(int S = 0;S < 2;++S){//loop over three particle spin: 1/2 and 3/2

               ward = (2.0*(S + 0.5) + 1.0)/(2.0*Z + 1.0);

               for(int gl = 0;gl < M/2;++gl){

                  m_l = SPM::gg2ms(gl,0);
                  l = SPM::gg2ms(gl,1);

                  (*this)(B,i,j) += ward * pphm(S,m_a + m_b + m_l,Z,m_a,a,m_b,b,m_l,l,Z,m_c,c,m_d,d,m_l,l);

               }

            }

         }
      }

   }

   this->symmetrize();

}


/**
 * The spincoupled and axially symmetric T2-down map that maps a PPHM on a TPM object.
 * @param pphm input PPHM object
 */
void TPM::T(const PPHM &pphm){

   //first make the bar tpm
   TPM tpm;
   tpm.bar(pphm);

   //then make the bar phm
   PHM phm;
   phm.bar(pphm);

   //also make the bar spm with the correct scale factor
   SPM spm;
   spm.bar(0.5/(N - 1.0),pphm);

   int ga,gb,gc,gd;

   int a,b,c,d;
   int m_a,m_b,m_c,m_d;

   int S;

   int sign;

   double norm;

   for(int B = 0;B < gnr();++B){

      S = B2SM[B][0];

      sign = 1 - 2*S;

      for(int i = 0;i < gdim(B);++i){

         ga = t2s[B][i][0];
         gb = t2s[B][i][1];

         m_a = SPM::gg2ms(ga,0);
         a = SPM::gg2ms(ga,1);

         m_b = SPM::gg2ms(gb,0);
         b = SPM::gg2ms(gb,1);

         for(int j = i;j < gdim(B);++j){

            gc = t2s[B][j][0];
            gd = t2s[B][j][1];

            m_c = SPM::gg2ms(gc,0);
            c = SPM::gg2ms(gc,1);

            m_d = SPM::gg2ms(gd,0);
            d = SPM::gg2ms(gd,1);

            //determine the norm for the basisset
            norm = 1.0;

            if(S == 0){

               if(ga == gb)
                  norm /= std::sqrt(2.0);

               if(gc == gd)
                  norm /= std::sqrt(2.0);

            }

            //first the tp part
            (*this)(B,i,j) = tpm(B,i,j);

            //sp part, 4 terms:
            if(gb == gd)
               (*this)(B,i,j) += norm * spm(m_a+l_max,a,c);

            if(ga == gd)
               (*this)(B,i,j) += sign * norm * spm(m_b+l_max,b,c);

            if(gb == gc)
               (*this)(B,i,j) += sign * norm * spm(m_a+l_max,a,d);

            if(ga == gc)
               (*this)(B,i,j) += norm * spm(m_b+l_max,b,d);

            for(int Z = 0;Z < 2;++Z){

               (*this)(B,i,j) -= norm * (2.0 * Z + 1.0) * Tools::g6j(0,0,S,Z) * ( phm(Z,m_d-m_a,m_d,d,-m_a,a,m_b,b,-m_c,c) 

                     + sign * phm(Z,m_d-m_b,m_d,d,-m_b,b,m_a,a,-m_c,c) + sign * phm(Z,m_c-m_a,m_c,c,-m_a,a,m_b,b,-m_d,d)

                     + phm(Z,m_c-m_b,m_c,c,-m_b,b,m_a,a,-m_d,d) );

            }

         }
      }

   }

   this->symmetrize();

}

/**
 * Collaps a SUP matrix S onto a TPM matrix like this:\n\n
 * sum_i Tr (S u^i)f^i = this
 * @param option = 0, project onto full symmetric matrix space, = 1 project onto traceless symmetric matrix space
 * @param S input SUP
 */
void TPM::collaps(int option,const SUP &S){

   *this = S.gI();

#ifdef __Q_CON
   TPM hulp;

   hulp.Q(1,S.gQ());

   *this += hulp;
#endif

#ifdef __G_CON
   hulp.G(S.gG());

   *this += hulp;
#endif

#ifdef __T1_CON
   hulp.T(S.gT1());

   *this += hulp;
#endif

#ifdef __T2_CON
   hulp.T(S.gT2());

   *this += hulp;
#endif

   if(option == 1)
      this->proj_Tr();

}

/**
 * ( Overlapmatrix of the U-basis ) - map, maps a TPM onto a different TPM, this map is actually a Q-like map
 * for which the paramaters a,b and c are calculated in primal_dual.pdf. Since it is a Q-like map the inverse
 * can be taken as well.
 * @param option = 1 direct overlapmatrix-map is used , = -1 inverse overlapmatrix map is used
 * @param tpm_d the input TPM
 */

void TPM::S(int option,const TPM &tpm_d){

   this->Q(option,Sa,0.0,Sc,tpm_d);

}

int TPM::SaveToFile(const char *filename)
{
   hid_t       file_id, group_id, dataset_id, attribute_id, dataspace_id, strtype;
   hsize_t     dims = 1;
   herr_t      status;

   // new file
   file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

   // make group for rdm
   group_id = H5Gcreate(file_id, "/RDM", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   for(int i=0;i<this->gnr();i++)
   {
      hsize_t dimblock = this->gdim(i)*this->gdim(i);

      // make dataspace for a block
      dataspace_id = H5Screate_simple(1, &dimblock, NULL);

      char name[16];
      sprintf(name,"/RDM/%d",i);

      // make dataset
      dataset_id = H5Dcreate(file_id, name, H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      double **matrix = (*this)[i].gMatrix();

      // fill dataset
      status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix[0] );
      HDF5_STATUS_CHECK(status);

      /* Terminate access to the data space. */
      status = H5Sclose(dataspace_id);
      HDF5_STATUS_CHECK(status);

      // add as attribute the S and M quantum numbers to each block
      dataspace_id = H5Screate_simple(1, &dims, NULL);

      attribute_id = H5Acreate (dataset_id, "S", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite (attribute_id, H5T_NATIVE_INT, &B2SM[i][0] );
      HDF5_STATUS_CHECK(status);

      status = H5Aclose(attribute_id);
      HDF5_STATUS_CHECK(status);

      attribute_id = H5Acreate (dataset_id, "M", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite (attribute_id, H5T_NATIVE_INT, &B2SM[i][1] );
      HDF5_STATUS_CHECK(status);

      status = H5Aclose(attribute_id);
      HDF5_STATUS_CHECK(status);
      status = H5Sclose(dataspace_id);
      HDF5_STATUS_CHECK(status);

      /* End access to the dataset and release resources used by it. */
      status = H5Dclose(dataset_id);
      HDF5_STATUS_CHECK(status);

   }

   int nr_blocks = this->gnr();
   dataspace_id = H5Screate_simple(1, &dims, NULL);

   attribute_id = H5Acreate (group_id, "Blocks", H5T_STD_I64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite (attribute_id, H5T_NATIVE_INT, &nr_blocks );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   attribute_id = H5Acreate (group_id, "M", H5T_STD_I64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite (attribute_id, H5T_NATIVE_INT, &M );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   attribute_id = H5Acreate (group_id, "N", H5T_STD_I64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite (attribute_id, H5T_NATIVE_INT, &N );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   std::ifstream specs("start.stp");
   std::stringstream specs_buffer;
   specs_buffer << specs.rdbuf();
   std::string specs_string = specs_buffer.str();

   // make string type of correct size
   strtype = H5Tcopy(H5T_C_S1);
   H5Tset_size(strtype,specs_string.size());

   attribute_id = H5Acreate (group_id, "start.stp", strtype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite (attribute_id, strtype, specs_string.c_str() );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   status = H5Tclose(strtype);
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);

   dims = 6;
   dataspace_id = H5Screate_simple(1, &dims, NULL);

   int typeofcalculation[6];

   typeofcalculation[0] = 1; // P

#ifdef __Q_CON
   typeofcalculation[1] = 1; // Q
#else
   typeofcalculation[1] = 0; // Q
#endif

#ifdef __G_CON
   typeofcalculation[2] = 1; // G
#else
   typeofcalculation[2] = 0; // Q
#endif

#ifdef __T1_CON
   typeofcalculation[3] = 1; // T1
#else
   typeofcalculation[3] = 0; // Q
#endif

#ifdef __T2_CON
   typeofcalculation[4] = 1; // T2
#else
   typeofcalculation[4] = 0; // Q
#endif

#ifdef __T2P_CON
   typeofcalculation[5] = 1; // T2P
#else
   typeofcalculation[5] = 0; // Q
#endif

   attribute_id = H5Acreate (group_id, "Type", H5T_STD_I64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite (attribute_id, H5T_NATIVE_INT, &typeofcalculation[0] );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);

   /* Close the group. */
   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);

   /* Terminate access to the file. */
   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);

   return 0;
}

int TPM::ReadInitfromFile(const char *filename, string &setupdata)
{
   hid_t file_id, group_id, attribute_id, strtype, dataspace_id;
   herr_t status;
   size_t sdim;
   hsize_t dims = 6;

   // open file
   file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
   HDF5_STATUS_CHECK(file_id);

   // open group for rdm
   group_id = H5Gopen(file_id, "/RDM", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   status = H5Aexists(group_id,"start.stp");
   HDF5_STATUS_CHECK(status);

   if(status == 0)
   {
      std::cerr << "HDF5 input file '" << filename << "' is not of the correct type." << std::endl;
      return -1;
   }

   attribute_id = H5Aopen(group_id, "start.stp", H5P_DEFAULT);
   HDF5_STATUS_CHECK(attribute_id);

   strtype = H5Aget_type(attribute_id);
   sdim = H5Tget_size(strtype);

   char *specs = new char[sdim+1];
   specs[sdim] = '\0';

   status = H5Aread(attribute_id, strtype, specs);
   HDF5_STATUS_CHECK(status);

   setupdata = specs;

   delete [] specs;

   status = H5Tclose(strtype);
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   status = H5Aexists(group_id,"Type");
   HDF5_STATUS_CHECK(status);

   if(status == 0)
      std::cerr << "HDF5 input file '" << filename << "' has no information about the type of calculation (P,Q,G,T1 or T2)." << std::endl;
   else
   {
      dataspace_id = H5Screate_simple(1, &dims, NULL);

      int typeofcalculation[6];

      attribute_id = H5Aopen(group_id, "Type", H5P_DEFAULT);
      HDF5_STATUS_CHECK(attribute_id);

      status = H5Aread(attribute_id, H5T_NATIVE_INT, &typeofcalculation[0]);
      HDF5_STATUS_CHECK(status);

      status = H5Aclose(attribute_id);
      HDF5_STATUS_CHECK(status);

      status = H5Sclose(dataspace_id);
      HDF5_STATUS_CHECK(status);

#ifdef __Q_CON
      if(typeofcalculation[1] == 0)
         std::cerr << "HDF5 input file '" << filename << "' hasn't got the Q condition active while the program has" << std::endl;
#else
      if(typeofcalculation[1] == 1)
         std::cerr << "HDF5 input file '" << filename << "' has got the Q condition active while the program has not" << std::endl;
#endif

#ifdef __G_CON
      if(typeofcalculation[2] == 0)
         std::cerr << "HDF5 input file '" << filename << "' hasn't got the G condition active while the program has" << std::endl;
#else
      if(typeofcalculation[2] == 1)
         std::cerr << "HDF5 input file '" << filename << "' has got the G condition active while the program has not" << std::endl;
#endif

#ifdef __T1_CON
      if(typeofcalculation[3] == 0)
         std::cerr << "HDF5 input file '" << filename << "' hasn't got the T1 condition active while the program has" << std::endl;
#else
      if(typeofcalculation[3] == 1)
         std::cerr << "HDF5 input file '" << filename << "' has got the T1 condition active while the program has not" << std::endl;
#endif

#ifdef __T2_CON
      if(typeofcalculation[4] == 0)
         std::cerr << "HDF5 input file '" << filename << "' hasn't got the T2 condition active while the program has" << std::endl;
#else
      if(typeofcalculation[4] == 1)
         std::cerr << "HDF5 input file '" << filename << "' has got the T2 condition active while the program has not" << std::endl;
#endif

#ifdef __T2P_CON
      if(typeofcalculation[5] == 0)
         std::cerr << "HDF5 input file '" << filename << "' hasn't got the T2P condition active while the program has" << std::endl;
#else
      if(typeofcalculation[5] == 1)
         std::cerr << "HDF5 input file '" << filename << "' has got the T2P condition active while the program has not" << std::endl;
#endif

   }

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);

   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);

   return 0;
}

int TPM::ReadfromFile(const char *filename)
{
   hid_t       file_id, group_id, dataset_id;
   herr_t      status;

   // open file
   file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
   HDF5_STATUS_CHECK(file_id);

   // open group for rdm
   group_id = H5Gopen(file_id, "/RDM", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   for(int i=0;i<this->gnr();i++)
   {
      char name[16];
      sprintf(name,"/RDM/%d",i);

      // make dataset
      dataset_id = H5Dopen(group_id, name, H5P_DEFAULT);
      HDF5_STATUS_CHECK(dataset_id);

      double **matrix = (*this)[i].gMatrix();

      // fill dataset
      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix[0]);
      HDF5_STATUS_CHECK(status);

      status = H5Dclose(dataset_id);
      HDF5_STATUS_CHECK(status);
   }

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);

   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);

   return 0;
}

/* vim: set ts=3 sw=3 expandtab :*/
