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

vector< vector< vector<int> > > DPM::d2s;
int *****DPM::s2d;

vector< vector<int> > DPM::B2SM;
int **DPM::SM2B;

int DPM::M;
int DPM::N;

int DPM::l_max;

/**
 * initialize the static lists and variables
 * @param M_in input dimension of sp space
 * @param N_in input nr of particles
 */
void DPM::init(int M_in,int N_in){

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

      //first S_ab = 0 and ga == gb != gc
      for(int ga = 0;ga < M/2;++ga){

         for(int gc = 0;gc < ga;++gc){

            if(L_z == 2*SPM::gg2ms(ga,0) + SPM::gg2ms(gc,0)){//correct L_z projection

               v[0] = 0;//S_ab
               v[1] = ga;
               v[2] = ga;
               v[3] = gc;

               bv.push_back(v);

            }

         }

         for(int gc = ga + 1;gc < M/2;++gc){

            if(L_z == 2*SPM::gg2ms(ga,0) + SPM::gg2ms(gc,0)){//correct L_z projection

               v[0] = 0;//S_ab
               v[1] = ga;
               v[2] = ga;
               v[3] = gc;

               bv.push_back(v);

            }

         }

      }

      //then the rest, i.e. S_ab = 0/1 ga < gb < gc
      for(int S_ab = 0;S_ab < 2;++S_ab)
         for(int ga = 0;ga < M/2;++ga)
            for(int gb = ga + 1;gb < M/2;++gb)
               for(int gc = gb + 1;gc < M/2;++gc){

                  if(L_z == SPM::gg2ms(ga,0) + SPM::gg2ms(gb,0) + SPM::gg2ms(gc,0)){//correct L_z projection

                     v[0] = S_ab;
                     v[1] = ga;
                     v[2] = gb;
                     v[3] = gc;

                     bv.push_back(v);

                  }

               }

      if(bv.size() != 0){

         d2s.push_back(bv);

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
            for(int gc = gb + 1;gc < M/2;++gc){

               if(L_z == SPM::gg2ms(ga,0) + SPM::gg2ms(gb,0) + SPM::gg2ms(gc,0)){//correct L_z projection

                  v[0] = 1;
                  v[1] = ga;
                  v[2] = gb;
                  v[3] = gc;

                  bv.push_back(v);

               }

            }

      if(bv.size() != 0){

         d2s.push_back(bv);

         bq[0] = 1;
         bq[1] = L_z;

         B2SM.push_back(bq);

         SM2B[1][L_z + 3*l_max] = B;

         ++B;

      }

   }

   s2d = new int **** [B2SM.size()];

   for(unsigned int B = 0;B < B2SM.size();++B){

      s2d[B] = new int *** [2];

      for(int S_ab = 0;S_ab < 2;++S_ab){

         s2d[B][S_ab] = new int ** [M/2];

         for(int ga = 0;ga < M/2;++ga){

            s2d[B][S_ab][ga] = new int * [M/2];

            for(int gb = 0;gb < M/2;++gb)
               s2d[B][S_ab][ga][gb] = new int [M/2];

         }
      }
   }

   for(unsigned int B = 0;B < B2SM.size();++B)
      for(unsigned int d = 0;d < d2s[B].size();++d)
         s2d[B][d2s[B][d][0]][d2s[B][d][1]][d2s[B][d][2]][d2s[B][d][3]] = d;

}

/**
 * deallocate the static lists
 */
void DPM::clear(){

   for(unsigned int B = 0;B < B2SM.size();++B){

      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int ga = 0;ga < M/2;++ga){

            for(int gb = 0;gb < M/2;++gb)
               delete [] s2d[B][S_ab][ga][gb];

            delete [] s2d[B][S_ab][ga];

         }

         delete [] s2d[B][S_ab];

      }

      delete [] s2d[B];

   }

   delete [] s2d;

   for(int S = 0;S < 2;++S)
      delete [] SM2B[S];

   delete [] SM2B;

}

/**
 * standard constructor for a spin and axial symmetrical tp matrix.
 */
DPM::DPM() : BlockMatrix(B2SM.size()) {

   for(int B = 0;B < gnr();++B)
      this->setMatrixDim(B,d2s[B].size(),2*B2SM[B][0] + 2);

}

/**
 * copy constructor:
 * @param dpm_c object that will be copied into this.
 */
DPM::DPM(const DPM &dpm_c) : BlockMatrix(dpm_c){ }

/**
 * destructor
 */
DPM::~DPM(){ }

/**
 * @return number of particles
 */
int DPM::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int DPM::gM() const{

   return M;

}

ostream &operator<<(ostream &output,const DPM &dpm_p){

   for(int B = 0;B < dpm_p.gnr();++B){

      output << endl;
      output << "Block " << B << " with S = " << dpm_p.B2SM[B][0] << " and M = " << dpm_p.B2SM[B][1] << endl;
      output << endl;

      int ga,gb,gc,gd,ge,gz;
      int S_ab,S_de;

      for(int d_i = 0;d_i < dpm_p.gdim(B);++d_i){

         S_ab = dpm_p.d2s[B][d_i][0];

         ga = dpm_p.d2s[B][d_i][1];
         gb = dpm_p.d2s[B][d_i][2];
         gc = dpm_p.d2s[B][d_i][3];

         for(int d_j = d_i;d_j < dpm_p.gdim(B);++d_j){

            S_de = dpm_p.d2s[B][d_j][0];

            gd = dpm_p.d2s[B][d_j][1];
            ge = dpm_p.d2s[B][d_j][2];
            gz = dpm_p.d2s[B][d_j][3];

            output << "[ " << S_ab << "\t(" << SPM::gg2ms(ga,0) << "," << SPM::gg2ms(ga,1) << ")\t"

               << "(" << SPM::gg2ms(gb,0) << "," << SPM::gg2ms(gb,1) << ")\t(" <<  SPM::gg2ms(gc,0) << "," << SPM::gg2ms(gc,1) << ") ]"

               << "\t|\t[ " << S_de << "\t(" << SPM::gg2ms(gd,0) << "," << SPM::gg2ms(gd,1) << ")\t"

               << "(" << SPM::gg2ms(ge,0) << "," << SPM::gg2ms(ge,1) << ")\t(" << SPM::gg2ms(gz,0) << "," << SPM::gg2ms(gz,1) << ") ]"

               << "\t||\t" << dpm_p(B,d_i,d_j) << endl;

         }
      }

   }

   return output;

}

/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * DPM(S,S_ab,a,b,c,S_de,d,e,f) = sum_S_ac (some terms dependent on spin) DPM(S,S_ac,a,c,b,S_de,d,e,f) etc...
 * @param S The block index, when == 0 then access the block S = 1/2, for block == 1 we access the S = 3/2.
 * @param Lz the three particle angular momentum z-projection index
 * @param S_ab The intermediate spinquantumnumber of a and b.
 * @param m_a z-projectino of the first index
 * @param a first sp index that forms the dp row index i together with b, c and S_ab in block S
 * @param m_b z-projection of the second index
 * @param b second sp index that forms the dp row index i together with a, c and S_ab in block S
 * @param m_c z-projection of the third index
 * @param c third sp index that forms the dp row index i together with a, b and S_ab in block S
 * @param S_de The intermediate spinquantumnumber of d and e.
 * @param m_d z-projection of the first column index
 * @param d first sp index that forms the dp column index j together with e, z and S_de in block S
 * @param m_e z-projection of the second column index
 * @param e second sp index that forms the dp column index j together with d, z and S_de in block S
 * @param m_z z-projection of the third column index
 * @param z third sp index that forms the dp column index j together with d, e and S_de in block S
 * @return the number on place DPM(S,i,j) with the right phase and forefactor.
 */
double DPM::operator()(int S,int Lz,int S_ab,int m_a,int a,int m_b,int b,int m_c,int c,int S_de,int m_d,int d,int m_e,int e,int m_z,int z) const {

   if(m_a + m_b + m_c != m_d + m_e + m_z)
      return 0.0;

   if(Lz != m_a + m_b + m_c)
      return 0.0;

   int ga = SPM::gms2g(m_a,a);
   int gb = SPM::gms2g(m_b,b);
   int gc = SPM::gms2g(m_c,c);

   int *i = new int [2];
   double *coef_i = new double [2];

   int dim_i = get_inco(S,Lz,S_ab,ga,gb,gc,i,coef_i);

   if(dim_i == 0){

      delete [] i;
      delete [] coef_i;

      return 0.0;

   }

   int gd = SPM::gms2g(m_d,d);
   int ge = SPM::gms2g(m_e,e);
   int gz = SPM::gms2g(m_z,z);

   int *j = new int [2];
   double *coef_j = new double [2];

   int dim_j = get_inco(S,Lz,S_de,gd,ge,gz,j,coef_j);

   if(dim_j == 0){

      delete [] i;
      delete [] j;

      delete [] coef_i;
      delete [] coef_j;

      return 0.0;

   }

   double ward = 0.0;

   for(int I = 0;I < dim_i;++I)
      for(int J = 0;J < dim_j;++J)
         ward += coef_i[I] * coef_j[J] * (*this)(SM2B[S][Lz + 3*l_max],i[I],j[J]);

   delete [] i;
   delete [] j;

   delete [] coef_i;
   delete [] coef_j;

   return ward;

}

/** 
 * Static member function that gets the dp-indices and their coefficients in terms of S,S_ab,ga,gb,gc.
 * @param S three-particle spin
 * @param Lz three-particle momentum projection
 * @param S_ab intermediate spincoupling of a and b.
 * @param ga first sp orbital
 * @param gb second sp orbital
 * @param gc third sp orbital
 * @param i pointer of dim 1 or 2 containing the indices occuring in the expansion of this particular dp state in the normal basis (a==b,c a < b < c).
 * @param coef pointer of dim 1 or 2 containing the coefficients occuring in the expansion.
 * @return the number of terms in the expansion (1 or 2), also the dim of pointers i and coef. When zero is returned this is not a valid element.
 */
int DPM::get_inco(int S,int Lz,int S_ab,int ga,int gb,int gc,int *i,double *coef) const{

   //they cannot all be equal
   if(ga == gb && gb == gc)
      return 0;

   if(S == 0){//spin 1/2 block:

      //if normal basis:
      if(ga == gb){

         if(S_ab == 1)//spin has to be zero for a == b
            return 0;

         i[0] = s2d[SM2B[S][Lz + 3*l_max]][0][ga][gb][gc];
         coef[0] = 1;

         return 1;

      }
      else if (ga < gb && gb < gc){

         i[0] = s2d[SM2B[S][Lz + 3*l_max]][S_ab][ga][gb][gc];
         coef[0] = 1;

         return 1;

      }
      else{//anomal basis:

         int min,max,phase;

         //first order a and b for code saving reasons
         if(ga < gb){

            min = ga;
            max = gb;

            phase = 1;

         }
         else{

            min = gb;
            max = ga;

            phase = 1 - 2*S_ab;

            if(gc > max){//we still have one simple dim = 1 term left: b < a < c

               i[0] = s2d[SM2B[S][Lz + 3*l_max]][S_ab][gb][ga][gc];
               coef[0] = phase;

               return 1;

            }

         }

         //now we have four possibilities left:
         //don't forget to multiply every result by phase to get the right a and b for min and max!
         // 1) gc < min < max
         // 2) gc == min < max
         // 3) min < gc < max
         // 4) min < max == gc
         if(gc < min){//gc < min < max

            //the S_ca == 0 part:
            i[0] = s2d[SM2B[S][Lz + 3*l_max]][0][gc][min][max];
            coef[0] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * Tools::g6j(0,0,0,S_ab);

            //the S_ca == 1 part:
            i[1] = s2d[SM2B[S][Lz + 3*l_max]][1][gc][min][max];
            coef[1] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * std::sqrt(3.0) * Tools::g6j(0,0,1,S_ab);

            return 2;

         }
         else if(gc == min){//gc == min < max: this will also be a 1 dim list, because S_ac can only be 0 if ga == gc.

            i[0] = s2d[SM2B[S][Lz + 3*l_max]][0][gc][min][max];
            coef[0] = std::sqrt(2.0) * phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * Tools::g6j(0,0,0,S_ab);

            return 1;

         }
         else if(gc < max){//min < gc < max

            //S_ac == 0 part:
            i[0] = s2d[SM2B[S][Lz + 3*l_max]][0][min][gc][max];
            coef[0] = phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * Tools::g6j(0,0,0,S_ab);

            //S_ac == 1 part:
            i[1] = s2d[SM2B[S][Lz + 3*l_max]][1][min][gc][max];
            coef[1] = - phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * std::sqrt(3.0) * Tools::g6j(0,0,1,S_ab);

            return 2;

         }
         else{// min < gc == max: also a 1 dim list, s_bc can only be 0 if gb == gc

            i[0] = s2d[SM2B[S][Lz + 3*l_max]][0][max][gc][min];
            coef[0] = phase * std::sqrt(2.0) * std::sqrt(2.0*S_ab + 1.0) * Tools::g6j(0,0,0,S_ab);

            return 1;

         }

      }

   }
   else{//spin 3/2 block, totally antisymmetrical in the spatial sp orbs.

      //only S_ab == 1 can couple to 3/2's.
      if(S_ab == 0)
         return 0;

      //if any of the sp orbs are equal, antisymmetry leads to zero:
      if(ga == gb || gb == gc || gc == ga)
         return 0;

      if(ga < gb){

         if(gb < gc){//ga < gb < gc

            i[0] = s2d[SM2B[S][Lz + 3*l_max]][1][ga][gb][gc];
            coef[0] = 1;

         }
         else if(gc < ga){//gc < ga < gb

            i[0] = s2d[SM2B[S][Lz + 3*l_max]][1][gc][ga][gb];
            coef[0] = 1;

         }
         else{//ga < gc < gb

            i[0] = s2d[SM2B[S][Lz + 3*l_max]][1][ga][gc][gb];
            coef[0] = -1;

         }

      }
      else{//gb < ga

         if(ga < gc){//gb < ga < gc

            i[0] = s2d[SM2B[S][Lz + 3*l_max]][1][gb][ga][gc];
            coef[0] = -1;

         }
         else if(gc < gb){//gc < gb < ga

            i[0] = s2d[SM2B[S][Lz + 3*l_max]][1][gc][gb][ga];
            coef[0] = -1;

         }
         else{//gb < gc < ga

            i[0] = s2d[SM2B[S][Lz + 3*l_max]][1][gb][gc][ga];
            coef[0] = 1;

         }

      }

      return 1;

   }

}

/**
 * The spincoupled T1-like (generalized T1) map: maps a TPM object (tpm) on a DPM object (*this)
 * @param A term before the tp part of the map
 * @param B term before the np part of the map
 * @param C term before the sp part of the map
 * @param tpm input TPM
 */
void DPM::T(double A,double B,double C,const TPM &tpm){

   //make sp matrix out of tpm
   SPM spm;
   spm.bar(C,tpm);

   double ward = 2.0*B*tpm.trace();

   int ga,gb,gc,gd,ge,gz;
   int S_ab,S_de;

   int a,b,c,d,e,z;
   int m_a,m_b,m_c,m_d,m_e,m_z;

   int S;

   int sign_ab,sign_de;

   double norm_ab,norm_de;

   double hard;

   for(int B = 0;B < gnr();++B){

      S = B2SM[B][0];

      //start with the S = 1/2 block, this is the most difficult one:
      if(S == 0){

         for(int i = 0;i < gdim(B);++i){

            S_ab = d2s[B][i][0];

            ga = d2s[B][i][1];
            gb = d2s[B][i][2];
            gc = d2s[B][i][3];

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

               S_de = d2s[B][j][0];

               gd = d2s[B][j][1];
               ge = d2s[B][j][2];
               gz = d2s[B][j][3];

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

               hard = std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * Tools::g6j(0,0,S_ab,S_de);

               //init
               (*this)(B,i,j) = 0.0;

               //the np part
               if(i == j)
                  (*this)(B,i,j) = ward;

               //other parts are a bit more difficult.
               if(gc == gz){

                  if(S_ab == S_de){

                     //tp(1)
                     (*this)(B,i,j) += A * tpm(S_ab,m_a + m_b,m_a,a,m_b,b,m_d,d,m_e,e);

                     //sp(1) first term
                     if(gb == ge)
                        (*this)(B,i,j) -= norm_ab * norm_de * spm(m_a+l_max,a,d);

                     //sp(2) first term
                     if(ga == ge)
                        (*this)(B,i,j) -= sign_ab * norm_ab * norm_de * spm(m_b+l_max,b,d);

                     //sp(4) first term
                     if(gb == gd)
                        (*this)(B,i,j) -= sign_de * norm_ab * norm_de * spm(m_a+l_max,a,e);

                     //sp(5) first term
                     if(ga == gd)
                        (*this)(B,i,j) -= norm_ab * norm_de * spm(m_b+l_max,b,e);

                  }

               }

               if(gb == gz){

                  //tp(2)
                  if(ga == gc)
                     (*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,m_a+m_c,m_a,a,m_c,c,m_d,d,m_e,e);
                  else
                     (*this)(B,i,j) += A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,m_a+m_c,m_a,a,m_c,c,m_d,d,m_e,e);

                  //sp(1) second term
                  if(gc == ge)
                     (*this)(B,i,j) -= sign_ab * sign_de * norm_ab * norm_de * hard * spm(m_a+l_max,a,d);

                  //sp(3)
                  if(ga == ge)
                     (*this)(B,i,j) -= sign_ab * norm_ab * norm_de * hard * spm(m_c+l_max,c,d);

                  //sp(4) second term
                  if(gc == gd)
                     (*this)(B,i,j) -= sign_ab * norm_ab * norm_de * hard * spm(m_a+l_max,a,e);

                  //sp(6)
                  if(ga == gd)
                     (*this)(B,i,j) -= sign_ab * sign_de * norm_ab * norm_de * hard * spm(m_c+l_max,c,e);

               }

               if(ga == gz){

                  //tp(3)
                  if(gb == gc)
                     (*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_de * hard * tpm(S_de,m_b+m_c,m_b,b,m_c,c,m_d,d,m_e,e);
                  else
                     (*this)(B,i,j) += A * norm_ab * sign_de * hard * tpm(S_de,m_b+m_c,m_b,b,m_c,c,m_d,d,m_e,e);

                  //sp(2) second term
                  if(gc == ge)
                     (*this)(B,i,j) -= sign_de * norm_ab * norm_de * hard * spm(m_b+l_max,b,d);

                  //sp(5) second term
                  if(gc == gd)
                     (*this)(B,i,j) -= norm_ab * norm_de * hard * spm(m_b+l_max,b,e);

               }

               if(gc == ge){

                  //tp(4)
                  if(gd == gz)
                     (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,m_a+m_b,m_a,a,m_b,b,m_d,d,m_z,z);
                  else
                     (*this)(B,i,j) += A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,m_a+m_b,m_a,a,m_b,b,m_d,d,m_z,z);

                  //sp(7) first term
                  if(gb == gd)
                     (*this)(B,i,j) -= norm_ab * norm_de * sign_de * hard * spm(m_a+l_max,a,z);

                  //sp(8) first term
                  if(ga == gd)
                     (*this)(B,i,j) -= norm_ab * norm_de * sign_ab * sign_de * hard * spm(m_b+l_max,b,z);

               }

               if(gb == ge){

                  //tp(5)
                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,m_a+m_c,m_a,a,m_c,c,m_d,d,m_z,z);

                  //correct for norms of the tpm
                  if(ga == gc)
                     hulp *= std::sqrt(2.0);

                  if(gd == gz)
                     hulp *= std::sqrt(2.0);

                  (*this)(B,i,j) += A * norm_ab * norm_de * sign_ab * sign_de * std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * hulp;

                  //sp(7) second term
                  if(gc == gd)
                     (*this)(B,i,j) -= norm_ab * norm_de * hard * spm(m_a+l_max,a,z);

                  //sp(9) first term
                  if(ga == gd)
                     if(S_ab == S_de)
                        (*this)(B,i,j) -= norm_ab * norm_de * spm(m_c+l_max,c,z);

               }

               if(ga == ge){

                  //tp(6)
                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,m_b+m_c,m_b,b,m_c,c,m_d,d,m_z,z);

                  if(gb == gc)
                     hulp *= std::sqrt(2.0);

                  if(gd == gz)
                     hulp *= std::sqrt(2.0);

                  (*this)(B,i,j) += A * sign_de * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

                  //sp(8) second term
                  if(gc == gd)
                     (*this)(B,i,j) -= sign_ab * norm_ab * norm_de * hard * spm(m_b+l_max,b,z);

                  //sp(9) second term
                  if(gb == gd)
                     if(S_ab == S_de)
                        (*this)(B,i,j) -= sign_ab * norm_ab * norm_de * spm(m_c+l_max,c,z);

               }

               if(gc == gd){

                  //tp(7)
                  if(ge == gz)
                     (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * hard * tpm(S_ab,m_a+m_b,m_a,a,m_b,b,m_e,e,m_z,z);
                  else
                     (*this)(B,i,j) += A * norm_de * sign_ab * hard * tpm(S_ab,m_a+m_b,m_a,a,m_b,b,m_e,e,m_z,z);

               }

               if(gb == gd){

                  //tp(8)
                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,m_a+m_c,m_a,a,m_c,c,m_e,e,m_z,z);

                  if(ga == gc)
                     hulp *= std::sqrt(2.0);

                  if(ge == gz)
                     hulp *= std::sqrt(2.0);

                  (*this)(B,i,j) += A * sign_ab * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

               }

               if(ga == gd){

                  //tp(8)
                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,m_b+m_c,m_b,b,m_c,c,m_e,e,m_z,z);

                  if(gb == gc)
                     hulp *= std::sqrt(2.0);

                  if(ge == gz)
                     hulp *= std::sqrt(2.0);

                  (*this)(B,i,j) += A * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

               }

            }
         }

      }
      else{//then the S = 3/2 block, this should be easy, totally antisymmetrical 

         for(int i = 0;i < gdim(B);++i){

            ga = d2s[B][i][1];
            gb = d2s[B][i][2];
            gc = d2s[B][i][3];

            m_a = SPM::gg2ms(ga,0);
            a = SPM::gg2ms(ga,1);

            m_b = SPM::gg2ms(gb,0);
            b = SPM::gg2ms(gb,1);

            m_c = SPM::gg2ms(gc,0);
            c = SPM::gg2ms(gc,1);

            for(int j = i;j < gdim(B);++j){

               gd = d2s[B][j][1];
               ge = d2s[B][j][2];
               gz = d2s[B][j][3];

               m_d = SPM::gg2ms(gd,0);
               d = SPM::gg2ms(gd,1);

               m_e = SPM::gg2ms(ge,0);
               e = SPM::gg2ms(ge,1);

               m_z = SPM::gg2ms(gz,0);
               z = SPM::gg2ms(gz,1);

               (*this)(B,i,j) = 0.0;

               if(i == j)
                  (*this)(B,i,j) += ward;

               if(gc == gz){

                  //tp(1)
                  (*this)(B,i,j) += A * tpm(1,m_a+m_b,m_a,a,m_b,b,m_d,d,m_e,e);

                  //sp(1) first part
                  if(gb == ge)
                     (*this)(B,i,j) -= spm(m_a+l_max,a,d);

                  //sp(4) first part
                  if(gb == gd)
                     (*this)(B,i,j) += spm(m_a+l_max,a,e);

                  //sp(5)
                  if(ga == gd)
                     (*this)(B,i,j) -= spm(m_b+l_max,b,e);

               }

               if(gb == gz){

                  //tp(2)
                  (*this)(B,i,j) -= A * tpm(1,m_a+m_c,m_a,a,m_c,c,m_d,d,m_e,e);

                  //sp(1) second part
                  if(gc == ge)
                     (*this)(B,i,j) += spm(m_a+l_max,a,d);

                  //sp(4) second part
                  if(gc == gd)
                     (*this)(B,i,j) -= spm(m_a+l_max,a,e);

                  //sp(6)
                  if(ga == gd)
                     (*this)(B,i,j) += spm(m_c+l_max,c,e);

               }

               if(gc == ge){

                  //tp(4)
                  (*this)(B,i,j) -= A * tpm(1,m_a+m_b,m_a,a,m_b,b,m_d,d,m_z,z);

                  //sp(7) first part
                  if(gb == gd)
                     (*this)(B,i,j) -= spm(m_a+l_max,a,z);

                  //sp(8) first part
                  if(ga == gd)
                     (*this)(B,i,j) += spm(m_b+l_max,b,z);

               }

               if(gb == ge){

                  //tp(5)
                  (*this)(B,i,j) += A * tpm(1,m_a+m_c,m_a,a,m_c,c,m_d,d,m_z,z);

                  //sp(7) second part
                  if(gc == gd)
                     (*this)(B,i,j) += spm(m_a+l_max,a,z);

                  //sp(9) first part
                  if(ga == gd)
                     (*this)(B,i,j) -= spm(m_c+l_max,c,z);

               }

               //tp(7)
               if(gc == gd)
                  (*this)(B,i,j) += A * tpm(1,m_a+m_b,m_a,a,m_b,b,m_e,e,m_z,z);

               //tp(8)
               if(gb == gd)
                  (*this)(B,i,j) -= A * tpm(1,m_a+m_c,m_a,a,m_c,c,m_e,e,m_z,z);

               //tp(9)
               if(ga == gd)
                  (*this)(B,i,j) += A * tpm(1,m_b+m_c,m_b,b,m_c,c,m_e,e,m_z,z);

            }
         }

      }

   }

   this->symmetrize();

}

/**
 * The T1-map: maps a TPM object (tpm) on a DPM object (*this). 
 * @param tpm input TPM
 */
void DPM::T(const TPM &tpm){

   double a = 1.0;
   double b = 1.0/(N*(N - 1.0));
   double c = 1.0/(N - 1.0);

   this->T(a,b,c,tpm);

}

/** 
 * The hat function maps a TPM object tpm to a DPM object (*this) so that bar(this) = tpm,
 * The inverse of the TPM::bar function. It is a T1-like map.
 * @param tpm input TPM
 */
void DPM::hat(const TPM &tpm){

   double a = 1.0/(M - 4.0);
   double b = 1.0/((M - 4.0)*(M - 3.0)*(M - 2.0));
   double c = 1.0/((M - 4.0)*(M - 3.0));

   this->T(a,b,c,tpm);

}
