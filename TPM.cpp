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

            (*this)(B,i,j) += si.gV(a_me,b_me,c_me,d_me) + sign * si.gV(a_me,b_me,d_me,c_me);

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
 * Construct the right hand side of the Newton equation for the determination of the search direction, 
 * the gradient of the potential:
 * @param t scaling factor of the potential
 * @param ham Hamiltonian of the current problem
 * @param P SUP matrix containing the inverse of the constraint matrices (carrier space matrices).
 */
void TPM::constr_grad(double t,const TPM &ham,const SUP &P,const LinIneq &li){

   //eerst P conditie 
   *this = P.gI();

#ifdef __Q_CON
   //de Q conditie toevoegen
   TPM hulp;

   hulp.Q(1,P.gQ());

   *this += hulp;
#endif

#ifdef __G_CON
   hulp.G(P.gG());

   *this += hulp;
#endif

#ifdef __T1_CON
   hulp.T(P.gT1());

   *this += hulp;
#endif

#ifdef __T2_CON
   hulp.T(P.gT2());

   *this +=hulp;
#endif

   //and the linear constraints
   for(int i = 0;i < li.gnr();++i)
      this->daxpy(1.0/li.constraint(i),li[i].gI());

   this->dscal(t);

   *this -= ham;

   this->proj_Tr();

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
 * solve the Newton equations for the determination of the search direction,
 * @param t scaling factor of the potential
 * @param P SUP matrix containing the inverse of the constraint matrices (carrier space matrices).
 * @param b right hand side (the gradient constructed int TPM::constr_grad)
 * @return nr of iterations needed to converge to the desired accuracy
 */
int TPM::solve(double t,const SUP &P,TPM &b,const LinIneq &li){

   int iter = 0;

   //delta = 0
   *this = 0;

   //residu:
   TPM r(b);

   //norm van het residu
   double rr = r.ddot(r);

   //enkele variabelen
   double rr_old,ward;

   TPM Hb;

   while(rr > 1.0e-10){ 

      Hb.H(t,b,P,li);

      ward = rr/b.ddot(Hb);

      //delta += ward*b
      this->daxpy(ward,b);

      //r -= ward*Hb
      r.daxpy(-ward,Hb);

      //nieuwe variabelen berekenen en oude overdragen
      rr_old = rr;
      rr = r.ddot(r);

      //nieuwe b nog:
      b.dscal(rr/rr_old);

      b += r;

      ++iter;

   }

   return iter;

}

/**
 * The hessian-map of the Newton system:
 * @param t potential scaling factor
 * @param b the TPM on which the hamiltonian will work, the image will be put in (*this)
 * @param P the SUP matrix containing the constraints, (can be seen as the metric).
 */
void TPM::H(double t,const TPM &b,const SUP &P,const LinIneq &li){

   //eerst de P conditie:
   this->L_map(P.gI(),b);

#ifdef __Q_CON
   TPM hulp;

   //maak Q(b)
   TPM Q_b;
   Q_b.Q(1,b);

   //stop Q(rdm)^{-1}Q(b)Q(rdm)^{-1} in hulp
   hulp.L_map(P.gQ(),Q_b);

   //maak Q(hulp) en stop in Q_b
   Q_b.Q(1,hulp);

   //en tel op bij this
   *this += Q_b;
#endif

#ifdef __G_CON
   //hulpje voor het PHM stuk
   PHM hulp_ph;
   PHM G_b;

   //stop G(b) in G_b
   G_b.G(b);

   //bereken G(rdm)^{-1}G(b)G(rdm)^{-1} en stop in hulp_ph
   hulp_ph.L_map(P.gG(),G_b);

   //tenslotte nog de antisymmetrische G hierop:
   hulp.G(hulp_ph);

   //en optellen bij this
   *this += hulp;
#endif
   
#ifdef __T1_CON
   //hulpjes voor het DPM stuk
   DPM hulp_dp;
   DPM T1_b;

   //stop T1(b) in T1_b
   T1_b.T(b);

   hulp_dp.L_map(P.gT1(),T1_b);

   hulp.T(hulp_dp);

   *this += hulp;
#endif

#ifdef __T2_CON
   PPHM hulp_pph;
   PPHM T2_b;

   T2_b.T(b);

   hulp_pph.L_map(P.gT2(),T2_b);

   hulp.T(hulp_pph);

   *this+=hulp;
#endif

   //linear constraints
   for(int i = 0;i < li.gnr();++i)
      this->daxpy( b.ddot(li[i].gI()) / (li.constraint(i)*li.constraint(i)) , li[i].gI() );

   //nog schalen met t:
   this->dscal(t);

   //en projecteren op spoorloze ruimte
   this->proj_Tr();

}

/**
 * perform a line search what step size in along the Newton direction is ideal.
 * @param t potential scaling factor
 * @param P SUP matrix containing the inverse of the constraints (carrier space matrices)
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,SUP &P,const TPM &ham,const LinIneq &li){

   double tolerance = 1.0e-5*t;

   if(tolerance < 1.0e-12)
      tolerance = 1.0e-12;

   //neem de wortel uit P
   P.sqrt(1);

   //maak eerst een SUP van delta
   SUP S_delta;

   S_delta.fill(*this);

   //hulpje om dingskes in te steken:
   SUP hulp;

   hulp.L_map(P,S_delta);

   EIG eigen(hulp);

   double a = 0;

   double b = -1.0/eigen.min();

   double c(0);

   double ham_delta = ham.ddot(*this);

   LinIneq li_epsi;
   li_epsi.fill(*this);

   double lin_end = li.min_end(li_epsi);

   if(lin_end < b)
      b = lin_end;

   while(b - a > tolerance){

      c = (b + a)/2.0;

      if( (ham_delta - t*eigen.lsfunc(c) - t*li.lsfunc(c,li_epsi)) < 0.0)
         a = c;
      else
         b = c;

   }

   return c;

}

/**
 * perform a line search what step size in along the Newton direction is ideal, this one is used for extrapolation.
 * @param t potential scaling factor
 * @param rdm TPM containing the current approximation of the rdm.
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,const TPM &rdm,const TPM &ham,const LinIneq &li){

   SUP P;

   P.fill(rdm);

   P.invert();

   return this->line_search(t,P,ham,li);

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
 * fill the TPM object with the S^2 matrix
 */
void TPM::set_S_2(){

   *this = 0.0;

   for(int B = 0;B < gnr();++B){

      if(B2SM[B][0] == 0){

         for(int i = 0;i < gdim(B);++i)
            (*this)(B,i,i) = -1.5 * (N - 2.0)/(N - 1.0);

      }
      else{

         for(int i = 0;i < gdim(B);++i)
            (*this)(B,i,i) = -1.5 * (N - 2.0)/(N - 1.0) + 2.0;

      }

   }

}
