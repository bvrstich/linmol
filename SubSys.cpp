#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

int SubSys::M;
int SubSys::N;

/**
 * initialize the static lists and variables
 * @param M_in input dimension of sp space
 * @param N_in input nr of particles
 */
void SubSys::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * Constructor of a SubSys object, takes a subsystem on a certrain core
 * @param si_in input SphInt object
 * @param core number of the core
 */
SubSys::SubSys(int core,const SphInt &si_in){

   this->core = core;

   subham = new TPM();

   si = new SphInt(si_in);

   //set the "non-subsystem" matrixelements to zero:
   for(int s_i = 0;s_i < si->gdim();++s_i)
      for(int s_j = s_i;s_j < si->gdim();++s_j){

         if(SphInt::gs2inlm(s_i,0) != core || SphInt::gs2inlm(s_j,0) != core){

            si->gT()(s_i,s_j) = 0.0;
            si->gU()(s_i,s_j) = 0.0;

         }

      }

   int s_i,s_j,s_k,s_l;

   for(int t_i = 0;t_i < si->gdim() * si->gdim();++t_i){

      s_i = si->gt2s(t_i,0);
      s_j = si->gt2s(t_i,1);

      for(int t_j = t_i;t_j < si->gdim() * si->gdim();++t_j){

         s_k = si->gt2s(t_j,0);
         s_l = si->gt2s(t_j,1);

         if(SphInt::gs2inlm(s_i,0) != core || SphInt::gs2inlm(s_j,0) != core || SphInt::gs2inlm(s_k,0) != core || SphInt::gs2inlm(s_l,0) != core)
            si->gV()(t_i,t_j) = 0.0;

      }
   }

   si->gT().symmetrize();
   si->gU().symmetrize();
   si->gV().symmetrize();

   si->orthogonalize();

   subham->molecule(*si);

   //make the subsystem occupation operator
   subocc = new TPM();

   si->gS().invert();

   subocc->subocc_op(core,si->gS());

}


/**
 * copy constructor
 * @param ss_copy The SubSys object to be copied
 */
SubSys::SubSys(const SubSys &ss_copy){

   this->core = ss_copy.gcore();

   subham = new TPM(ss_copy.gsubham());

   E = ss_copy.gE();

}

/**
 * destructor
 */
SubSys::~SubSys(){

   delete subham;

   delete si;

}

/**
 * @param option if == 0 return N, if == 1 return E(N)
 * @return the energy corresponding to the index "particle number"
 */
double SubSys::gE(int index,int option) const{

   return E[index][option];

}

/**
 * @return the full vector< vector<int> > E object containing the subsystem energies and particle numbers
 */
const vector< vector<int> > &SubSys::gE() const {

   return E;

}

/**
 * @return the full vector< vector<int> > E object containing the subsystem energies and particle numbers
 */
vector< vector<int> > &SubSys::gE() {

   return E;

}

/**
 * @return the index of the core the subsystem belongs to
 */
int SubSys::gcore() const {

   return core;

}

/**
 * @return the subsystem Hamiltonian
 */
const TPM &SubSys::gsubham() const {

   return *subham;

}

/**
 * @return the subsystem Hamiltonian
 */
TPM &SubSys::gsubham() {

   return *subham;

}

/**
 * @return the subsystem Hamiltonian
 */
const TPM &SubSys::gsubocc() const {

   return *subocc;

}

/**
 * @return the subsystem Hamiltonian
 */
TPM &SubSys::gsubocc() {

   return *subocc;

}

/**
 * get the occupation of the subsystem: watch out, need overlapmatrix for this, and preferably the inverse sqrt !
 */
double SubSys::subocc_func(const TPM &tpm) const {

   SPM spm;
   spm.bar(1.0/(N - 1.0),tpm);

   double ward = 0.0;

   int i,m;
   int i_a,n_a,l_a,m_a;
   int i_b,n_b,l_b,m_b;

   for(int s_i = 0;s_i < SphInt::gdim();++s_i){

      i = SphInt::gs2inlm(s_i,0);
      m = SphInt::gs2inlm(s_i,3);

      if(i == core){

         for(int s_a = 0;s_a < SphInt::gdim();++s_a){

            i_a = SphInt::gs2inlm(s_a,0);
            n_a = SphInt::gs2inlm(s_a,1);
            l_a = SphInt::gs2inlm(s_a,2);
            m_a = SphInt::gs2inlm(s_a,3);

            for(int s_b = 0;s_b < SphInt::gdim();++s_b){

               i_b = SphInt::gs2inlm(s_b,0);
               n_b = SphInt::gs2inlm(s_b,1);
               l_b = SphInt::gs2inlm(s_b,2);
               m_b = SphInt::gs2inlm(s_b,3);

               if(m_a == m && m_b == m)
                  ward += si->gS()(s_i,s_a) * spm[m + SphInt::gl_max()](SPM::ginl2s(m_a,i_a,n_a,l_a),SPM::ginl2s(m_b,i_b,n_b,l_b)) * si->gS()(s_b,s_i);

            }
         }

      }

   }

   return ward;

}

/**
 * set the subsystem energies for Be with occupations 3, 4 and 5.
 */
void SubSys::setBe(){

   vector<int> v(2);

   //3
   v[0] = 3;
   v[1] = -14.27612024;

   E.push_back(v);

   //4
   v[0] = 4;
   v[1] = -14.61749376;

   E.push_back(v);

   //5
   v[0] = 5;
   v[1] = -14.58336127;

   E.push_back(v);

}

/**
 * set the subsystem energies for B with occupations 3, 4 and 5.
 */
void SubSys::setB(){

   vector<int> v(2);

   //3
   v[0] = 3;
   v[1] = -23.372643;

   E.push_back(v);

   //4
   v[0] = 4;
   v[1] = -24.29393637;

   E.push_back(v);

   //5
   v[0] = 5;
   v[1] = -24.60479038;

   E.push_back(v);

}
