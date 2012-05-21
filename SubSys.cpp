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

   E = new double [3];//I guess just three numbers is standard

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

}


/**
 * copy constructor
 * @param ss_copy The SubSys object to be copied
 */
SubSys::SubSys(const SubSys &ss_copy){

   this->core = ss_copy.gcore();

   subham = new TPM(ss_copy.gsubham());

   E = new double [3];

   for(int i = 0;i < 3;++i)
      E[i] = ss_copy.gE(i);

   delete si;

}

/**
 * destructor
 */
SubSys::~SubSys(){

   delete subham;

   delete [] E;

}

/**
 * @return the energy corresponding to the index "particle number"
 */
double SubSys::gE(int index) const{

   return E[index];

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
 * get the occupation of the subsystem: watch out, need overlapmatrix for this!
 */
double SubSys::subocc(const TPM &tpm) const {

   SPM spm;
   spm.bar(1.0/(N - 1.0),tpm);

   double ward = 0.0;

   int i,n,l,m;

   for(int s_i = 0;s_i < SphInt::gdim();++s_i){

      i = SphInt::gs2inlm(s_i,0);
      n = SphInt::gs2inlm(s_i,1);
      l = SphInt::gs2inlm(s_i,2);
      m = SphInt::gs2inlm(s_i,3);

      if(i == core)
         ward += spm[m + SphInt::gl_max()](SPM::ginl2s(m,i,n,l),);

   }

   return ward;

}
