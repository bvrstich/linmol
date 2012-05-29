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
 * @param si_in input SphInt object: Make sure that it is NOT orthogonalized yet!
 * @param core number of the core
 */
SubSys::SubSys(int core,const SphInt &si_in){

   this->core = core;

   //make the s2f list and get the subsystem dimension
   n = 0;

   for(int a = 0;a < si_in.gdim();++a){

      if(SphInt::gs2inlm(a,0) == core){

         s2f.push_back(a);
         ++n;

      }

   }

   si = new SphInt(si_in);

   //overlap matrix
   S = new Matrix(n);

   //copy the right elements 
   for(int i = 0;i < n;++i)
      for(int j = i;j < n;++j)
         (*S)(i,j) = si_in.gS()(s2f[i],s2f[j]);

   S->symmetrize();

   //set the correct elements to zero in the T,U and V matrices

   //first T and U
   for(int a = 0;a < si->gdim();++a)
      for(int b = a;b < si->gdim();++b){

         if(SphInt::gs2inlm(a,0) != core || SphInt::gs2inlm(b,0) != core){

            si->gT()(a,b) = 0.0;

            for(int i = 0;i < si->gN_Z();++i)
               si->gU(i)(a,b) = 0.0;

         }

      }

   si->gT().symmetrize();

   for(int i = 0;i < si->gdim();++i)
      si->gU(i).symmetrize();

   //then V
   int a,b,c,d;

   for(int i = 0;i < si->gdim()*si->gdim();++i){

      a = SphInt::gt2s(i,0);
      b = SphInt::gt2s(i,1);

      for(int j = i;j < si->gdim()*si->gdim();++j){

         c = SphInt::gt2s(j,0);
         d = SphInt::gt2s(j,1);

         
         if(SphInt::gs2inlm(a,0) != core || SphInt::gs2inlm(b,0) != core || SphInt::gs2inlm(c,0) != core || SphInt::gs2inlm(d,0) != core)
            si->gV()(i,j) = 0.0;

      }
   }

   si->gV().symmetrize();

   //construct the L matrix: linear transformation between the orthogonal and non-orthogonal basis
   si->gS().sqrt(1);

   //rectangular matrix! More rows than columns
   W = new double [n*M/2];

   //now construct it
   S->invert();

   for(int s = 0;s < n;++s)
      for(int f = 0;f < M/2;++f){

         W[s + n*f] = 0.0;

         for(int i = 0;i < n;++i)
            W[s + n*f] += si->gS()(f,s2f[i]) * (*S)(i,s);

      }

}

/**
 * copy constructor:
 * @param ss_copy The SubSys object to be copied
 */
SubSys::SubSys(const SubSys &ss_copy){

   this->core = ss_copy.gcore();

   n = ss_copy.gn();
   s2f = ss_copy.gs2f();

   //matrix already inverted!
   S = new Matrix(ss_copy.gS());

   //si already "cut down on subsystem space" and the overlap matrix has been sqrted
   si = new SphInt(ss_copy.gsi());

   //rectangular matrix! More rows than columns
   W = new double [n*M/2];

   for(int s = 0;s < n;++s)
      for(int f = 0;f < M/2;++f)
         W[s + n*f] = ss_copy.gW(f,s);

}

/**
 * destructor
 */
SubSys::~SubSys(){

   delete si;

   delete S;

   delete [] W;

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
 * get the occupation of the subsystem: watch out, need overlapmatrix for this, and preferably the inverse sqrt !
 */
double SubSys::subocc_func(const TPM &tpm) const {

   SPM spm;
   spm.bar(1.0/(N - 1.0),tpm);

   //rough test
   Matrix N_sub(M/2);

   for(int ga = 0;ga < M/2;++ga)
      for(int gb = 0;gb < M/2;++gb){

         N_sub(ga,gb) = 0.0;

         for(int a = 0;a < n;++a)
            for(int b = 0;b < n;++b)
               N_sub(ga,gb) += si->gS()(ga,s2f[a]) * (*S)(a,b) * si->gS()(gb,s2f[b]);

      }

   //convert to SPM
   SPM n_sub;
   n_sub.si21dm(N_sub);

   return spm.ddot(n_sub);

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

int SubSys::gn() const {

   return n;

}

const vector<int> &SubSys::gs2f() const {

   return s2f;

}

/** 
 * @return the overlapmatrix, const version
 */
const Matrix &SubSys::gS() const { 

   return *S;

}

/** 
 * @return the overlapmatrix
 */
Matrix &SubSys::gS() { 

   return *S;

}

/**
 * access to the W thingy
 * @param f full system index
 * @param s subsystem index
 */
double SubSys::gW(int f ,int s) const {

   return W[s + n*f];

}

/**
 * change the W thingy
 * @param f full system index
 * @param s subsystem index
 * @param value fill in the value you want at index W(f,s)
 */
void SubSys::sW(int f ,int s,double value) {

   W[s + n*f] = value;

}

/**
 * access to the s2f list from outside
 * @param s input subsystem index 
 * @return f the full system index corresponding to s
 */
int SubSys::gs2f(int s) const{

   return s2f[s];

}

/**
 * @return the subsystem SphInt object
 */
const SphInt &SubSys::gsi() const {

   return *si;

}

/**
 * @return the subsystem SphInt object
 */
SphInt &SubSys::gsi() {

   return *si;

}

/**
 * transform the T,U and V components of the SphInt object so that they are a good projector in an orthogonal basis on the subsystem space.
 */
void SubSys::orthogonalize(){

   //first T
   Matrix backup_T(si->gdim());

   for(int a = 0;a < si->gdim();++a)
      for(int b = 0;b < si->gdim();++b){

         backup_T(a,b) = 0.0;

         //loop over subsystem index
         for(int sa = 0;sa < n;++sa)
            backup_T(a,b) += gW(a,sa) * si->gT()(s2f[sa],b);

      }

   for(int a = 0;a < si->gdim();++a)
      for(int b = a;b < si->gdim();++b){

         si->gT()(a,b) = 0.0;

         //loop over subsystem index
         for(int sb = 0;sb < n;++sb)
            si->gT()(a,b) += backup_T(a,s2f[sb]) * gW(b,sb);

      }

   si->gT().symmetrize();

}
