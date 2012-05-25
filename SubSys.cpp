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
 * @param si input SphInt object: Make sure that it is NOT orthogonalized yet!
 * @param core number of the core
 */
SubSys::SubSys(int core,const SphInt &si){

   this->core = core;

   n = 0;

   for(int s_i = 0;s_i < si.gdim();++s_i){

      if(SphInt::gs2inlm(s_i,0) == core){

         s2f.push_back(s_i);
         ++n;

      }

   }

   S = new Matrix(n);
   T = new Matrix(n);
   U = new Matrix(n);

   V = new Matrix(n*n);


   //copy the right elements 
   for(int i = 0;i < n;++i)
      for(int j = i;j < n;++j){

         (*S)(i,j) = si.gS()(s2f[i],s2f[j]);
         (*T)(i,j) = si.gS()(s2f[i],s2f[j]);
         (*U)(i,j) = si.gS()(s2f[i],s2f[j]);

      }

   S->symmetrize();
   T->symmetrize();
   U->symmetrize();

   //construct "two-particle" lists on subspace
   vector<int> v(2);

   s2t = new int * [n];

   for(int i = 0;i < n;++i)
      s2t[i] = new int [n];

   int iter = 0;

   for(int i = 0;i < n;++i)
      for(int j = 0;j < n;++j){

         v[0] = i;
         v[1] = j;

         t2s.push_back(v);

         s2t[i][j] = iter;

         ++iter;

      }

   int s_i,s_j,s_k,s_l;

   for(int t_i = 0;t_i < n*n;++t_i){

      s_i = t2s[t_i][0];
      s_j = t2s[t_i][1];

      for(int t_j = t_i;t_j < n*n;++t_j){

         s_k = t2s[t_j][0];
         s_l = t2s[t_j][1];

         (*V)(t_i,t_j) = si.gV(s2f[s_i],s2f[s_j],s2f[s_k],s2f[s_l]);

      }
   }

   V->symmetrize();

   //construct the L matrix
   L = new Matrix(si.gS());
   
   L->sqrt(-1);

   //rectangular matrix! More rows than columns
   W = new double [n*M/2];

   //now construct it
   S->invert();

   for(int s = 0;s < n;++s)
      for(int f = 0;f < M/2;++f){

         W[s + n*f] = 0.0;

         for(int i = 0;i < n;++i)
            W[s + n*f] += (*L)(f,i) * (*S)(i,s);

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

   S = new Matrix(ss_copy.gS());
   T = new Matrix(ss_copy.gT());
   U = new Matrix(ss_copy.gU());
   V = new Matrix(ss_copy.gV());

   //construct "two-particle" lists on subspace
   vector<int> v(2);

   s2t = new int * [n];

   for(int i = 0;i < n;++i)
      s2t[i] = new int [n];

   int iter = 0;

   for(int i = 0;i < n;++i)
      for(int j = 0;j < n;++j){

         v[0] = i;
         v[1] = j;

         t2s.push_back(v);

         s2t[i][j] = iter;

         ++iter;

      }

   //construct the L matrix
   L = new Matrix(ss_copy.gL());

   //rectangular matrix! More rows than columns
   W = new double [n*M/2];

   for(int s = 0;s < n;++s)
      for(int f = 0;f < M/2;++f){

         W[s + n*f] = 0.0;

         for(int i = 0;i < n;++i)
            W[s + n*f] += (*L)(f,i) * (*S)(i,s);

      }

}

/**
 * destructor
 */
SubSys::~SubSys(){

   delete S;
   delete T;
   delete U;
   delete V;

   delete L;

   for(int i = 0;i < n;++i)
      delete [] s2t[i];

   delete [] s2t;

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

   //delete this line later: need positive S here
   S->invert();

   //rough test
   Matrix N_sub(M/2);

   for(int ga = 0;ga < M/2;++ga)
      for(int gb = 0;gb < M/2;++gb){

         N_sub(ga,gb) = 0.0;

         for(int a = 0;a < n;++a)
            for(int b = 0;b < n;++b)
               N_sub(ga,gb) += gW(ga,a) * gW(gb,b) * (*S)(a,b);

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
 * @return the kinetic energy matrix, const version
 */
const Matrix &SubSys::gT() const { 

   return *T; 
}

/** 
 * @return the kinetic energy matrix
 */
Matrix &SubSys::gT() { 

   return *T;

}

/** 
 * @return the nuclear attraction matrix, const version
 */
const Matrix &SubSys::gU() const { 

   return *U; 
}

/** 
 * @return the nuclear attraction matrix
 */
Matrix &SubSys::gU() { 

   return *U;

}

/** 
 * @return the electronic repulsion matrix
 */
const Matrix &SubSys::gV() const { 

   return *V; 
}

/** 
 * @return the electronic repulsion matrix
 */
Matrix &SubSys::gV() { 

   return *V;

}

/** 
 * @return the linear transformation matrix
 */
const Matrix &SubSys::gL() const { 

   return *V; 

}

/** 
 * @return the linear transformation matrix
 */
Matrix &SubSys::gL() { 

   return *V;

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
