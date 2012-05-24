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

   int n = 0;

   //find subsystem dimension
   for(int i = 0;i < M;++i){

   }

}


/**
 * copy constructor
 * @param ss_copy The SubSys object to be copied
 */
SubSys::SubSys(const SubSys &ss_copy){

   this->core = ss_copy.gcore();

}

/**
 * destructor
 */
SubSys::~SubSys(){

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

   return 0.0;

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
