#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_sf_coupling.h>

#include "include.h"

using std::cout;
using std::endl;

double *Tools::x6j;
double *Tools::x9j;

int Tools::M;
int Tools::N;


/**
 * initialize and allocate the static variables
 * @param M_in nr of sites
 * @param N_in nr of particles
 */
void Tools::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

   x6j = new double [16];

   x9j = new double [16];

   //fill the x6j list
   for(int S_ = 0;S_ < 2;++S_)
      for(int S = 0;S < 2;++S)
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd)
               x6j[8*S_ + 4*S + 2*S_ab + S_cd] = gsl_sf_coupling_6j(2*S_ + 1,1,2*S_ab,2*S + 1,1,2*S_cd);

   //fill the x9j list
   for(int S = 0;S < 2;++S)
      for(int Z = 0;Z < 2;++Z)
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd)
               x9j[8*S + 4*Z + 2*S_ab + S_cd] = gsl_sf_coupling_9j(1,1,2*Z,2*S + 1,2*S_ab,1,2*S_cd,1,1);

}

/**
 * deallocate the static lists and objects
 */
void Tools::clear(){

   delete [] x6j;

   delete [] x9j;

}

/**
 * store only the six j symbols needed
 * @param S_ 0 or 1 means 1/2 or 3/2
 * @param S 0 or 1 means 1/2 or 3/2
 * @param S_ab 0 or 1 means 0 or 1
 * @param S_cd 0 or 1 means 0 or 1
 * @return the 6j symbol
 */
double Tools::g6j(int S_,int S,int S_ab,int S_cd){

   return x6j[8*S_ + 4*S + 2*S_ab + S_cd];

}

/**
 * store only the nine j symbols needed
 * @param S 0 or 1 means 1/2 or 3/2
 * @param Z 0 or 1 means 0 or 1
 * @param S_ab 0 or 1 means 0 or 1
 * @param S_cd 0 or 1 means 0 or 1
 * @return the 6j symbol
 */
double Tools::g9j(int S,int Z,int S_ab,int S_cd){

   return x9j[8*S + 4*Z + 2*S_ab + S_cd];

}

/* vim: set ts=3 sw=3 expandtab :*/
