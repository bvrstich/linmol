#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

double *LibInt::F;

const int LibInt::MAXFAC = 100;
const double LibInt::EPS = 1.0e-17;
const int LibInt::L_MAX = 5;

/**
 * initialize the lists
 */
void LibInt::init(){

   //allocate the F array needed for the calc_f function
   F = new double [21];

   
}

/**
 * deallocate every static list
 */
void LibInt::clear(){

   delete [] F;

}

/**
 * @return the static variable defining the size of the double faculty array
 */
int LibInt::gMAXFAC(){

   return MAXFAC;

}

/**
 * @return the static variable defining the precision with which to calculate the boys function F_n(t)
 */
int LibInt::gEPS(){

   return EPS;

}

/**
 * @return the static constant variable defining the maximal angular momentum the program can handle
 */
int LibInt::gL_MAX(){

   return L_MAX;

}

/**
 * calculate the boys function F_n(t)
 */
void LibInt::calc_f(int n, double t){

   int i, m;
   int m2;
   double t2;
   double num;
   double sum;
   double term1;
   static double K = 1.0 / M_2_SQRTPI;
   double et;

   if (t > 20.0) { /* For big t's do upward recursion */
      t2 = 2 * t;
      et = exp(-t);
      t = sqrt(t);
      F[0] = K * erf(t) / t;
      for (m = 0; m <= n - 1; m++) {
         F[m + 1] = ((2 * m + 1) * F[m] - et) / (t2);
      }
   } else { /* For smaller t's compute F with highest n using
               asymptotic series (see I. Shavitt in
               Methods in Computational Physics, ed. B. Alder eta l,
               vol 2, 1963, page 8) */
      et = exp(-t);
      t2 = 2 * t;
      m2 = 2 * n;
      num = Tools::gdfac(m2);
      i = 0;
      sum = 1.0 / (m2 + 1);
      do {
         i++;
         num = num * t2;
         term1 = num / Tools::gdfac(m2 + 2 * i + 2);
         sum += term1;
      } while (fabs(term1) > EPS && i < MAXFAC);
      F[n] = sum * et;
      for (m = n - 1; m >= 0; m--) { /* And then do downward recursion */
         F[m] = (t2 * F[m + 1] + et) / (2 * m + 1);
      }
   }

}

/**
 * @return the  Boys function F_n(t), be carefull, first run calc_f, otherwise this array will contain rubbish!
 */
double LibInt::gF(int n){

   return F[n];

}
