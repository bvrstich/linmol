#ifndef CARTINT_H
#define CARTINT_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"
#include "input.h"

class CI_SPM;
class CI_TPM;

/**
 * @author Brecht Verstichel
 * @date 13-04-2012\n\n
 * This class CartInt is a class written for the storage of matrixelements in a cartesian Gaussian basis
 * In this class the different integrals are calculated using the libint library
 */
class CartInt {

   public:
      
      static const CI_SPM &gS();

      static const CI_SPM &gT();

      static const CI_SPM &gU(int);

      static const CI_TPM &gV();

      static void norm();

      static void orthogonalize();

      static void calc_integrals();

      static void print();

      static void init();

      static void clear();

   private:

      //!overlapmatrix
      static CI_SPM *S;

      //!kinetic energy matrix
      static CI_SPM *T;

      //!nuclear attraction energy matrix
      static CI_SPM **U;

      //!electronic repulsion matrix
      static CI_TPM *V;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
