#ifndef SPHINT_H
#define SPHINT_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"
#include "input.h"

/**
 * @author Brecht Verstichel
 * @date 05-07-2012\n\n
 * This class SphInt is a class written for the storage of matrixelements in a spherical Gaussian basis
 * First the matrixelements are constructed in Cartesian orbitals in CartInt, then they are transformed to a spherical basis in this class
 */
class SphInt {

   public:
      
      static const SI_SPM &gS();

      static const SI_SPM &gT();

      static const SI_SPM &gU(int);

      static const SI_TPM &gV();

      static void orthogonalize();

      static void print();

      static void fill();

      static void init();

      static void clear();

   private:

      //!overlapmatrix
      static SI_SPM *S;

      //!kinetic energy matrix
      static SI_SPM *T;

      //!nuclear attraction energy matrix
      static SI_SPM **U;

      //!electronic repulsion matrix
      static SI_TPM *V;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
