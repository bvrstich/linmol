#ifndef CI_TPPM_H
#define CI_TPPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"
#include "input.h"

/**
 * @author Brecht Verstichel
 * @date 03-07-2012\n\n
 * This class is written for the evaluation of two-electron repulsion integrals!
 * This matrix stores the integrals between all the primitives of the orbitals individually.
 * Warning, chemical notation is used here because it allows for a more straightforward compact storage.
 */
class CI_TPPM : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output the output stream e.g. cout
    * @param ci_p the CI_TPPM you want to print
    */
   friend ostream &operator<<(ostream &output,CI_TPPM &ci_p);

   public:
      
      //constructor
      CI_TPPM();

      //copy constructor
      CI_TPPM(const CI_TPPM &);

      //destructor
      virtual ~CI_TPPM();

      using Matrix::operator=;

      using Matrix::operator();

      void constr_eri();

      double operator()(int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int) const;

      double &operator()(int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int);

      static void init();

      static void clear();

   private:
      
      static vector< vector<int> > tomi2omi;

      static int **omi2tomi;

      static vector< vector<int> > g2ijk;

      static int ***ijk2g;

      static vector< vector<int> > t2s;

      static int **s2t;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
