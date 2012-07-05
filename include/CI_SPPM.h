#ifndef CI_SPPM_H
#define CI_SPPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"
#include "input.h"

/**
 * @author Brecht Verstichel
 * @date 19-04-2012\n\n
 * This class is written for single-particle integrals of contracted gaussians.
 * This matrix stores the integrals between all the primitives of the orbitals individually
 */
class CI_SPPM : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output the output stream e.g. cout
    * @param ci_p the CI_SPPM you want to print
    */
   friend ostream &operator<<(ostream &output,CI_SPPM &ci_p);

   public:
      
      //constructor
      CI_SPPM();

      //copy constructor
      CI_SPPM(const CI_SPPM &);

      //destructor
      virtual ~CI_SPPM();

      using Matrix::operator=;

      using Matrix::operator();

      double operator()(int,int,int,int,int,int,int,int,int,int,int,int) const;

      void constr_overlap();

      void constr_kinetic();

      static int gomi2ijkxyz(int,int);

      static int gijkxyz2omi(int,int,int,int,int,int);

      static int gamlist(int);

      static int gammax(int,int);

      static int gdim();

      static void init();

      static void clear();

   private:

      //!static list relating the larger matrix index to the single-particle index and the primitive index
      static vector< vector<int> > omi2ijkxyz;

      //!reverse list
      static int ******ijkxyz2omi;

      //!list containing the maximal angular momentum corresponding to ij
      static int **ammax;

      //!list containing the angular momentum corresponding to a ijkxyz or omi index.
      static vector<int> amlist;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
