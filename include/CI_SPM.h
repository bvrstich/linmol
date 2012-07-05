#ifndef CI_SPM_H
#define CI_SPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"
#include "input.h"

/**
 * @author Brecht Verstichel
 * @date 13-04-2012\n\n
 * This class CI_SPM is a class written for the one-particle integrals in a Cartesian basis set, i.e. S,T and U.
 * For use in the CartInt class
 */
class CI_SPM : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output the output stream e.g. cout
    * @param ci_p the CI_SPM you want to print
    */
   friend ostream &operator<<(ostream &output,CI_SPM &ci_p);

   public:
      
      //constructor
      CI_SPM();

      //copy constructor
      CI_SPM(const CI_SPM &);

      //destructor
      virtual ~CI_SPM();

      using Matrix::operator=;

      using Matrix::operator();

      double operator()(int,int,int,int,int,int,int,int,int,int,int,int) const;

      void calc_overlap();

      void calc_kinetic();

      void calc_nucl_attrac(int);

      static int gs2gauss(int,int);

      static int gs2inlxyz(int,int);

      static int ginlxyz2s(int,int,int,int,int,int);

      static int gdim();

      static int gn_max();

      static int gl_max();

      static void init();

      static void clear();

   private:

      //!static objects needed to construct and destruct all the lists
      static int l_max,n_max;

      //!list to switch between matrix index and physical quantum numbers
      static vector< vector<int> > s2inlxyz;

      //!list to switch between matrix index and physical quantum numbers
      static int ******inlxyz2s;

      //!list that converts CI_SPM matrix index to Gaussian object index
      static vector< vector<int> > s2gauss;

      //!dimension of the basisset
      static int dim;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
