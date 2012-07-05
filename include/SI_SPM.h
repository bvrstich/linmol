#ifndef SI_SPM_H
#define SI_SPM_H

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
 * This class SI_SPM is a class written for the one-particle integrals in a Cartesian basis set, i.e. S,T and U.
 * For use in the CartInt class
 */
class SI_SPM : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output the output stream e.g. cout
    * @param si_p the SI_SPM you want to print
    */
   friend ostream &operator<<(ostream &output,SI_SPM &si_p);

   public:
      
      //constructor
      SI_SPM();

      //copy constructor
      SI_SPM(const SI_SPM &);

      //destructor
      virtual ~SI_SPM();

      using Matrix::operator=;

      using Matrix::operator();

      double operator()(int,int,int,int,int,int,int,int) const;

      void transform(const CI_SPM &);

      static int gs2inlm(int,int);

      static int ginlm2s(int,int,int,int);

      static int gdim();

      static int gn_max();

      static int gl_max();

      static void init();

      static void clear();

   private:

      //!static objects needed to construct and destruct all the lists
      static int l_max,n_max;

      //!list keeping the relationship between a single-particle index s and physical quantumnumbers inlm.
      static vector< vector<int> > s2inlm;

      //!reverse list keeping the relationship between a single-particle index s and physical quantumnumbers inlm.
      static int ****inlm2s;

      //!dimension of the basisset
      static int dim;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
