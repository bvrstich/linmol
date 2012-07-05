#ifndef CI_TPM_H
#define CI_TPM_H

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
 * This class CI_TPM is a class written for the two-particle electron repulsion integrals in a Cartesian basis set.
 * For use in the CartInt class
 */
class CI_TPM : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output the output stream e.g. cout
    * @param ci_p the CI_TPM you want to print
    */
   friend ostream &operator<<(ostream &output,CI_TPM &ci_p);

   public:
      
      //constructor
      CI_TPM();

      //copy constructor
      CI_TPM(const CI_TPM &);

      //destructor
      virtual ~CI_TPM();

      using Matrix::operator=;

      using Matrix::operator();

      double operator()(int,int,int,int) const;

      double operator()(int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int) const;

      void calc_eri();

      static int gtp_dim();

      static int gt2s(int,int);

      static int gs2t(int,int);

      static int gn_max();

      static int gl_max();

      static void init();

      static void clear();

   private:

      //!dimension of the tp matrices
      static int tp_dim;

      //!list relating tp to sp indices
      static vector< vector<int> > t2s;

      //!list relating tp to sp indices
      static int **s2t;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
