#ifndef SPHINT_H
#define SPHINT_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"
#include "input.h"
#include "CartInt.h"

/**
 * @author Brecht Verstichel
 * @date 13-04-2012\n\n
 * This class SphInt is a class written for the storage of matrixelements in a spherical Gaussian basis
 * First the matrixelements are constructed using ThING, then they are transformed to a spherical basis in this class
 */
class SphInt {

   /**
    * Output stream operator overloaded
    * @param si_p the SphInt you want to print
    */
   friend ostream &operator<<(ostream &output,SphInt &si_p);

   public:
      
      //constructor
      SphInt(const CartInt &);

      //copy constructor
      SphInt(const SphInt &);

      //destructor
      virtual ~SphInt();

      const Matrix &gS() const;

      Matrix &gS();

      const Matrix &gT() const;

      Matrix &gT();

      const Matrix &gU() const;

      Matrix &gU();

      const Matrix &gh() const;

      Matrix &gh();

      const Matrix &gV() const;

      Matrix &gV();

      double gS(int,int) const;

      double gT(int,int) const;

      double gU(int,int) const;

      double gh(int,int) const;

      double gV(int,int,int,int) const;

      void orthogonalize();

      static int gn_max();

      static int gN_Z();

      static int gs2inlm(int,int);

      static int ginlm2s(int,int,int,int);

      static int gg2s(int);

      static int gl_max();

      static int gdim();

      static int gN();

      static void init();

      static void clear();

   private:

      //!static objects needed to construct and destruct all the lists
      static int l_max,n_max,N_Z;

      //!nuclear repulsion energy
      static double NucRepEn;

      //!nr of particles
      static int N;

      static vector< vector<int> > s2inlm;

      static int ****inlm2s;

      //!list relating tp to sp indices
      static vector< vector<int> > t2s;

      //!list relating tp to sp indices
      static int **s2t;

      //!dimension of the basisset
      static int dim;
      
      //!overlapmatrix
      Matrix *S;

      //!kinetic energy matrix
      Matrix *T;

      //!nuclear attraction energy matrix
      Matrix *U;

      //!mean field of the "other" core
      Matrix *h;

      //!electronic repulsion matrix
      Matrix *V;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
