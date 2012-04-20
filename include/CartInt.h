#ifndef CARTINT_H
#define CARTINT_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"
#include "input.h"

/**
 * @author Brecht Verstichel
 * @date 12-04-2012\n\n
 * This class CartInt is a class written for the storage of matrixelements in a cartesian Gaussian basis
 * The goal is to translate ThING input to a framework where I can transform it to a spherical basis.
 */
class CartInt {

   /**
    * Output stream operator overloaded
    * @param ci_p the CartInt you want to print
    */
   friend ostream &operator<<(ostream &output,CartInt &ci_p);

   public:
      
      //constructor
      CartInt();

      //copy constructor
      CartInt(const CartInt &);

      //destructor
      virtual ~CartInt();

      const Matrix &gS() const;

      Matrix &gS();

      const Matrix &gT() const;

      Matrix &gT();

      const Matrix &gU() const;

      Matrix &gU();

      const Matrix &gV() const;

      Matrix &gV();

      void norm();

      void orthogonalize();

      double gS(int,int) const;

      double gT(int,int) const;

      double gU(int,int) const;

      double gV(int,int,int,int) const;

      static int gs2inlxyz(int,int);

      static int ginlxyz2s(int,int,int,int,int,int);

      static int gdim();

      static int gN_Z();

      static int gn_max();

      static int gl_max();

      static int gN();

      static double gNucRepEn();

      static void init();

      static void clear();

   private:

      //!static objects needed to construct and destruct all the lists
      static int l_max,n_max,N_Z;

      //!nuclear repulsion energy
      static double NucRepEn;

      //!nr of electrons
      static int N;

      //!input object contains all info about system
      static input *readin;

      //!list to switch between matrix index and physical quantum numbers
      static vector< vector<int> > s2inlxyz;

      //!list to switch between matrix index and physical quantum numbers
      static int ******inlxyz2s;

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

      //!electronic repulsion matrix
      Matrix *V;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
