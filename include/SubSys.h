#ifndef SUBSYS_H
#define SUBSYS_H

#include <iostream>
#include <cstdlib>

using std::ostream;

class TPM;
class SPM;

/**
 * @author Brecht Verstichel
 * @date 21-05-2012\n\n
 * This is a class that contains the information about a subsystem, and an operator and energies of that subsystem
 */

class SubSys{

   /**
    * output stream operator overloaded, will print the constraint matrix, the value of the minimal projection, and the current projection.
    * @param output The stream to which you are writing (e.g. cout)
    * @param ss_p the SubSys object you want to print
    */
   friend ostream &operator<<(ostream &output,const SubSys &ss_p);

   public:

      static void init(int,int);

      //constructor
      SubSys(int,const SphInt &);

      //copy constructor
      SubSys(const SubSys &);

      //destructor
      virtual ~SubSys();

      int gcore() const;

      double gE(int,int) const;

      const vector<int> &gs2f() const;

      int gs2f(int) const;

      const Matrix &gS() const;

      Matrix &gS();

      double gW(int,int) const;

      void sW(int,int,double);

      int gn() const;

      void setB();

      void setBe();

      void setN();

      void setO();

      double subocc_func(const TPM &) const;

      const SphInt &gsi() const;

      SphInt &gsi();

      void orthogonalize();

   private:

      //!SphInt object containing the matrixelements for the subsystem
      SphInt *si;
      
      //!subsystem version of the overlap-matrix: different because the inverse in on the subsystem block ONLY!
      Matrix *S;

      //!transformation matrix between the non-orthogonal subsystem basis and the orthogonal full system on: Wath out, rectangular
      double *W;
      
      //!array of energies for different occupations of the subsystem
      double **E;

      //!the index of the subsystem core
      int core;

      //!list relating the full system index to the subsystem index (non-orthogonal basis!)
      vector<int> s2f;

      //!dimension of the subspace
      int n;

      //!nr of orbitals
      static int M;

      //!nr ofparticles
      static int N;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
