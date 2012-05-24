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

      const vector< vector<int> > &gE() const;

      vector< vector<int> > &gE();

      void setB();

      void setBe();

      double subocc_func(const TPM &) const;

   private:
      
      //!subsystem version of the T matrix
      Matrix *T;

      //!subsystem version of the U matrix
      Matrix *U;

      //!subsystem version of the V matrix
      Matrix *V;

      //!subsystem version of the S matrix
      Matrix *S;

      //!transformation matrix between the non-orthogonal subsystem basis and the orthogonal full system on: Wath out, rectangular
      double *W;
      
      //!array of energies for different occupations of the subsystem
      vector< vector<int> > E;

      //!the index of the subsystem core
      int core;

      //!dimension of the subspace
      int n;

      //!nr of orbitals
      static int M;

      //!nr ofparticles
      static int N;

};

#endif
