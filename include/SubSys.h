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

      const TPM &gsubham() const;

      TPM &gsubham();

      double gE(int) const;

      double subocc(const TPM &) const;

   private:
      
      //!hamiltonian on the subsystem
      TPM *subham;

      //!array of energies for different occupations of the subsystem
      double *E;

      //!the index of the subsystem core
      int core;

      static int M;

      static int N;

      SphInt *si;

};

#endif
