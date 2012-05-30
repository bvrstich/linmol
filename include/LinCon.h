#ifndef LINCON_H
#define LINCON_H

#include <iostream>
#include <cstdlib>

using std::ostream;

class TPM;
class SPM;
class SubSys;

/**
 * @author Brecht Verstichel
 * @date 21-05-2012\n\n
 * This is a class that contains the information about a single linear constraint.
 */

class LinCon{

   /**
    * output stream operator overloaded, will print the constraint matrix, the value of the minimal projection, and the current projection.
    * @param output The stream to which you are writing (e.g. cout)
    * @param lc_p the LinCon object you want to print
    */
   friend ostream &operator<<(ostream &output,const LinCon &lc_p);

   public:

      //constructor
      LinCon();

      //copy constructor
      LinCon(const LinCon &);

      //destructor
      virtual ~LinCon();

      const TPM &gI() const;

      const SPM &gI_bar() const;

      double gI_tr() const;

      double gi() const;

      void sI(const TPM &);

      void si(double);

      int gM() const;

      int gN() const;

      void fill_Random();

      void spincon(double);

      void subcon(const SubSys &,int);

      static void init(int,int);

   private:

      //!Traceless Constraint matrix: Watch out, shifted with i_c unity, so that Tr Gamma I_c > 0
      TPM *I_c;

      //!Partially traced, traceless constraint matrix
      SPM *I_c_bar;

      //!scaled trace of the contraint
      double I_c_tr;

      //!minimal projection on the constraint matrix, such that Tr(Gamma I_C) geq i_c
      double i_c;

      static int M;

      static int N;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/