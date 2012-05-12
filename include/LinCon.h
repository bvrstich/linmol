#ifndef LINCON_H
#define LINCON_H

#include <iostream>
#include <cstdlib>

using std::ostream;

class TPM;

/**
 * @author Brecht Verstichel
 * @date 08-12-2010\n\n
 * This is a class that contains the information about a single linear constraint.
 */
class LinCon {

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

      TPM &gI() const;

      double gi() const;

      void sI(const TPM &);

      void si(double);

      int gM() const;

      int gN() const;

      void diag_T(int);

      void spincon(double);

      static void init(int,int);

   private:

      //!Constraint matrix
      TPM *I_c;

      //!minimal projection on the constraint matrix, such that Tr(Gamma I_C) geq i_c
      double i_c;

      static int M;

      static int N;

};

#endif
