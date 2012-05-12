#ifndef SUP_H
#define SUP_H

#include <iostream> 
#include <fstream> 

using std::ostream;
using std::ofstream;
using std::ifstream;

#include "TPM.h"
#include "LinIneq.h"

#ifdef PQ

#define __Q_CON

#endif

#ifdef PQG

#define __Q_CON
#define __G_CON

#endif

#ifdef PQGT1

#define __Q_CON
#define __G_CON
#define __T1_CON

#endif

#ifdef PQGT2

#define __Q_CON
#define __G_CON
#define __T2_CON

#endif

#ifdef PQGT

#define __Q_CON
#define __G_CON
#define __T1_CON
#define __T2_CON

#endif

class EIG;

/**
 * @author Brecht Verstichel
 * @date 23-04-2012\n\n
 * This class, SUP is a blockmatrix over the carrierspace's of active N-representability conditions. 
 * This class contains two TPM objects, and if compiled with the right option a PHM, DPM or PPHM object.
 * You have to remember that these matrices are independent of each other (by which I mean that TPM::Q(SUP_PQ::tpm (0))
 * is not neccesarily equal to SUP_PQ::tpm (1)) etc. .
 */
class SUP{
  
   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param sup_p the SUP you want to print
    */
   friend ostream &operator<<(ostream &output,const SUP &sup_p);

   public:

      //constructor
      SUP();

      //copy constructor
      SUP(const SUP &);

      //destructor
      ~SUP();

      //overload += operator
      SUP &operator+=(const SUP &);

      //overload -= operator
      SUP &operator-=(const SUP &);

      //overload equality operator
      SUP &operator=(const SUP &);

      //overload equality operator
      SUP &operator=(double);

      TPM &gI();

      const TPM &gI() const;

      int gN() const;

      int gM() const;

      int gdim() const;

      double ddot(const SUP &) const;

      void invert();

      void dscal(double alpha);

      void sqrt(int option);

      void L_map(const SUP &,const SUP &);

      void daxpy(double alpha,const SUP &);

      void fill(const TPM &);

      void fill();

      void fill_Random();

#ifdef __Q_CON
      TPM &gQ();

      const TPM &gQ() const;
#endif

#ifdef __G_CON
      PHM &gG();

      const PHM &gG() const;
#endif

#ifdef __T1_CON
      DPM &gT1();

      const DPM &gT1() const;
#endif

#ifdef __T2_CON
      PPHM &gT2();

      const PPHM &gT2() const;
#endif
      
      LinIneq &gli();

      const LinIneq &gli() const;

      void sep_pm(SUP &,SUP &);

      static void init(int,int);

   private:

      //!number of sp orbitals
      static int M;

      //!nr of particles
      static int N;

      //!total dimension of the SUP matrix
      static int dim;

      //!pointer to TPM object, will contain the P block of the SUP object
      TPM *I;

#ifdef __Q_CON
      //!pointer to TPM object, will contain the Q block of the SUP object
      TPM *Q;
#endif

#ifdef __G_CON
      //!pointer to the particle hole matrix
      PHM *G;
#endif

#ifdef __T1_CON
      //!pointer tot he three particles matrix DPM
      DPM *T1;
#endif

#ifdef __T2_CON
      //!pointer tot he three particles matrix DPM
      PPHM *T2;
#endif
      
      //!the object containing the linear constraint info.
      LinIneq *li;


};

#endif
