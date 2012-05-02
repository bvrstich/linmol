#ifndef EIG_H
#define EIG_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "BlockVector.h"
#include "SUP.h"

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

/**
 * @author Brecht Verstichel
 * @date 24-04-2012\n\n
 * This class, EIG is a "block"-vector over the carrierspace's of the active condtions. It contains room
 * to store the eigenvalues and special member function that work with these eigenvalues.
 * This class should only be used when a SUP matrix has been diagonalized, some functions could give strange results when the EIG object is filled
 * with random numbers.\n\n
 */
class EIG{

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param eig_p the EIG you want to print
    */
   friend ostream &operator<<(ostream &output,const EIG &eig_p);

   public:

   //constructor met initialisatie op 
   EIG(SUP &);

   //copy constructor
   EIG(const EIG &);

   //destructor
   ~EIG();

   void diagonalize(SUP &);

   int gN() const;

   int gM() const;

   int gdim() const;

   //overload equality operator
   EIG &operator=(const EIG &);

   BlockVector<TPM> &gv_I();

   const BlockVector<TPM> &gv_I() const;

#ifdef __Q_CON
   BlockVector<TPM> &gv_Q();

   const BlockVector<TPM> &gv_Q() const;
#endif

#ifdef __G_CON
   BlockVector<PHM> &gv_G();

   const BlockVector<PHM> &gv_G() const;
#endif

#ifdef __T1_CON
   BlockVector<DPM> &gv_T1();

   const BlockVector<DPM> &gv_T1() const;
#endif

#ifdef __T2_CON
   BlockVector<PPHM> &gv_T2();

   const BlockVector<PPHM> &gv_T2() const;
#endif

   double min() const;

   double max() const;

   static void init(int,int);

   private:

   //!pointer to a BlockVector<TPM> object, the eigenvalues of the I part of a SUP matrix will be stored here.
   BlockVector<TPM> *v_I;

#ifdef __Q_CON
   //!pointer to a BlockVector<TPM> object, the eigenvalues of the Q part of a SUP matrix will be stored here.
   BlockVector<TPM> *v_Q;
#endif

#ifdef __G_CON
   //!single pointer to a BlockVector<PHM> object, the eigenvalues of G part of a SUP matrix will be stored here.
   BlockVector<PHM> *v_G;
#endif

#ifdef __T1_CON
   //!single pointer to a BlockVector<DPM> object, the eigenvalues of T1 part of a SUP matrix will be stored here.
   BlockVector<DPM> *v_T1;
#endif

#ifdef __T2_CON
   //!single pointer to a BlockVector<PPHM> object, the eigenvalues of T2 part of a SUP matrix will be stored here.
   BlockVector<PPHM> *v_T2;

#endif

   //!number of particles
   static int N;

   //!dimension of sp space
   static int M;

   //!total dimension of the EIG object
   static int dim;

};

#endif
