#ifndef DPM_H
#define DPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "BlockMatrix.h"

class TPM;

/**
 * @author Brecht Verstichel
 * @date 25-04-2012\n\n
 * This class DPM is a class written for three-particle matrices with spin and axial-rotation symmetry included.
 * It inherits alle the function from its mother BlockMatrix, some special member functions 
 * and two lists that give the relationship between the sp and the dp basis.
 */
class DPM : public BlockMatrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param dpm_p the DPM you want to print
    */
   friend ostream &operator<<(ostream &output,const DPM &dpm_p);

   public:
      
      //constructor
      DPM();

      //copy constructor
      DPM(const DPM &);

      //destructor
      virtual ~DPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      double operator()(int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int) const;

      int get_inco(int S,int Lz,int S_ab,int a,int b,int c,int *i,double *coef) const;

      //generalized T1 map
      void T(double,double,double,const TPM &);

      //maak een DPM van een TPM via de T1 conditie
      void T(const TPM &);

      //maak een DPM van een TPM via de hat functie
      void hat(const TPM &);

      static void init(int,int);

      static void clear();

   private:

      //!list that relates the blockindex to physical three-particle quantumnumbers
      static vector< vector<int> > B2SM;

      //!list that relates the physical three-particle quantumnumbers S and M to the blockindex B
      static int **SM2B;

      //!static lists that translates the three-particle indices to single-particle ones
      static vector< vector< vector<int> > > d2s;

      //!inverse of the t2s list
      static int *****s2d;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!maximal orbital angular momentum in basisset
      static int l_max;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
