#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "BlockMatrix.h"

/**
 * @author Brecht Verstichel
 * @date 20-04-2012\n\n
 * This class TPM is a class written for two-particle matrices with spin and axial rotation symmetry included.
 * It inherits alle the function from its mother BlockMatrix, some special member functions 
 * and two lists that give the relationship between the sp and the tp basis.
 */
class TPM : public BlockMatrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param tpm_p the TPM you want to print
    */
   friend ostream &operator<<(ostream &output,const TPM &tpm_p);

   public:
      
      //constructor
      TPM();

      //copy constructor
      TPM(const TPM &);

      //destructor
      virtual ~TPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      double operator()(int B,int a,int b,int c,int d) const;

      double operator()(int,int,int,int,int,int,int,int,int,int) const;

      void molecule();

      static void init(int,int);

      static void clear();

   private:

      //!list that relates the blockindex to physical two-particle quantumnumbers
      static vector< vector<int> > B2SM;

      //!list that relates the physical two-particle quantumnumbers S and M to the blockindex B
      static int **SM2B;

      //!static lists that translates the two-particle indices to single-particle ones
      static vector< vector< vector<int> > > t2s;

      //!inverse of the t2s list
      static int ***s2t;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!maximal orbital angular momentum in basisset
      static int l_max;

};

#endif
