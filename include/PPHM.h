#ifndef PPHM_H
#define PPHM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "BlockMatrix.h"

class TPM;

/**
 * @author Brecht Verstichel
 * @date 27-04-2012\n\n
 * This class PPHM is a class written for two-particle-one-hole matrices with spin and axial-rotation symmetry included.
 * It inherits alle the function from its mother BlockMatrix, some special member functions 
 * and two lists that give the relationship between the sp and the pph basis.
 */
class PPHM : public BlockMatrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param pphm_p the PPHM you want to print
    */
   friend ostream &operator<<(ostream &output,const PPHM &pphm_p);

   public:
      
      //constructor
      PPHM();

      //copy constructor
      PPHM(const PPHM &);

      //destructor
      virtual ~PPHM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      void printnax(const char *) const;

      double operator()(int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int) const;

      static int get_inco(int S,int Lz,int S_ab,int a,int b,int c,int &i);

      void T(const TPM &);

      static void init(int,int);

      static void clear();

   private:

      //!list that relates the blockindex to physical three-particle quantumnumbers
      static vector< vector<int> > B2SM;

      //!list that relates the physical three-particle quantumnumbers S and M to the blockindex B
      static int **SM2B;

      //!static lists that translates the three-particle indices to single-particle ones
      static vector< vector< vector<int> > > pph2s;

      //!inverse of the t2s list
      static int *****s2pph;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!maximal orbital angular momentum in basisset
      static int l_max;

};

#endif
