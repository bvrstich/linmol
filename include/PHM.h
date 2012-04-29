#ifndef PHM_H
#define PHM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "BlockMatrix.h"

class TPM;
class PPHM;

/**
 * @author Brecht Verstichel
 * @date 25-04-2012\n\n
 * This class PHM is a class written for particle-hole matrices with spin and axial rotation symmetry included.
 * It inherits alle the function from its mother BlockMatrix, some special member functions 
 * and two lists that give the relationship between the sp and the ph basis.
 */
class PHM : public BlockMatrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param phm_p the PHM you want to print
    */
   friend ostream &operator<<(ostream &output,const PHM &phm_p);

   public:
      
      //constructor
      PHM();

      //copy constructor
      PHM(const PHM &);

      //destructor
      virtual ~PHM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      double operator()(int,int,int,int,int) const;

      double operator()(int,int,int,int,int,int,int,int,int,int) const;

      void G(const TPM &);

      void bar(const PPHM &);

      static void init(int,int);

      static void clear();

   private:

      //!list that relates the blockindex to physical particle-hole quantumnumbers
      static vector< vector<int> > B2SM;

      //!list that relates the physical particle-hole quantumnumbers S and M to the blockindex B
      static int **SM2B;

      //!static lists that translates the particle-particle indices to single-particle ones
      static vector< vector< vector<int> > > ph2s;

      //!inverse of the t2s list
      static int ***s2ph;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!maximal orbital angular momentum in basisset
      static int l_max;

};

#endif
