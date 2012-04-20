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

      static void init(int,int);

      static void clear();

   private:

      //!static list of dimension [2][dim[i]][2] that takes in a tp index i and a spinquantumnumber S, and returns two sp indices: a = t2s[S][i][0] and b = t2s[S][i][1]
      static int ***t2s;

      //!static list of dimension [2][M/2][M/2] that takes two sp indices a,b and a spinquantumnumber S, and returns a tp index i: i = s2t[S][a][b]
      static int ***s2t;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

};

#endif
