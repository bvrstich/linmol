#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>

using std::ostream;

class SphInt;
class SUP;
class DPM;
class PHM;
class PPHM;

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

      double operator()(int,int,int,int,int) const;

      double operator()(int,int,int,int,int,int,int,int,int,int) const;

      void molecule(const SphInt &);

      void unit();

      void proj_Tr();
      
      void min_unit(double scale);

      //Q afbeelding en zijn inverse
      void Q(int option,const TPM &);

      //Q like afbeelding Q(A,B,C,tpm_d)
      void Q(int option,double A,double B,double C,const TPM &);

      //G down
      void G(const PHM &);

      void bar(const DPM &dpm);

      //T1 down
      void T(const DPM &);

      void bar(const PPHM &);

      void T(const PPHM &);

      void printnax(const char *) const;

      void collaps(int,const SUP &);

      //overlapmatrix afbeelding en zijn inverse
      void S(int,const TPM &);

      static void init_overlap();

      static void init(int,int);

      static void clear();

      int SaveToFile(const char *);

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

      //!overlapmatrix parameters
      static double Sa,Sc;

      //!dimension of sp hilbert space
      static int M;

      //!maximal orbital angular momentum in basisset
      static int l_max;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
