#ifndef SPM_H
#define SPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

class TPM;
class PHM;
class PPHM;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 20-04-2012\n\n
 * This class SPM was written for single-particle matrices in a system with spin and axial rotational symmetry.
 * It inherits from the class BlockMatrix and expands it with
 * specific memberfunction and a knowledge of the nr of sp orbitals and particles, and with lists translating the 
 * single-particle index to real quantumnumber information. s --> am
 */
class SPM : public BlockMatrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param spm_p de SPM you want to print
    */
   friend ostream &operator<<(ostream &output,const SPM &spm_p);

   public:
      
      //constructor
      SPM();

      //copy constructor
      SPM(const SPM &);

      //destructor
      virtual ~SPM();

      using BlockMatrix::operator=;

      int gN() const;

      int gM() const;

      void bar(double,const TPM &);

      void bar(double,const PHM &);
      
      void bar(double,const PPHM &);

      static int gl_max();

      static int gs2inl(int,int,int);

      static int ginl2s(int,int,int,int);

      static int gms2g(int,int);

      static int gg2ms(int,int);

      static void init(int,int);

      static void clear();

   private:

      //!list translating the single-particle index s and the block index m to a global single-particle index g
      static int **ms2g;

      //! and the inverse list
      static vector< vector<int> > g2ms;

      //!list translating the single-particle index to the quantumnumbers inl
      static vector< vector<int> > *s2inl;

      //!inverse list of s2am, translates two quantumnumbers a and m to the single-particle index s
      static int ****inl2s;

      //!maximal l value
      static int l_max;

      //!dimension of single particle space
      static int M;

      //!nr of particles
      static int N;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
