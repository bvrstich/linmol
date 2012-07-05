#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <hdf5.h>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::vector;

#include "include.h"

SI_SPM *SphInt::S;
SI_SPM *SphInt::T;
SI_SPM **SphInt::U;

SI_TPM *SphInt::V;

/** 
 * static function that allocates and constructs the matrixelements
 */
void SphInt::init(){

   S = new SI_SPM();
   T = new SI_SPM();

   U = new SI_SPM * [input::gN_Z()];

   for(int i = 0;i < input::gN_Z();++i)
      U[i] = new SI_SPM();

   V = new SI_TPM();

}

/** 
 * function that deallocates the static members
 */
void SphInt::clear(){

   delete S;
   delete T;

   for(int i = 0;i < input::gN_Z();++i)
      delete U[i];

   delete [] U;

   delete V;
   
}

/** 
 * @return the overlapmatrix
 */
const SI_SPM &SphInt::gS() { 

   return *S;

}

/** 
 * @return the kinetic energy matrix
 */
const SI_SPM &SphInt::gT() { 

   return *T; 

}


/** 
 * @param i index of the core
 * @return the nuclear attraction matrix of core number i
 */
const SI_SPM &SphInt::gU(int i) { 

   return *U[i]; 

}

/** 
 * @return the electronic repulsion matrix
 */
const SI_TPM &SphInt::gV() { 

   return *V; 

}

void SphInt::print(){

   cout << endl;
   cout << "Overlap Matrix:" << endl;
   cout << endl;

   cout << *S << endl;

   cout << endl;
   cout << "Kinetic energy:" << endl;
   cout << endl;

   cout << *T << endl;

   for(int i = 0;i < input::gN_Z();++i){

      cout << endl;
      cout << "Nuclear attraction energy of core " << i << endl;
      cout << endl;
      cout << *U[i] << endl;

   }

   cout << endl;
   cout << "Electronic repulsion energy:" << endl;
   cout << endl;
   cout << *V;

}

/**
 * orthogonalizes the basis: inverse sqrt of S
 */
void SphInt::orthogonalize(){

   //first inverse sqrt of S
   S->sqrt(-1);

   SI_SPM T_copy;
   SI_SPM **U_copy = new SI_SPM * [input::gN_Z()];

   for(int i = 0;i < input::gN_Z();++i)
      U_copy[i] = new SI_SPM();

   T_copy = 0.0;

   for(int i = 0;i < input::gN_Z();++i)
      *U_copy[i] = 0.0;

   //transform T
   for(int a = 0;a < SI_SPM::gdim();++a)
      for(int b = 0;b < SI_SPM::gdim();++b){

         for(int c = 0;c < SI_SPM::gdim();++c){

            T_copy(a,b) += (*S)(a,c) * (*T)(c,b);

            for(int i = 0;i < input::gN_Z();++i)
               (*U_copy[i])(a,b) += (*S)(a,c) * (*U[i])(c,b);

         }

      }

   *T = 0.0;

   for(int i = 0;i < input::gN_Z();++i)
      *U[i] = 0.0;

   for(int a = 0;a < SI_SPM::gdim();++a)
      for(int b = a;b < SI_SPM::gdim();++b){

         for(int c = 0;c < SI_SPM::gdim();++c){

            (*T)(a,b) += T_copy(a,c) * (*S)(c,b);

            for(int i = 0;i < input::gN_Z();++i)
               (*U[i])(a,b) += (*U_copy[i])(a,c) * (*S)(c,b);

         }

      }

   for(int i = 0;i < input::gN_Z();++i)
      delete U_copy[i];

   delete [] U_copy;

   SI_TPM V_copy;

   V_copy = 0.0;

   int a,b,c,d;

   //contract a
   for(int i = 0;i < SI_TPM::gtp_dim();++i){

      a = SI_TPM::gt2s(i,0);
      b = SI_TPM::gt2s(i,1);

      for(int j = 0;j < SI_TPM::gtp_dim();++j){

         c = SI_TPM::gt2s(j,0);
         d = SI_TPM::gt2s(j,1);

         for(int a_ = 0;a_ < SI_SPM::gdim();++a_)
            V_copy(i,j) += (*S)(a,a_) * (*V)(SI_TPM::gs2t(a_,b),j); 

      }
   }

   *V = 0.0;

   //contract b
   for(int i = 0;i < SI_TPM::gtp_dim();++i){

      a = SI_TPM::gt2s(i,0);
      b = SI_TPM::gt2s(i,1);

      for(int j = 0;j < SI_TPM::gtp_dim();++j){

         c = SI_TPM::gt2s(j,0);
         d = SI_TPM::gt2s(j,1);

         for(int b_ = 0;b_ < SI_SPM::gdim();++b_)
            (*V)(i,j) += (*S)(b,b_) * V_copy(SI_TPM::gs2t(a,b_),j); 

      }
   }

   V_copy = 0.0;

   //contract c
   for(int i = 0;i < SI_TPM::gtp_dim();++i){

      a = SI_TPM::gt2s(i,0);
      b = SI_TPM::gt2s(i,1);

      for(int j = 0;j < SI_TPM::gtp_dim();++j){

         c = SI_TPM::gt2s(j,0);
         d = SI_TPM::gt2s(j,1);

         for(int c_ = 0;c_ < SI_SPM::gdim();++c_)
            V_copy(i,j) += (*V)(i,SI_TPM::gs2t(c_,d)) * (*S)(c_,c); 

      }
   }

   *V = 0.0;

   //contract d
   for(int i = 0;i < SI_TPM::gtp_dim();++i){

      a = SI_TPM::gt2s(i,0);
      b = SI_TPM::gt2s(i,1);

      for(int j = 0;j < SI_TPM::gtp_dim();++j){

         c = SI_TPM::gt2s(j,0);
         d = SI_TPM::gt2s(j,1);

         for(int d_ = 0;d_ < SI_SPM::gdim();++d_)
            (*V)(i,j) += V_copy(i,SI_TPM::gs2t(c,d_)) * (*S)(d_,d); 

      }
   }

   T->symmetrize();

   for(int i = 0;i < input::gN_Z();++i)
      U[i]->symmetrize();

}

/**
 * function that transforms the cartesian integrals in CartInt to spherical ones to be stored here.
 */
void SphInt::fill(){

   S->transform(CartInt::gS());

   T->transform(CartInt::gT());

   for(int i = 0;i < input::gN_Z();++i)
      U[i]->transform(CartInt::gU(i));

   V->transform(CartInt::gV());

}

/* vim: set ts=3 sw=3 expandtab :*/
