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

CI_SPM *CartInt::S;
CI_SPM *CartInt::T;
CI_SPM **CartInt::U;

CI_TPM *CartInt::V;

/** 
 * static function that allocates and constructs the matrixelements
 */
void CartInt::init(){

   S = new CI_SPM();
   T = new CI_SPM();

   U = new CI_SPM * [input::gN_Z()];

   for(int i = 0;i < input::gN_Z();++i)
      U[i] = new CI_SPM();

   V = new CI_TPM();

}

/** 
 * function that deallocates the static members
 */
void CartInt::clear(){

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
const CI_SPM &CartInt::gS() { 

   return *S;

}

/** 
 * @return the kinetic energy matrix
 */
const CI_SPM &CartInt::gT() { 

   return *T; 

}


/** 
 * @param i index of the core
 * @return the nuclear attraction matrix of core number i
 */
const CI_SPM &CartInt::gU(int i) { 

   return *U[i]; 

}

/** 
 * @return the electronic repulsion matrix
 */
const CI_TPM &CartInt::gV() { 

   return *V; 

}

/**
 * normalize the cartesian wavefunctions, needed for transformation to spherical
 */
void CartInt::norm() {

   double norm[CI_SPM::gdim()];

   for(int a = 0;a < CI_SPM::gdim();++a)
      norm[a] = sqrt( (*S)(a,a) );

   //first the one body matrices
   for(int a = 0;a < CI_SPM::gdim();++a)
      for(int b = a;b < CI_SPM::gdim();++b){

         (*S)(a,b) /= norm[a] * norm[b];
         (*T)(a,b) /= norm[a] * norm[b];

         for(int i = 0;i < input::gN_Z();++i)
            (*U[i])(a,b) /= norm[a] * norm[b];

      }

   int a,b,c,d;

   //then the two body matrices
   for(int i = 0;i < CI_TPM::gtp_dim();++i){

      a = CI_TPM::gt2s(i,0);
      b = CI_TPM::gt2s(i,1);

      for(int j = i;j < CI_TPM::gtp_dim();++j){

         c = CI_TPM::gt2s(j,0);
         d = CI_TPM::gt2s(j,1);

         (*V)(i,j) /= norm[a] * norm[b] * norm[c] * norm[d];

      }

   }

   S->symmetrize();
   T->symmetrize();

   for(int i = 0;i < input::gN_Z();++i)
      U[i]->symmetrize();

   V->symmetrize();

}

void CartInt::print(){

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
void CartInt::orthogonalize() {

   //first inverse sqrt of S
   S->sqrt(-1);

   CI_SPM T_copy;
   CI_SPM **U_copy = new CI_SPM * [input::gN_Z()];

   for(int i = 0;i < input::gN_Z();++i)
      U_copy[i] = new CI_SPM();

   T_copy = 0.0;

   for(int i = 0;i < input::gN_Z();++i)
      *U_copy[i] = 0.0;

   //transform T
   for(int a = 0;a < CI_SPM::gdim();++a)
      for(int b = 0;b < CI_SPM::gdim();++b){

         for(int c = 0;c < CI_SPM::gdim();++c){

            T_copy(a,b) += (*S)(a,c) * (*T)(c,b);

            for(int i = 0;i < input::gN_Z();++i)
               (*U_copy[i])(a,b) += (*S)(a,c) * (*U[i])(c,b);

         }

      }

   *T = 0.0;

   for(int i = 0;i < input::gN_Z();++i)
      *U[i] = 0.0;

   for(int a = 0;a < CI_SPM::gdim();++a)
      for(int b = a;b < CI_SPM::gdim();++b){

         for(int c = 0;c < CI_SPM::gdim();++c){

            (*T)(a,b) += T_copy(a,c) * (*S)(c,b);

            for(int i = 0;i < input::gN_Z();++i)
               (*U[i])(a,b) += (*U_copy[i])(a,c) * (*S)(c,b);

         }

      }

   for(int i = 0;i < input::gN_Z();++i)
      delete U_copy[i];

   delete [] U_copy;

   CI_TPM V_copy;

   V_copy = 0.0;

   int a,b,c,d;

   //contract a
   for(int i = 0;i < CI_TPM::gtp_dim();++i){

      a = CI_TPM::gt2s(i,0);
      b = CI_TPM::gt2s(i,1);

      for(int j = 0;j < CI_TPM::gtp_dim();++j){

         c = CI_TPM::gt2s(j,0);
         d = CI_TPM::gt2s(j,1);

         for(int a_ = 0;a_ < CI_SPM::gdim();++a_)
            V_copy(i,j) += (*S)(a,a_) * (*V)(CI_TPM::gs2t(a_,b),j); 

      }
   }

   *V = 0.0;

   //contract b
   for(int i = 0;i < CI_TPM::gtp_dim();++i){

      a = CI_TPM::gt2s(i,0);
      b = CI_TPM::gt2s(i,1);

      for(int j = 0;j < CI_TPM::gtp_dim();++j){

         c = CI_TPM::gt2s(j,0);
         d = CI_TPM::gt2s(j,1);

         for(int b_ = 0;b_ < CI_SPM::gdim();++b_)
            (*V)(i,j) += (*S)(b,b_) * V_copy(CI_TPM::gs2t(a,b_),j); 

      }
   }

   V_copy = 0.0;

   //contract c
   for(int i = 0;i < CI_TPM::gtp_dim();++i){

      a = CI_TPM::gt2s(i,0);
      b = CI_TPM::gt2s(i,1);

      for(int j = 0;j < CI_TPM::gtp_dim();++j){

         c = CI_TPM::gt2s(j,0);
         d = CI_TPM::gt2s(j,1);

         for(int c_ = 0;c_ < CI_SPM::gdim();++c_)
            V_copy(i,j) += (*V)(i,CI_TPM::gs2t(c_,d)) * (*S)(c_,c); 

      }
   }

   *V = 0.0;

   //contract d
   for(int i = 0;i < CI_TPM::gtp_dim();++i){

      a = CI_TPM::gt2s(i,0);
      b = CI_TPM::gt2s(i,1);

      for(int j = 0;j < CI_TPM::gtp_dim();++j){

         c = CI_TPM::gt2s(j,0);
         d = CI_TPM::gt2s(j,1);

         for(int d_ = 0;d_ < CI_SPM::gdim();++d_)
            (*V)(i,j) += V_copy(i,CI_TPM::gs2t(c,d_)) * (*S)(d_,d); 

      }
   }

   T->symmetrize();

   for(int i = 0;i < input::gN_Z();++i)
      U[i]->symmetrize();

}

/**
 * function that calculates the integrals
 */
void CartInt::calc_integrals(){

   S->calc_overlap();

   T->calc_kinetic();

   for(int i = 0;i < input::gN_Z();++i)
      U[i]->calc_nucl_attrac(i);

   V->calc_eri();

}

/* vim: set ts=3 sw=3 expandtab :*/
