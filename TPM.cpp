#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

int ***TPM::t2s;
int ***TPM::s2t;

int TPM::M;
int TPM::N;

/**
 * standard constructor for a spinsymmetrical tp matrix: constructs BlockMatrix object with 2 blocks, for S = 0 or 1,
 */
TPM::TPM() : BlockMatrix(2) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpm_c
 * if counter == 0, the lists containing the relationship between sp and tp basis.
 * @param tpm_c object that will be copied into this.
 */
TPM::TPM(const TPM &tpm_c) : BlockMatrix(tpm_c){ }

/**
 * destructor: if counter == 1 the memory for the static lists t2s en s2t will be deleted.
 * 
 */
TPM::~TPM(){ }

/**
 * @return number of particles
 */
int TPM::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int TPM::gM() const{

   return M;

}
