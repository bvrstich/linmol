/**
 * @mainpage 
 * This is an implementation of the dual only, potential reduction interior point methodfor optimizing the second order density matrix
 * for linear molecules, in which the axial and spin symmetry present in these systems is exploited. The available N-representability conditions are
 * the P, Q, G, T1 and T2 conditions.
 * Compiling is done with the options PQ, PQG, PQGT1, PQGT2 and PQGT (for all conditions active) with logical consequences for the program.
 * @author Brecht Verstichel, Ward Poelmans
 * @date 19-04-2012
 */

#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

//includes all important headers and defines which conditions are
//going to be used:
#include "include.h"

/**
 * In the main the actual program is run.\n 
 * We start from the unity density matrix normed on the particle number and minimize the 
 * ojective function:\n\n
 * Tr (Gamma H) - t * ln(det P(Gamma)) \n\n
 * Once the minimum is found the parameter t is reduced and a new search is initiated,
 * this goes on until convergence is reached.\n
 * The potential is minimized using the Newton-Raphson method and the resulting linear system
 * is solved via the linear conjugate gradient method.
 */
int main(void){

   //initialize the random nr generator
   srand(time(NULL));

   cout.precision(10);

   CartInt::init();
   SphInt::init();

   const int M = 2*SphInt::gdim();//dim sp hilbert space
   const int N = SphInt::gN();//nr of particles

   Tools::init(M,N);

   SPM::init(M,N);
   TPM::init(M,N);
   PHM::init(M,N);
   DPM::init(M,N);
   PPHM::init(M,N);

   SubSys::init(M,N);

   CartInt ci;
   ci.norm();

   SphInt si(ci);

   SubSys ss_Be(0,si);
   ss_Be.setBe();

   SubSys ss_B(1,si);
   ss_B.setB();

   LinCon::init(M,N);
   LinIneq::init(M,N,si);

   SUP::init(M,N);
   EIG::init(M,N);

   //hamiltoniaan
   si.orthogonalize();

   TPM ham;
   ham.molecule(si);

   ifstream in("/home/bright/bestanden/results/linmol/BeB/DM_out/BeB-20.rdm");

   TPM rdm;

   for(int B = 0;B < rdm.gnr();++B)
      for(int i = 0;i < rdm.gdim(B);++i)
         for(int j = i;j < rdm.gdim(B);++j)
            in >> B >> i >> j >> rdm(B,i,j);

   rdm.symmetrize();

   cout << N*(N - 1)/2 << "\t" << rdm.trace() << "\t" << rdm.ddot(ham) + CartInt::gNucRepEn() << endl;

   LinIneq::clear();

   PPHM::clear();
   DPM::clear();
   PHM::clear();
   TPM::clear();
   SPM::clear();

   Tools::clear();

   SphInt::clear();
   CartInt::clear();

   return 0;

}
