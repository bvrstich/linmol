/**
 * @mainpage 
 * This is an implementation of the dual only, potential reduction interior point methodfor optimizing the second order density matrix
 * for linear molecules, in which the axial and spin symmetry present in these systems is exploited. The available N-representability conditions are
 * the P, Q, G, T1 and T2 conditions.
 * Compiling is done with the options PQ, PQG, PQGT1, PQGT2 and PQGT (for all conditions active) with logical consequences for the program.
 * @author Orecht Verstichel, Ward Poelmans
 * @date 19-04-2012
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>

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
int main(int argc, char **argv){

   //initialize the random nr generator
   srand(time(NULL));

   cout.precision(10);

   CartInt::init();
   SphInt::init();

   const int M = 2*SphInt::gdim();//dim sp hilbert space
   int N = SphInt::gN();//nr of particles

   struct option long_options[] =
   {
      {"particles",  required_argument, 0, 'n'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;
   while( (j = getopt_long (argc, argv, "hn:", long_options, &i)) != -1)
      switch(j)
      {
         case 'h':
         case '?':
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
               "\n"
               "    -n, --particles=particles    Set the number of particles\n"
               "    -h, --help                   Display this help\n"
               "\n";
            return 0;
            break;
         case 'n':
            N = atoi(optarg);
            if( N <= 0)
            {
               std::cerr << "Invalid particle number!" << endl;
               return -1;
            }
            break;
      }


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

   SubSys ss_N(0,si);
   ss_N.setN();

   ss_N.orthogonalize();

   SubSys ss_O(1,si);
   ss_O.setO();

   ss_O.orthogonalize();

   LinCon::init(M,N);
   LinIneq::init(M,N,si);

   SUP::init(M,N);
   EIG::init(M,N);

   //hamiltoniaan
   si.orthogonalize();

   TPM ham;
   ham.molecule(si);

   TPM rdm;

   rdm.ReadfromFile("/home/bright/bestanden/results/linmol/NO/sub/DM_out/NO-10.h5");

   cout << N*(N - 1)/2 << "\t" << rdm.trace() << "\t" << rdm.ddot(ham) + CartInt::gNucRepEn() << endl;

   TPM subham_N;
   subham_N.subham(ss_N);

   TPM subocc_N;
   subocc_N.subocc_op(ss_N);

   TPM subham_O;
   subham_O.subham(ss_O);

   TPM subocc_O;
   subocc_O.subocc_op(ss_O);

   cout << rdm.ddot(subocc_N) << "\t" << ss_N.subocc_func(rdm) << "\t" << rdm.ddot(subham_N) << endl;
   cout << rdm.ddot(subocc_O) << "\t" << ss_O.subocc_func(rdm) << "\t" << rdm.ddot(subham_O) << endl;

   LinIneq li;
   li.fill(rdm);

   cout << li << endl;

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

/* vim: set ts=3 sw=3 expandtab :*/
