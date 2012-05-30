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

   SUP::init(M,N);
   EIG::init(M,N);

   CartInt ci;
   ci.norm();

   SphInt si(ci);
   si.orthogonalize();

   //hamiltoniaan
   TPM ham;
   ham.molecule(si);

   TPM ham_copy(ham);

   //only traceless hamiltonian needed in program.
   ham.proj_Tr();

   //primal
   SUP X;

   //dual
   SUP Z;

   //Lagrange multiplier
   SUP V;

   //just dubya
   SUP W;

   SUP u_0;

   //little help
   TPM hulp;

   u_0.gI().unit();

   u_0.fill();

   X = 0.0;
   Z = 0.0;

   //what does this do?
   double sigma = 1.0;

   double tolerance = 1.0e-5;

   double D_conv(1.0),P_conv(1.0),convergence(1.0);

   // mazziotti uses 1.6 for this
   double mazzy = 1.6;

   int iter_dual,iter_primal(0);
   int max_iter = 1;

   int tot_iter = 0;

   while(P_conv > tolerance || D_conv > tolerance || fabs(convergence) > tolerance){

      ++iter_primal;

      D_conv = 1.0;

      iter_dual = 0;

      while(D_conv > tolerance  && iter_dual <= max_iter)
      {
         ++tot_iter;

         ++iter_dual;

         //solve system
         SUP B(Z);

         B -= u_0;

         B.daxpy(mazzy/sigma,X);

         TPM b;

         b.collaps(1,B);

         b.daxpy(-mazzy/sigma,ham);

         hulp.S(-1,b);

         //hulp is the matrix containing the gamma_i's
         hulp.proj_Tr();

         //construct W
         W.fill(hulp);

         W += u_0;

         W.daxpy(-1.0/sigma,X);

         //update Z and V with eigenvalue decomposition:
         W.sep_pm(Z,V);

         V.dscal(-sigma);

         //check infeasibility of the primal problem:
         TPM v;

         v.collaps(1,V);

         v -= ham;

         D_conv = sqrt(v.ddot(v));

     }

      //update primal:
      X = V;

      //check dual feasibility (W is a helping variable now)
      W.fill(hulp);

      W += u_0;

      W -= Z;

      P_conv = sqrt(W.ddot(W));

      if(D_conv < P_conv)
         sigma *= 1.01;
      else
         sigma /= 1.01;

      convergence = ham.ddot(Z.gI()) + u_0.ddot(X);

      cout << P_conv << "\t" << D_conv << "\t" << sigma << "\t" << convergence << "\t" << ham_copy.ddot(Z.gI()) + CartInt::gNucRepEn() << endl;

   }

   cout << endl;
   cout << "Energy: " << ham_copy.ddot(Z.gI()) + CartInt::gNucRepEn() << endl;
   cout << "pd gap: " << Z.ddot(X) << endl;
   cout << "dual conv: " << D_conv << endl;
   cout << "primal conv: " << P_conv << endl;

   cout << endl;
   cout << "total nr of iterations = " << tot_iter << endl;

   Z.gI().SaveToFile("output.h5");

   ci.SaveToFile("output.h5",true);

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
