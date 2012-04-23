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

   SUP::init(M,N);

   CartInt ci;
   ci.norm();

   SphInt si(ci);
   si.orthogonalize();

   //hamiltoniaan
   TPM ham;
   ham.molecule(si);

/*
   TPM rdm(M,N);
   rdm.unit();

   TPM backup_rdm(rdm);

   double t = 1.0;
   double tolerance = 1.0e-5;

   //outer iteration: scaling of the potential barrier
   while(t > 1.0e-12){

      cout << t << "\t" << rdm.trace() << "\t" << rdm.ddot(ham) << endl;

      double convergence = 1.0;

      //inner iteration: 
      //Newton's method for finding the minimum of the current potential
      while(convergence > tolerance){

         SUP P(M,N);

         P.fill(rdm);

         P.invert();

         //eerst -gradient aanmaken:
         TPM grad(M,N);

         grad.constr_grad(t,ham,P);

         //dit wordt de stap:
         TPM delta(M,N);

         //los het hessiaan stelsel op:
         cout << delta.solve(t,P,grad) << endl;

         //line search
         double a = delta.line_search(t,P,ham);

         //rdm += a*delta;
         rdm.daxpy(a,delta);

         convergence = a*a*delta.ddot(delta);

      }

      t /= 2.0;

      //what is the tolerance for the newton method?
      tolerance = 1.0e-5*t;

      if(tolerance < 1.0e-12)
         tolerance = 1.0e-12;

      //extrapolatie:
      TPM extrapol(rdm);

      extrapol -= backup_rdm;

      //overzetten voor volgende stap
      backup_rdm = rdm;

      double a = extrapol.line_search(t,rdm,ham);

      rdm.daxpy(a,extrapol);

   }
*/
   TPM::clear();
   SPM::clear();

   Tools::clear();
   
   SphInt::clear();
   CartInt::clear();

   return 0;
}
