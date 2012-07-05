#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_sf_coupling.h>

#include "include.h"

using std::cout;
using std::endl;

double *Tools::x6j;
double *Tools::x9j;

double **Tools::norm;
double ****Tools::gamma;

int Tools::gaussdim;

vector< vector<int> > Tools::gauss2ij;
int **Tools::ij2gauss;

vector< vector<int> > Tools::tc2i;
int **Tools::i2tc;

R ******Tools::P;
double *****Tools::prefac_overlap;

double *Tools::dfac;
double *Tools::fac;
double **Tools::comb;

CI_SPPM *Tools::S;
CI_SPPM *Tools::T;

CI_SPPM_m **Tools::U;

CI_TPPM *Tools::V;

/**
 * initialize and allocate the static variables
 */
void Tools::init(){

   x6j = new double [16];

   x9j = new double [16];

   //fill the x6j list
   for(int S_ = 0;S_ < 2;++S_)
      for(int S = 0;S < 2;++S)
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd)
               x6j[8*S_ + 4*S + 2*S_ab + S_cd] = gsl_sf_coupling_6j(2*S_ + 1,1,2*S_ab,2*S + 1,1,2*S_cd);

   //fill the x9j list
   for(int S_ = 0;S_ < 2;++S_)
      for(int Z = 0;Z < 2;++Z)
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd)
               x9j[8*S_ + 4*Z + 2*S_ab + S_cd] = gsl_sf_coupling_9j(1,1,2*Z,2*S_ + 1,2*S_ab,1,2*S_cd,1,1);

   //construct the double faculty list
   dfac = new double [2 * LibInt::gMAXFAC()];

   dfac[0] = 1.0;
   dfac[1] = 1.0;
   dfac[2] = 1.0;

   for (int i = 3;i < LibInt::gMAXFAC() * 2; i++) 
      dfac[i] = (i - 1) * dfac[i - 2];

   //allocate and construct the faculty list
   fac = new double [LibInt::gMAXFAC()];

   fac[0] = 1.0;
   fac[1] = 1.0;

   for(int i = 2;i < LibInt::gMAXFAC();++i)
      fac[i] = i * fac[i - 1];

   //allocate and construct the combination array
   comb = new double * [LibInt::gL_MAX() + 1];

   for(int i = 0;i <= LibInt::gL_MAX();++i)
      comb[i] = new double [LibInt::gL_MAX() + 1];

   for(int n = 0;n <= LibInt::gL_MAX();++n)
      for(int k = 0;k <= n;++k)
         comb[k][n] = fac[n]/(fac[k]*fac[n - k]);

   //make the norms of the primitives
   norm = new double * [CI_SPM::gdim()];

   int i,x,y,z;
   int j;

   for(int s = 0;s < CI_SPM::gdim();++s){

      i = CI_SPM::gs2inlxyz(s,0);
      j = CI_SPM::gs2gauss(s,1);

      x = CI_SPM::gs2inlxyz(s,3);
      y = CI_SPM::gs2inlxyz(s,4);
      z = CI_SPM::gs2inlxyz(s,5);

      norm[s] = new double [input::gGaussInfo(i).gNcontr(j)];

      for(int n_p = 0;n_p < input::gGaussInfo(i).gNcontr(j);++n_p)
         norm[s][n_p] = Tools::norm_const(x,y,z,input::gGaussInfo(i).galpha(j,n_p));

   }

   vector<int> v(2);

   //what is the gaussian dimension?
   gaussdim = 0;

   for(int i = 0;i < input::gN_Z();++i){

      for(int j = 0;j < input::gGaussInfo(i).gNtypes();++j){

         v[0] = i;
         v[1] = j;

         gauss2ij.push_back(v);

         gaussdim++;

      }

   }

   //make the reverse list
   ij2gauss = new int * [input::gN_Z()];

   for(int i = 0;i < input::gN_Z();++i)
      ij2gauss[i] = new int [input::gGaussInfo(i).gNtypes()];

   for(int g = 0;g < gaussdim;++g)
      ij2gauss[gauss2ij[g][0]][gauss2ij[g][1]] = g;

   //make the gamma matrix
   gamma = new double *** [gaussdim];

   for(int ga = 0;ga < gaussdim;++ga){

      int i_a = gauss2ij[ga][0];
      int j_a = gauss2ij[ga][1];

      gamma[ga] = new double ** [gaussdim];

      //per "Gaussian Matrix" element construct a rectangular matrix of the primitives
      for(int gb = 0;gb < gaussdim;++gb){

         int i_b = gauss2ij[gb][0];
         int j_b = gauss2ij[gb][1];

         int n_p_a = input::gGaussInfo(i_a).gNcontr(j_a);

         gamma[ga][gb] = new double * [n_p_a];

         for(int k = 0;k < n_p_a;++k)
            gamma[ga][gb][k] = new double [input::gGaussInfo(i_b).gNcontr(j_b)];

      }

   }

   //fill gamma:
   for(int ga = 0;ga < gaussdim;++ga){

      int i_a = gauss2ij[ga][0];
      int j_a = gauss2ij[ga][1];

      int n_p_a = input::gGaussInfo(i_a).gNcontr(j_a);

      //per "Gaussian Matrix" element construct a rectangular matrix of the primitives
      for(int gb = 0;gb < gaussdim;++gb){

         int i_b = gauss2ij[gb][0];
         int j_b = gauss2ij[gb][1];

         int n_p_b = input::gGaussInfo(i_b).gNcontr(j_b);

         for(int k_a = 0;k_a < n_p_a;++k_a)
            for(int k_b = 0;k_b < n_p_b;++k_b)
               gamma[ga][gb][k_a][k_b] = input::gGaussInfo(i_a).galpha(j_a,k_a) + input::gGaussInfo(i_b).galpha(j_b,k_b);

      }
   }

   //construct some more lists, relating the off-diagonal cores to a larger "two-core" index
   i2tc = new int * [input::gN_Z()];

   for(int i = 0;i < input::gN_Z();++i)
      i2tc[i] = new int [input::gN_Z()];

   int iter = 0;

   for(int i = 0;i < input::gN_Z();++i)
      for(int j = i + 1;j < input::gN_Z();++j){

         v[0] = i;
         v[1] = j;

         tc2i.push_back(v);

         i2tc[i][j] = iter;
         i2tc[j][i] = iter;

         ++iter;

      }

   //construct the P object
   P = new R ***** [tc2i.size()];

   for(unsigned int tc = 0;tc < tc2i.size();++tc){

      int i_a = tc2i[tc][0];
      int i_b = tc2i[tc][1];

      P[tc] = new R **** [input::gGaussInfo(i_a).gNtypes()];

      for(int j_a = 0;j_a < input::gGaussInfo(i_a).gNtypes();++j_a){

         int n_p_a = input::gGaussInfo(i_a).gNcontr(j_a);

         P[tc][j_a] = new R *** [input::gGaussInfo(i_b).gNtypes()];

         for(int j_b = 0;j_b < input::gGaussInfo(i_b).gNtypes();++j_b){

            int n_p_b = input::gGaussInfo(i_b).gNcontr(j_b);

            P[tc][j_a][j_b] = new R ** [n_p_a];

            for(int k_a = 0;k_a < n_p_a;++k_a){

               P[tc][j_a][j_b][k_a] = new R * [input::gGaussInfo(i_b).gNcontr(j_b)];

               for(int k_b = 0;k_b < n_p_b;++k_b){

                  double x_P = input::gGaussInfo(i_a).galpha(j_a,k_a) * input::gR(i_a)[0] + input::gGaussInfo(i_b).galpha(j_b,k_b) * input::gR(i_b)[0];
                  double y_P = input::gGaussInfo(i_a).galpha(j_a,k_a) * input::gR(i_a)[1] + input::gGaussInfo(i_b).galpha(j_b,k_b) * input::gR(i_b)[1];
                  double z_P = input::gGaussInfo(i_a).galpha(j_a,k_a) * input::gR(i_a)[2] + input::gGaussInfo(i_b).galpha(j_b,k_b) * input::gR(i_b)[2];

                  x_P /= ggamma(i_a,j_a,i_b,j_b,k_a,k_b);
                  y_P /= ggamma(i_a,j_a,i_b,j_b,k_a,k_b);
                  z_P /= ggamma(i_a,j_a,i_b,j_b,k_a,k_b);

                  P[tc][j_a][j_b][k_a][k_b] = new R(x_P,y_P,z_P);

               }

            }

         }

      }

   }

   //construct the prefac_overlap
   prefac_overlap = new double **** [tc2i.size()];

   for(unsigned int tc = 0;tc < tc2i.size();++tc){

      int i_a = tc2i[tc][0];
      int i_b = tc2i[tc][1];

      R ab(input::gR(i_a));

      ab -= input::gR(i_b);

      double ab_dist = ab.ddot(ab);

      prefac_overlap[tc] = new double *** [input::gGaussInfo(i_a).gNtypes()];

      for(int j_a = 0;j_a < input::gGaussInfo(i_a).gNtypes();++j_a){

         int n_p_a = input::gGaussInfo(i_a).gNcontr(j_a);

         prefac_overlap[tc][j_a] = new double ** [input::gGaussInfo(i_b).gNtypes()];

         for(int j_b = 0;j_b < input::gGaussInfo(i_b).gNtypes();++j_b){

            int n_p_b = input::gGaussInfo(i_b).gNcontr(j_b);

            prefac_overlap[tc][j_a][j_b] = new double * [n_p_a];

            for(int k_a = 0;k_a < n_p_a;++k_a){

               prefac_overlap[tc][j_a][j_b][k_a] = new double [n_p_b];

               for(int k_b = 0;k_b < n_p_b;++k_b){

                  double tmp = -input::gGaussInfo(i_a).galpha(j_a,k_a) * input::gGaussInfo(i_b).galpha(j_b,k_b) * ab_dist / ggamma(i_a,j_a,i_b,j_b,k_a,k_b);

                  prefac_overlap[tc][j_a][j_b][k_a][k_b] = exp( tmp );

               }
            }
         }
      }
   }

   S = new CI_SPPM();
   T = new CI_SPPM();

   S->constr_overlap();
   T->constr_kinetic();

   U = new CI_SPPM_m * [input::gN_Z()];

   for(int i = 0;i < input::gN_Z();++i){

      U[i] = new CI_SPPM_m();

      U[i]->constr_nucl_attrac(i);

   }

   V = new CI_TPPM();

   V->constr_eri();

}

/**
 * deallocate the static lists and objects
 */
void Tools::clear(){

   delete [] x6j;

   delete [] x9j;

      for(int s = 0;s < CI_SPM::gdim();++s)
      delete [] norm[s];

   delete [] norm;

   //delete the gamma matrix
   for(int ga = 0;ga < gaussdim;++ga){

      int i_a = gauss2ij[ga][0];
      int j_a = gauss2ij[ga][1];

      for(int gb = 0;gb < gaussdim;++gb){

         for(int k = 0;k < input::gGaussInfo(i_a).gNcontr(j_a);++k)
            delete [] gamma[ga][gb][k];

         delete [] gamma[ga][gb];

      }

      delete [] gamma[ga];

   }

   delete [] gamma;

   for(int i = 0;i < input::gN_Z();++i)
      delete [] ij2gauss[i];

   delete [] ij2gauss;

   for(int i = 0;i < input::gN_Z();++i)
      delete [] i2tc[i];

   delete [] i2tc;

   //delete the P object
   for(unsigned int tc = 0;tc < tc2i.size();++tc){

      int i_a = tc2i[tc][0];
      int i_b = tc2i[tc][1];

      for(int j_a = 0;j_a < input::gGaussInfo(i_a).gNtypes();++j_a){

         int n_p_a = input::gGaussInfo(i_a).gNcontr(j_a);

         for(int j_b = 0;j_b < input::gGaussInfo(i_b).gNtypes();++j_b){

            int n_p_b = input::gGaussInfo(i_b).gNcontr(j_b);

            for(int k_a = 0;k_a < n_p_a;++k_a){

               for(int k_b = 0;k_b < n_p_b;++k_b)
                  delete P[tc][j_a][j_b][k_a][k_b];

               delete [] P[tc][j_a][j_b][k_a];

            }

            delete [] P[tc][j_a][j_b];

         }

         delete [] P[tc][j_a];

      }

      delete [] P[tc];

   }

   delete [] P;

   for(unsigned int tc = 0;tc < tc2i.size();++tc){

      int i_a = tc2i[tc][0];
      int i_b = tc2i[tc][1];

      for(int j_a = 0;j_a < input::gGaussInfo(i_a).gNtypes();++j_a){

         int n_p_a = input::gGaussInfo(i_a).gNcontr(j_a);

         for(int j_b = 0;j_b < input::gGaussInfo(i_b).gNtypes();++j_b){

            for(int k_a = 0;k_a < n_p_a;++k_a)
               delete [] prefac_overlap[tc][j_a][j_b][k_a];

            delete [] prefac_overlap[tc][j_a][j_b];

         }

         delete [] prefac_overlap[tc][j_a];

      }

      delete [] prefac_overlap[tc];

   }

   delete [] prefac_overlap;

   delete [] dfac;

   delete [] fac;

   for(int i = 0;i <= LibInt::gL_MAX();++i)
      delete [] comb[i];

   delete [] comb;

   delete S;
   delete T;

   for(int i = 0;i < input::gN_Z();++i)
      delete U[i];

   delete [] U;

   delete V;

}

/**
 * store only the six j symbols needed
 * @param S_ 0 or 1 means 1/2 or 3/2
 * @param S 0 or 1 means 1/2 or 3/2
 * @param S_ab 0 or 1 means 0 or 1
 * @param S_cd 0 or 1 means 0 or 1
 * @return the 6j symbol
 */
double Tools::g6j(int S_,int S,int S_ab,int S_cd){

   return x6j[8*S_ + 4*S + 2*S_ab + S_cd];

}

/**
 * store only the nine j symbols needed
 * @param S 0 or 1 means 1/2 or 3/2
 * @param Z 0 or 1 means 0 or 1
 * @param S_ab 0 or 1 means 0 or 1
 * @param S_cd 0 or 1 means 0 or 1
 * @return the 6j symbol
 */
double Tools::g9j(int S,int Z,int S_ab,int S_cd){

   return x9j[8*S + 4*Z + 2*S_ab + S_cd];

}

/**
 * @return the norm of the primitive with index k corresponding to the gaussian contraction with CI_SPM index s
 */
double Tools::gnorm(int s,int k){

   return norm[s][k];

}

/**
 * @param i_a location of the first core
 * @param j_a gaussian index of the orbital on the first core
 * @param i_b location of the second core
 * @param j_b gaussian index of the orbital on the second core 
 * @param k_a index of the primitive belonging to A
 * @param k_b index of the primitive belonging to B
 * @return the gamma corresponding to the input parameters
 */
double Tools::ggamma(int i_a,int j_a,int i_b,int j_b,int k_a,int k_b){

   int ga = ij2gauss[i_a][j_a];
   int gb = ij2gauss[i_b][j_b];

   return gamma[ga][gb][k_a][k_b];

}

/**
 * @param i_a location of the first core
 * @param j_a gaussian index of the orbital on the first core
 * @param i_b location of the second core
 * @param j_b gaussian index of the orbital on the second core 
 * @param k_a index of the primitive belonging to A
 * @param k_b index of the primitive belonging to B
 * @return the P-vector corresponding to the input parameters
 */
const R &Tools::gP(int i_a,int j_a,int i_b,int j_b,int k_a,int k_b){

   if(i_a == i_b)
      cout << "What the fuck are you doing man!" << endl;

   int tc = i2tc[i_a][i_b];

   return *P[tc][j_a][j_b][k_a][k_b];

}

/**
 * @param i_a location of the first core
 * @param j_a gaussian index of the orbital on the first core
 * @param i_b location of the second core
 * @param j_b gaussian index of the orbital on the second core 
 * @param k_a index of the primitive belonging to A
 * @param k_b index of the primitive belonging to B
 * @return the prefactor appearing in every overlap integral: ONLY off-diagonal integrals!
 */
double Tools::gprefac_overlap(int i_a,int j_a,int i_b,int j_b,int k_a,int k_b){

   if(i_a == i_b)
      cout << "What the fuck are you doing man!" << endl;

   int tc = i2tc[i_a][i_b];

   return prefac_overlap[tc][j_a][j_b][k_a][k_b];

}

/**
 * get the double faculty:
 * @return gdfac(2n) gives back (2n - 1)!!
 */
double Tools::gdfac(int n){

   return dfac[n];

}

/**
 * get the faculty: 
 * @return n!
 */
double Tools::gfac(int n){

   return fac[n];

}

/**
 * get the combinations
 * @return n!/k!(n-k)!
 */
double Tools::gcomb(int k,int n){

   return comb[k][n];

}

/**
 * @return the norm of a single primitive with input nx,ny,nz and alpha
 */
double Tools::norm_const(int nx,int ny,int nz, double alpha) {

   return pow(2 * alpha / M_PI, 0.75) * pow(4 * alpha, 0.5 * (nx + ny + nz)) / sqrt(Tools::gdfac(2 * nx) * Tools::gdfac(2 * ny) * Tools::gdfac(2 * nz));

}

/**
 * access to the overlapmatrix of the primitives
 */
const CI_SPPM &Tools::gS(){

   return *S;

}

/**
 * access to the kinetic energy matrix of the primitives
 */
const CI_SPPM &Tools::gT(){

   return *T;

}

/**
 * @param core the index of the core
 * access to the nuclear attraction energy matrix of the primitives
 */
const CI_SPPM_m &Tools::gU(int core){

   return *U[core];

}

/**
 * access to the the eri object of the primitives
 */
const CI_TPPM &Tools::gV(){

   return *V;

}

/* vim: set ts=3 sw=3 expandtab :*/
