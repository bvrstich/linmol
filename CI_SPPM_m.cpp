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

int **CI_SPPM_m::ab2m;

/** 
 * static function that reads in the input data and makes the corresponding lists and dimensions
 */
void CI_SPPM_m::init(){

   ab2m = new int * [CI_SPPM::gdim()];

   int i_a,j_a;
   int i_b,j_b;

   for(int a = 0;a < CI_SPPM::gdim();++a){

      i_a = CI_SPPM::gomi2ijkxyz(a,0);
      j_a = CI_SPPM::gomi2ijkxyz(a,1);

      ab2m[a] = new int [CI_SPPM::gdim()];

      for(int b = a;b < CI_SPPM::gdim();++b){

         i_b = CI_SPPM::gomi2ijkxyz(b,0);
         j_b = CI_SPPM::gomi2ijkxyz(b,1);

         ab2m[a][b] = CI_SPPM::gammax(i_a,j_a) + CI_SPPM::gammax(i_b,j_b) - CI_SPPM::gamlist(a) - CI_SPPM::gamlist(b);

      }

   }

}

/** 
 * function that deallocates the static members
 */
void CI_SPPM_m::clear(){

   for(int a = 0;a < CI_SPPM::gdim();++a)
      delete [] ab2m[a];

   delete [] ab2m;

}

/** 
 * Standard constructor
 */
CI_SPPM_m::CI_SPPM_m() {

   U = new double ** [CI_SPPM::gdim()];

   for(int a = 0;a < CI_SPPM::gdim();++a){

      U[a] = new double * [CI_SPPM::gdim() - a];

      for(int b = a;b < CI_SPPM::gdim();++b)
         U[a][b - a] = new double [ab2m[a][b] + 1];

   }

}

/** 
 * copy constructor
 * @param ci_c CI_SPPM_m object to be copied in the newly constructed object
 */
CI_SPPM_m::CI_SPPM_m(const CI_SPPM_m &ci_c) {

   U = new double ** [CI_SPPM::gdim()];

   for(int a = 0;a < CI_SPPM::gdim();++a){

      U[a] = new double * [CI_SPPM::gdim()];

      for(int b = a;b < CI_SPPM::gdim();++b)
         U[a][b - a] = new double [ab2m[a][b] + 1];

   }

   for(int a = 0;a < CI_SPPM::gdim();++a)
      for(int b = a;b < CI_SPPM::gdim();++b)
         for(int m = 0;m <= ab2m[a][b];++m)
            U[a][b - a][m] = ci_c(a,b,m);

}

/**
 * standard destructor
 */
CI_SPPM_m::~CI_SPPM_m(){

   for(int a = 0;a < CI_SPPM::gdim();++a){

      for(int b = a;b < CI_SPPM::gdim();++b)
         delete [] U[a][b - a];

      delete [] U[a];

   }

   delete [] U;

}

/**
 * access to the elements
 * @param a row vector
 * @param b column vector
 * @param m number of moments needed
 */
double CI_SPPM_m::operator()(int a,int b,int m) const {

   if(a > b)
      return U[b][a - b][m];
   else
      return U[a][b - a][m];

}

/**
 * access to the elements using quantumnumbers instead of matrixindex
 */
double CI_SPPM_m::operator()(int i_a,int j_a,int k_a,int x_a,int y_a,int z_a,int i_b,int j_b,int k_b,int x_b,int y_b,int z_b,int m) const {

   return (*this)(CI_SPPM::gijkxyz2omi(i_a,j_a,k_a,x_a,y_a,z_a),CI_SPPM::gijkxyz2omi(i_b,j_b,k_b,x_b,y_b,z_b),m);

}

/**
 * calculate the nuclear attraction energy with core "core"
 */
void CI_SPPM_m::constr_nucl_attrac(int core){

   int i_a,j_a,k_a,x_a,y_a,z_a;
   int i_b,j_b,k_b,x_b,y_b,z_b;

   for(int omi_a = 0;omi_a < CI_SPPM::gdim();++omi_a){

      i_a = CI_SPPM::gomi2ijkxyz(omi_a,0);
      j_a = CI_SPPM::gomi2ijkxyz(omi_a,1);
      k_a = CI_SPPM::gomi2ijkxyz(omi_a,2);
      x_a = CI_SPPM::gomi2ijkxyz(omi_a,3);
      y_a = CI_SPPM::gomi2ijkxyz(omi_a,4);
      z_a = CI_SPPM::gomi2ijkxyz(omi_a,5);

      for(int omi_b = omi_a;omi_b < CI_SPPM::gdim();++omi_b){

         i_b = CI_SPPM::gomi2ijkxyz(omi_b,0);
         j_b = CI_SPPM::gomi2ijkxyz(omi_b,1);
         k_b = CI_SPPM::gomi2ijkxyz(omi_b,2);
         x_b = CI_SPPM::gomi2ijkxyz(omi_b,3);
         y_b = CI_SPPM::gomi2ijkxyz(omi_b,4);
         z_b = CI_SPPM::gomi2ijkxyz(omi_b,5);

         if(i_a == i_b){//on core

            if(core == i_a){//attraction is from the same core: A = B = C

               if( (x_a + x_b)%2 != 0 || (y_a + y_b)%2 != 0 || (z_a + z_b)%2 != 0 ){

                  for(int m = 0;m <= ab2m[omi_a][omi_b];++m)
                     U[omi_a][omi_b - omi_a][m] = 0.0;

               }
               else if(CI_SPPM::gamlist(omi_a) == 0 && CI_SPPM::gamlist(omi_b) == 0){//(0_a|U(C)|0_b)

                  for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                     U[omi_a][omi_b - omi_a][m] = 2.0 * sqrt(Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b)/M_PI)

                        * Tools::gS()(i_a,j_a,k_a,0,0,0,i_b,j_b,k_b,0,0,0) / (2.0 * m + 1.0);

                  }

               }
               else{//higher angular momenta

                  //which has largest angular momentum, a or b?
                  if(CI_SPPM::gamlist(omi_b) >= CI_SPPM::gamlist(omi_a)){

                     //am_b >= am_a, find the maximal x_b,y_b or z_b
                     int flag = 0;
                     int max = x_b;

                     if(y_b > max){

                        flag = 1;
                        max = y_b;

                     }

                     if(z_b > max){

                        flag = 2;
                        max = z_b;

                     }

                     if(flag == 0){//now calculate if x_b is largest

                        if(x_b > 1){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] = (x_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 2,y_b,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 2,y_b,z_b,m + 1) );

                           }

                        }
                        else{

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m)
                              U[omi_a][omi_b - omi_a][m] = 0.0;

                        }

                        if(x_a > 0){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += x_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m + 1) );

                           }

                        }

                     }
                     else if(flag == 1){//y_b is largest

                        if(y_b > 1){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] = (y_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 2,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 2,z_b,m + 1) );

                           }

                        }
                        else{

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m)
                              U[omi_a][omi_b - omi_a][m] = 0.0;

                        }

                        if(y_a > 0){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += y_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m + 1) );

                           }

                        }

                     }
                     else{//z_b is largest

                        if(z_b > 1){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] = (z_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 2,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 2,m + 1) );

                           }

                        }
                        else{

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m)
                              U[omi_a][omi_b - omi_a][m] = 0.0;

                        }

                        if(z_a > 0){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += z_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1,m + 1) );
                           }

                        }

                     }

                  }
                  else{//am of a is larger than that of b

                     //am_a > am_b, find the maximal x_a,y_a or z_a
                     int flag = 0;
                     int max = x_a;

                     if(y_a > max){

                        flag = 1;
                        max = y_a;

                     }

                     if(z_a > max){

                        flag = 2;
                        max = z_a;

                     }

                     if(flag == 0){//now calculate if x_a is largest

                        if(x_a > 1){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] = (x_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a - 2,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a - 2,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m + 1) );

                           }

                        }
                        else{

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m)
                              U[omi_a][omi_b - omi_a][m] = 0.0;

                        }

                        if(x_b > 0){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += x_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m + 1) );

                           }

                        }

                     }
                     else if(flag == 1){//y_a is largest

                        if(y_a > 1){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] = (y_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a - 2,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a - 2,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m + 1) );

                           }

                        }
                        else{

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m)
                              U[omi_a][omi_b - omi_a][m] = 0.0;

                        }

                        if(y_b > 0){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += y_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m + 1) );

                           }

                        }

                     }
                     else{//z_a is largest

                        if(z_a > 1){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] = (z_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 2,i_b,j_b,k_b,x_b,y_b,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 2,i_b,j_b,k_b,x_b,y_b,z_b,m + 1) );

                           }

                        }
                        else{

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m)
                              U[omi_a][omi_b - omi_a][m] = 0.0;

                        }

                        if(z_b > 0){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += z_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1,m + 1) );

                           }

                        }

                     }

                  }

               }

            }
            else{//attraction is from different core: A = B but != C

               //distance between a and c
               R ac(input::gR(i_a));
               ac -= input::gR(core);

               if(CI_SPPM::gamlist(omi_a) == 0 && CI_SPPM::gamlist(omi_b) == 0){//(0_a|U(C)|0_b)

                  LibInt::calc_f(ab2m[omi_a][omi_b], Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b) * ac.ddot(ac) );

                  for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                     U[omi_a][omi_b - omi_a][m] = 2.0 * sqrt(Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b)/M_PI)

                        * Tools::gS()(i_a,j_a,k_a,0,0,0,i_b,j_b,k_b,0,0,0) * LibInt::gF(m);

                  }

               }
               else{//higher angular momenta

                  //which has largest angular momentum, a or b?
                  if(CI_SPPM::gamlist(omi_b) >= CI_SPPM::gamlist(omi_a)){

                     //am_b >= am_a, find the maximal x_b,y_b or z_b
                     int flag = 0;
                     int max = x_b;

                     if(y_b > max){

                        flag = 1;
                        max = y_b;

                     }

                     if(z_b > max){

                        flag = 2;
                        max = z_b;

                     }

                     if(flag == 0){//now calculate if x_b is largest

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m)
                           U[omi_a][omi_b - omi_a][m] = - ac[0] * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m + 1);

                        if(x_b > 1){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += (x_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 2,y_b,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 2,y_b,z_b,m + 1) );

                           }

                        }

                        if(x_a > 0){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += x_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m + 1) );

                           }

                        }

                     }
                     else if(flag == 1){//y_b is largest

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m)
                           U[omi_a][omi_b - omi_a][m] = - ac[1] * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m + 1);

                        if(y_b > 1){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += (y_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 2,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 2,z_b,m + 1) );

                           }

                        }

                        if(y_a > 0){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += y_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m + 1) );

                           }

                        }

                     }
                     else{//z_b is largest

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m)
                           U[omi_a][omi_b - omi_a][m] = - ac[2] * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 1,m + 1);

                        if(z_b > 1){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += (z_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 2,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 2,m + 1) );

                           }

                        }

                        if(z_a > 0){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += z_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1,m + 1) );
                           }

                        }

                     }

                  }
                  else{//am of a is larger than that of b

                     //am_a > am_b, find the maximal x_a,y_a or z_a
                     int flag = 0;
                     int max = x_a;

                     if(y_a > max){

                        flag = 1;
                        max = y_a;

                     }

                     if(z_a > max){

                        flag = 2;
                        max = z_a;

                     }

                     if(flag == 0){//now calculate if x_a is largest

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m)
                           U[omi_a][omi_b - omi_a][m] = - ac[0] * (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m + 1);

                        if(x_a > 1){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += (x_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a - 2,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a - 2,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m + 1) );

                           }

                        }

                        if(x_b > 0){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += x_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m + 1) );

                           }

                        }

                     }
                     else if(flag == 1){//y_a is largest

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m)
                           U[omi_a][omi_b - omi_a][m] = - ac[1] * (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m + 1);

                        if(y_a > 1){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += (y_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a - 2,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a - 2,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m + 1) );

                           }

                        }

                        if(y_b > 0){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += y_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m + 1) );

                           }

                        }

                     }
                     else{//z_a is largest

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m)
                           U[omi_a][omi_b - omi_a][m] = - ac[2] * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b,m + 1);

                        if(z_a > 1){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += (z_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 2,i_b,j_b,k_b,x_b,y_b,z_b,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 2,i_b,j_b,k_b,x_b,y_b,z_b,m + 1) );

                           }

                        }

                        if(z_b > 0){

                           for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                              U[omi_a][omi_b - omi_a][m] += z_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                                 * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1,m)

                                       - (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1,m + 1) );

                           }

                        }

                     }

                  }

               }

            }

         }
         else{//off-core

            //distance between P and c
            R pc(Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b));
            pc -= input::gR(core);

            if(CI_SPPM::gamlist(omi_a) == 0 && CI_SPPM::gamlist(omi_b) == 0){//(0_a|U(C)|0_b)

               LibInt::calc_f(ab2m[omi_a][omi_b], Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b) * pc.ddot(pc) );

               for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                  U[omi_a][omi_b - omi_a][m] = 2.0 * sqrt(Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b)/M_PI)

                     * Tools::gS()(i_a,j_a,k_a,0,0,0,i_b,j_b,k_b,0,0,0) * LibInt::gF(m);

               }

            }
            else{//higher angular momenta

               //which has largest angular momentum, a or b?
               if(CI_SPPM::gamlist(omi_b) >= CI_SPPM::gamlist(omi_a)){

                  //am_b >= am_a, find the maximal x_b,y_b or z_b
                  int flag = 0;
                  int max = x_b;

                  if(y_b > max){

                     flag = 1;
                     max = y_b;

                  }

                  if(z_b > max){

                     flag = 2;
                     max = z_b;

                  }

                  if(flag == 0){//now calculate if x_b is largest

                     for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                        U[omi_a][omi_b - omi_a][m] = ( Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[0] - input::gR(i_b)[0] ) 
                        
                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m)
                        
                           - pc[0] * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m + 1);

                     }

                     if(x_b > 1){

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                           U[omi_a][omi_b - omi_a][m] += (x_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                              * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 2,y_b,z_b,m)

                                    - (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 2,y_b,z_b,m + 1) );

                        }

                     }

                     if(x_a > 0){

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                           U[omi_a][omi_b - omi_a][m] += x_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                              * ( (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m)

                                    - (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m + 1) );

                        }

                     }

                  }
                  else if(flag == 1){//y_b is largest

                     for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                        U[omi_a][omi_b - omi_a][m] = ( Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[1] - input::gR(i_b)[1] ) 
                        
                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m)
                        
                           - pc[1] * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m + 1);

                     }

                     if(y_b > 1){

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                           U[omi_a][omi_b - omi_a][m] += (y_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                              * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 2,z_b,m)

                                    - (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 2,z_b,m + 1) );

                        }

                     }

                     if(y_a > 0){

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                           U[omi_a][omi_b - omi_a][m] += y_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                              * ( (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m)

                                    - (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m + 1) );

                        }

                     }

                  }
                  else{//z_b is largest

                     for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                        U[omi_a][omi_b - omi_a][m] = ( Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[2] - input::gR(i_b)[2] ) 
                        
                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 1,m)
                        
                           - pc[2] * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 1,m + 1);

                     }

                     if(z_b > 1){

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                           U[omi_a][omi_b - omi_a][m] += (z_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                              * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 2,m)

                                    - (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 2,m + 1) );

                        }

                     }

                     if(z_a > 0){

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                           U[omi_a][omi_b - omi_a][m] += z_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                              * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1,m)

                                    - (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1,m + 1) );
                        }

                     }

                  }

               }
               else{//am of a is larger than that of b

                  //am_a > am_b, find the maximal x_a,y_a or z_a
                  int flag = 0;
                  int max = x_a;

                  if(y_a > max){

                     flag = 1;
                     max = y_a;

                  }

                  if(z_a > max){

                     flag = 2;
                     max = z_a;

                  }

                  if(flag == 0){//now calculate if x_a is largest

                     for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                        U[omi_a][omi_b - omi_a][m] = ( Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[0] - input::gR(i_a)[0] ) 
                        
                           * (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m)
                        
                           - pc[0] * (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m + 1);

                     }

                     if(x_a > 1){

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                           U[omi_a][omi_b - omi_a][m] += (x_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                              * ( (*this)(i_a,j_a,k_a,x_a - 2,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m)

                                    - (*this)(i_a,j_a,k_a,x_a - 2,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m + 1) );

                        }

                     }

                     if(x_b > 0){

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                           U[omi_a][omi_b - omi_a][m] += x_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                              * ( (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m)

                                    - (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b,m + 1) );

                        }

                     }

                  }
                  else if(flag == 1){//y_a is largest

                     for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                        U[omi_a][omi_b - omi_a][m] = ( Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[1] - input::gR(i_a)[1] ) 

                           * (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m)

                           - pc[1] * (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m + 1);

                     }

                     if(y_a > 1){

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                           U[omi_a][omi_b - omi_a][m] += (y_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                              * ( (*this)(i_a,j_a,k_a,x_a,y_a - 2,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m)

                                    - (*this)(i_a,j_a,k_a,x_a,y_a - 2,z_a,i_b,j_b,k_b,x_b,y_b,z_b,m + 1) );

                        }

                     }

                     if(y_b > 0){

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                           U[omi_a][omi_b - omi_a][m] += y_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                              * ( (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m)

                                    - (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b,m + 1) );

                        }

                     }

                  }
                  else{//z_a is largest

                     for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                        U[omi_a][omi_b - omi_a][m] = ( Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[2] - input::gR(i_a)[2] ) 

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b,m)

                           - pc[2] * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b,m + 1);

                     }

                     if(z_a > 1){

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                           U[omi_a][omi_b - omi_a][m] += (z_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                              * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 2,i_b,j_b,k_b,x_b,y_b,z_b,m)

                                    - (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 2,i_b,j_b,k_b,x_b,y_b,z_b,m + 1) );

                        }

                     }

                     if(z_b > 0){

                        for(int m = 0;m <= ab2m[omi_a][omi_b];++m){

                           U[omi_a][omi_b - omi_a][m] += z_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                              * ( (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1,m)

                                    - (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1,m + 1) );

                        }

                     }

                  }

               }

            }

         }

      }
   }

}
