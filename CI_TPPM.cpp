#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <hdf5.h>
#include <libint2.h>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::vector;

#include "include.h"

vector< vector<int> > CI_TPPM::tomi2omi;
int **CI_TPPM::omi2tomi;

vector< vector<int> > CI_TPPM::g2ijk;
int ***CI_TPPM::ijk2g;

vector< vector<int> > CI_TPPM::t2s;
int **CI_TPPM::s2t;

/** 
 * static function that reads in the input data and makes the corresponding lists and dimensions
 */
void CI_TPPM::init(){

   omi2tomi = new int * [CI_SPPM::gdim()];

   for(int omi_a = 0;omi_a < CI_SPPM::gdim();++omi_a)
      omi2tomi[omi_a] = new int [CI_SPPM::gdim()];

   vector<int> v(2);

   int iter = 0;

   for(int omi_a = 0;omi_a < CI_SPPM::gdim();++omi_a)
      for(int omi_b = omi_a;omi_b < CI_SPPM::gdim();++omi_b){

         v[0] = omi_a;
         v[1] = omi_b;

         tomi2omi.push_back(v);

         omi2tomi[omi_a][omi_b] = iter;
         omi2tomi[omi_b][omi_a] = iter;

         ++iter;

      }

   vector<int> v_g(3);

   ijk2g = new int ** [input::gN_Z()];

   iter = 0;

   for(int i = 0;i < input::gN_Z();++i){

      ijk2g[i] = new int * [input::gGaussInfo(i).gNtypes()];

      for(int j = 0;j < input::gGaussInfo(i).gNtypes();++j){

         ijk2g[i][j] = new int [input::gGaussInfo(i).gNcontr(j)];

         for(int k = 0;k < input::gGaussInfo(i).gNcontr(j);++k){

            v_g[0] = i;
            v_g[1] = j;
            v_g[2] = k;

            g2ijk.push_back(v_g);

            ijk2g[i][j][k] = iter;

            ++iter;

         }
      }
   }

   s2t = new int * [g2ijk.size()];

   for(unsigned int ga = 0;ga < g2ijk.size();++ga)
      s2t[ga] = new int [g2ijk.size()];

   iter = 0;

   for(unsigned int ga = 0;ga < g2ijk.size();++ga)
      for(unsigned int gb = ga;gb < g2ijk.size();++gb){

         v[0] = ga;
         v[1] = gb;

         t2s.push_back(v);

         s2t[ga][gb] = iter;
         s2t[gb][ga] = iter;

         ++iter;

      }

}

/** 
 * function that deallocates the static members
 */
void CI_TPPM::clear(){

   for(int omi_a = 0;omi_a < CI_SPPM::gdim();++omi_a)
      delete [] omi2tomi[omi_a];

   delete [] omi2tomi;

   for(unsigned int ga = 0;ga < g2ijk.size();++ga)
      delete [] s2t[ga];

   delete [] s2t;

   for(int i = 0;i < input::gN_Z();++i){

      for(int j = 0;j < input::gGaussInfo(i).gNtypes();++j)
         delete [] ijk2g[i][j];

      delete [] ijk2g[i];

   }

   delete [] ijk2g;

}

/** 
 * Standard constructor: allocates the matrices and fills them with the correct elements
 */
CI_TPPM::CI_TPPM() : Matrix(tomi2omi.size()) { }

/** 
 * copy constructor
 * @param ci_c CI_TPPM object to be copied in the newly constructed object
 */
CI_TPPM::CI_TPPM(const CI_TPPM &ci_c) : Matrix(ci_c){ }

/**
 * standard destructor
 */
CI_TPPM::~CI_TPPM(){ }

/**
 * fill the CI_TPPM object with the electron repulsion integrals
 */
void CI_TPPM::constr_eri(){

   LIBINT2_PREFIXED_NAME(libint2_static_init)();

   Libint_t inteval;

   LIBINT2_PREFIXED_NAME(libint2_init_eri)(&inteval,CI_SPM::gl_max(),0);

   int ga,gb,gc,gd;
   int l_a,l_b,l_c,l_d;

   int i_a,j_a,k_a;
   int i_b,j_b,k_b;
   int i_c,j_c,k_c;
   int i_d,j_d,k_d;

   //warning! chemical notation! (a b | V | c d ) = < a c | V | b d >
   for(unsigned int ti = 0;ti < t2s.size();++ti){

      ga = t2s[ti][0];
      gb = t2s[ti][1];

      i_a = g2ijk[ga][0];
      j_a = g2ijk[ga][1];
      k_a = g2ijk[ga][2];

      i_b = g2ijk[gb][0];
      j_b = g2ijk[gb][1];
      k_b = g2ijk[gb][2];

      for(unsigned int tj = ti;tj < t2s.size();++tj){

         gc = t2s[tj][0];
         gd = t2s[tj][1];

         i_c = g2ijk[gc][0];
         j_c = g2ijk[gc][1];
         k_c = g2ijk[gc][2];

         i_d = g2ijk[gd][0];
         j_d = g2ijk[gd][1];
         k_d = g2ijk[gd][2];

         l_a = CI_SPPM::gammax(i_a,j_a);
         l_b = CI_SPPM::gammax(i_b,j_b);
         l_c = CI_SPPM::gammax(i_c,j_c);
         l_d = CI_SPPM::gammax(i_d,j_d);

         if(l_a == 0 && l_b == 0 && l_c == 0 && l_d == 0){ //(ss|ss) do it yourself

            if(i_a == i_b && i_b == i_c && i_c == i_d){//same core

               (*this)(i_a,j_a,k_a,0,0,0,i_b,j_b,k_b,0,0,0,i_c,j_c,k_c,0,0,0,i_d,j_d,k_d,0,0,0) = 

                  2.0 * pow(M_PI,2.5)/( Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b) * Tools::ggamma(i_c,j_c,i_d,j_d,k_c,k_d) 

                        * std::sqrt( Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b) + Tools::ggamma(i_c,j_c,i_d,j_d,k_c,k_d) ) );

            }
            else{

               R pq;

               if(i_a == i_b){

                  pq = input::gR(i_a);

                  if(i_c == i_d)
                     pq -= input::gR(i_c);
                  else
                     pq -= Tools::gP(i_c,j_c,i_d,j_d,k_c,k_d);

               }
               else{

                  pq = Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b);

                  if(i_c == i_d)
                     pq -= input::gR(i_c);
                  else
                     pq -= Tools::gP(i_c,j_c,i_d,j_d,k_c,k_d);

               }

               double pq2 = pq.ddot(pq);

               if(std::sqrt(pq2) < 1.0e-10){

                  double rho = Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b) * Tools::ggamma(i_c,j_c,i_d,j_d,k_c,k_d) /

                     ( Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b) + Tools::ggamma(i_c,j_c,i_d,j_d,k_c,k_d) );

                  //finally fill
                  (*this)(i_a,j_a,k_a,0,0,0,i_b,j_b,k_b,0,0,0,i_c,j_c,k_c,0,0,0,i_d,j_d,k_d,0,0,0) = 2.0 * std::sqrt(rho/M_PI)

                     * Tools::gS()(i_a,j_a,k_a,0,0,0,i_b,j_b,k_b,0,0,0) * Tools::gS()(i_c,j_c,k_c,0,0,0,i_d,j_d,k_d,0,0,0);

               }
               else{

                  double rho = Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b) * Tools::ggamma(i_c,j_c,i_d,j_d,k_c,k_d) /

                     ( Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b) + Tools::ggamma(i_c,j_c,i_d,j_d,k_c,k_d) );

                  //finally fill
                  (*this)(i_a,j_a,k_a,0,0,0,i_b,j_b,k_b,0,0,0,i_c,j_c,k_c,0,0,0,i_d,j_d,k_d,0,0,0) = erf( std::sqrt( rho * pq2 ) )

                     * Tools::gS()(i_a,j_a,k_a,0,0,0,i_b,j_b,k_b,0,0,0) * Tools::gS()(i_c,j_c,k_c,0,0,0,i_d,j_d,k_d,0,0,0) / std::sqrt(pq2);

               }

            }

         }
         else{//use libint

            if(l_a + l_b <= l_c + l_d){

               if(l_a < l_b){

                  if(l_c < l_d){

                     LibInt::prep_libint2(&inteval,l_b,input::gGaussInfo(i_b).galpha(j_b,k_b),input::gR_nc(i_b).gVector(), l_a, 

                           input::gGaussInfo(i_a).galpha(j_a,k_a),input::gR_nc(i_a).gVector(),l_d,input::gGaussInfo(i_d).galpha(j_d,k_d),

                           input::gR_nc(i_d).gVector(), l_c, input::gGaussInfo(i_c).galpha(j_c,k_c), input::gR_nc(i_c).gVector(), 0);


                     inteval.contrdepth = 1;

                     LIBINT2_PREFIXED_NAME (libint2_build_eri)[l_b][l_a][l_d][l_c](&inteval);

                     //extract the integrals!
                     int x_a,y_a,z_a;
                     int x_b,y_b,z_b;
                     int x_c,y_c,z_c;
                     int x_d,y_d,z_d;

                     int iter = 0;

                     //b
                     for(int I_b = 0;I_b <= l_b;++I_b){

                        x_b = l_b - I_b;

                        for(int J_b = 0;J_b <= I_b;++J_b){

                           y_b = I_b - J_b;
                           z_b = J_b;

                           //a
                           for(int I_a = 0;I_a <= l_a;++I_a){

                              x_a = l_a - I_a;

                              for(int J_a = 0;J_a <= I_a;++J_a){

                                 y_a = I_a - J_a;
                                 z_a = J_a;

                                 //d
                                 for(int I_d = 0;I_d <= l_d;++I_d){

                                    x_d = l_d - I_d;

                                    for(int J_d = 0;J_d <= I_d;++J_d){

                                       y_d = I_d - J_d;
                                       z_d = J_d;

                                       //c
                                       for(int I_c = 0;I_c <= l_c;++I_c){

                                          x_c = l_c - I_c;

                                          for(int J_c = 0;J_c <= I_c;++J_c){

                                             y_c  = I_c - J_c;
                                             z_c = J_c;

                                             (*this)(i_b,j_b,k_b,x_b,y_b,z_b,i_a,j_a,k_a,x_a,y_a,z_a,i_d,j_d,k_d,x_d,y_d,z_d,i_c,j_c,k_c,x_c,y_c,z_c)

                                                = inteval.targets[0][iter]; 

                                             ++iter;

                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }

                  }
                  else{

                     LibInt::prep_libint2(&inteval,l_b,input::gGaussInfo(i_b).galpha(j_b,k_b),input::gR_nc(i_b).gVector(), l_a, 

                           input::gGaussInfo(i_a).galpha(j_a,k_a),input::gR_nc(i_a).gVector(),l_c,input::gGaussInfo(i_c).galpha(j_c,k_c),

                           input::gR_nc(i_c).gVector(), l_d, input::gGaussInfo(i_d).galpha(j_d,k_d), input::gR_nc(i_d).gVector(), 0);

                     inteval.contrdepth = 1;

                     LIBINT2_PREFIXED_NAME (libint2_build_eri)[l_b][l_a][l_c][l_d](&inteval);

                     //extract the integrals!
                     int x_a,y_a,z_a;
                     int x_b,y_b,z_b;
                     int x_c,y_c,z_c;
                     int x_d,y_d,z_d;

                     int iter = 0;

                     //b
                     for(int I_b = 0;I_b <= l_b;++I_b){

                        x_b = l_b - I_b;

                        for(int J_b = 0;J_b <= I_b;++J_b){

                           y_b = I_b - J_b;
                           z_b = J_b;

                           //a
                           for(int I_a = 0;I_a <= l_a;++I_a){

                              x_a = l_a - I_a;

                              for(int J_a = 0;J_a <= I_a;++J_a){

                                 y_a = I_a - J_a;
                                 z_a = J_a;

                                 //c
                                 for(int I_c = 0;I_c <= l_c;++I_c){

                                    x_c = l_c - I_c;

                                    for(int J_c = 0;J_c <= I_c;++J_c){

                                       y_c = I_c - J_c;
                                       z_c = J_c;

                                       //d
                                       for(int I_d = 0;I_d <= l_d;++I_d){

                                          x_d = l_d - I_d;

                                          for(int J_d = 0;J_d <= I_d;++J_d){

                                             y_d  = I_d - J_d;
                                             z_d = J_d;

                                             (*this)(i_b,j_b,k_b,x_b,y_b,z_b,i_a,j_a,k_a,x_a,y_a,z_a,i_c,j_c,k_c,x_c,y_c,z_c,i_d,j_d,k_d,x_d,y_d,z_d)

                                                = inteval.targets[0][iter]; 

                                             ++iter;

                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }

                  }

               }
               else{

                  if(l_c < l_d){

                     LibInt::prep_libint2(&inteval,l_a,input::gGaussInfo(i_a).galpha(j_a,k_a),input::gR_nc(i_a).gVector(), l_b, 

                           input::gGaussInfo(i_b).galpha(j_b,k_b),input::gR_nc(i_b).gVector(),l_d,input::gGaussInfo(i_d).galpha(j_d,k_d),

                           input::gR_nc(i_d).gVector(), l_c, input::gGaussInfo(i_c).galpha(j_c,k_c), input::gR_nc(i_c).gVector(), 0);


                     inteval.contrdepth = 1;

                     LIBINT2_PREFIXED_NAME (libint2_build_eri)[l_a][l_b][l_d][l_c](&inteval);

                     //extract the integrals!
                     int x_a,y_a,z_a;
                     int x_b,y_b,z_b;
                     int x_c,y_c,z_c;
                     int x_d,y_d,z_d;

                     int iter = 0;

                     //a
                     for(int I_a = 0;I_a <= l_a;++I_a){

                        x_a = l_a - I_a;

                        for(int J_a = 0;J_a <= I_a;++J_a){

                           y_a = I_a - J_a;
                           z_a = J_a;

                           //b
                           for(int I_b = 0;I_b <= l_b;++I_b){

                              x_b = l_b - I_b;

                              for(int J_b = 0;J_b <= I_b;++J_b){

                                 y_b = I_b - J_b;
                                 z_b = J_b;

                                 //d
                                 for(int I_d = 0;I_d <= l_d;++I_d){

                                    x_d = l_d - I_d;

                                    for(int J_d = 0;J_d <= I_d;++J_d){

                                       y_d = I_d - J_d;
                                       z_d = J_d;

                                       //c
                                       for(int I_c = 0;I_c <= l_c;++I_c){

                                          x_c = l_c - I_c;

                                          for(int J_c = 0;J_c <= I_c;++J_c){

                                             y_c  = I_c - J_c;
                                             z_c = J_c;

                                             (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b,i_d,j_d,k_d,x_d,y_d,z_d,i_c,j_c,k_c,x_c,y_c,z_c)

                                                = inteval.targets[0][iter]; 

                                             ++iter;

                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }

                  }
                  else{

                     LibInt::prep_libint2(&inteval,l_a,input::gGaussInfo(i_a).galpha(j_a,k_a),input::gR_nc(i_a).gVector(), l_b, 

                           input::gGaussInfo(i_b).galpha(j_b,k_b),input::gR_nc(i_b).gVector(),l_c,input::gGaussInfo(i_c).galpha(j_c,k_c),

                           input::gR_nc(i_c).gVector(), l_d, input::gGaussInfo(i_d).galpha(j_d,k_d), input::gR_nc(i_d).gVector(), 0);


                     inteval.contrdepth = 1;

                     LIBINT2_PREFIXED_NAME (libint2_build_eri)[l_a][l_b][l_c][l_d](&inteval);

                     //extract the integrals!
                     int x_a,y_a,z_a;
                     int x_b,y_b,z_b;
                     int x_c,y_c,z_c;
                     int x_d,y_d,z_d;

                     int iter = 0;

                     //a
                     for(int I_a = 0;I_a <= l_a;++I_a){

                        x_a = l_a - I_a;

                        for(int J_a = 0;J_a <= I_a;++J_a){

                           y_a = I_a - J_a;
                           z_a = J_a;

                           //b
                           for(int I_b = 0;I_b <= l_b;++I_b){

                              x_b = l_b - I_b;

                              for(int J_b = 0;J_b <= I_b;++J_b){

                                 y_b = I_b - J_b;
                                 z_b = J_b;

                                 //c
                                 for(int I_c = 0;I_c <= l_c;++I_c){

                                    x_c = l_c - I_c;

                                    for(int J_c = 0;J_c <= I_c;++J_c){

                                       y_c = I_c - J_c;
                                       z_c = J_c;

                                       //d
                                       for(int I_d = 0;I_d <= l_d;++I_d){

                                          x_d = l_d - I_d;

                                          for(int J_d = 0;J_d <= I_d;++J_d){

                                             y_d = I_d - J_d;
                                             z_d = J_d;

                                             (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b,i_c,j_c,k_c,x_c,y_c,z_c,i_d,j_d,k_d,x_d,y_d,z_d)

                                                = inteval.targets[0][iter]; 

                                             ++iter;

                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }

                  }

               }

            }//tot hier
            else{//la + lb > lc + ld

               if(l_c < l_d){

                  if(l_a < l_b){

                     LibInt::prep_libint2(&inteval,l_d,input::gGaussInfo(i_d).galpha(j_d,k_d),input::gR_nc(i_d).gVector(), l_c, 

                           input::gGaussInfo(i_c).galpha(j_c,k_c),input::gR_nc(i_c).gVector(),l_b,input::gGaussInfo(i_b).galpha(j_b,k_b),

                           input::gR_nc(i_b).gVector(), l_a, input::gGaussInfo(i_a).galpha(j_a,k_a), input::gR_nc(i_a).gVector(), 0);


                     inteval.contrdepth = 1;

                     LIBINT2_PREFIXED_NAME (libint2_build_eri)[l_d][l_c][l_b][l_a](&inteval);

                     //extract the integrals!
                     int x_a,y_a,z_a;
                     int x_b,y_b,z_b;
                     int x_c,y_c,z_c;
                     int x_d,y_d,z_d;

                     int iter = 0;

                     //d
                     for(int I_d = 0;I_d <= l_d;++I_d){

                        x_d = l_d - I_d;

                        for(int J_d = 0;J_d <= I_d;++J_d){

                           y_d = I_d - J_d;
                           z_d = J_d;

                           //c
                           for(int I_c = 0;I_c <= l_c;++I_c){

                              x_c = l_c - I_c;

                              for(int J_c = 0;J_c <= I_c;++J_c){

                                 y_c = I_c - J_c;
                                 z_c = J_c;

                                 //b
                                 for(int I_b = 0;I_b <= l_b;++I_b){

                                    x_b = l_b - I_b;

                                    for(int J_b = 0;J_b <= I_b;++J_b){

                                       y_b = I_b - J_b;
                                       z_b = J_b;

                                       //a
                                       for(int I_a = 0;I_a <= l_a;++I_a){

                                          x_a = l_a - I_a;

                                          for(int J_a = 0;J_a <= I_a;++J_a){

                                             y_a  = I_a - J_a;
                                             z_a = J_a;

                                             (*this)(i_d,j_d,k_d,x_d,y_d,z_d,i_c,j_c,k_c,x_c,y_c,z_c,i_b,j_b,k_b,x_b,y_b,z_b,i_a,j_a,k_a,x_a,y_a,z_a)

                                                = inteval.targets[0][iter]; 

                                             ++iter;

                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }

                  }
                  else{

                     LibInt::prep_libint2(&inteval,l_d,input::gGaussInfo(i_d).galpha(j_d,k_d),input::gR_nc(i_d).gVector(), l_c, 

                           input::gGaussInfo(i_c).galpha(j_c,k_c),input::gR_nc(i_c).gVector(),l_a,input::gGaussInfo(i_a).galpha(j_a,k_a),

                           input::gR_nc(i_a).gVector(), l_b, input::gGaussInfo(i_b).galpha(j_b,k_b), input::gR_nc(i_b).gVector(), 0);


                     inteval.contrdepth = 1;

                     LIBINT2_PREFIXED_NAME (libint2_build_eri)[l_d][l_c][l_a][l_b](&inteval);

                     //extract the integrals!
                     int x_a,y_a,z_a;
                     int x_b,y_b,z_b;
                     int x_c,y_c,z_c;
                     int x_d,y_d,z_d;

                     int iter = 0;

                     //d
                     for(int I_d = 0;I_d <= l_d;++I_d){

                        x_d = l_d - I_d;

                        for(int J_d = 0;J_d <= I_d;++J_d){

                           y_d = I_d - J_d;
                           z_d = J_d;

                           //c
                           for(int I_c = 0;I_c <= l_c;++I_c){

                              x_c = l_c - I_c;

                              for(int J_c = 0;J_c <= I_c;++J_c){

                                 y_c = I_c - J_c;
                                 z_c = J_c;

                                 //a
                                 for(int I_a = 0;I_a <= l_a;++I_a){

                                    x_a = l_a - I_a;

                                    for(int J_a = 0;J_a <= I_a;++J_a){

                                       y_a = I_a - J_a;
                                       z_a = J_a;

                                       //b
                                       for(int I_b = 0;I_b <= l_b;++I_b){

                                          x_b = l_b - I_b;

                                          for(int J_b = 0;J_b <= I_b;++J_b){

                                             y_b  = I_b - J_b;
                                             z_b = J_b;

                                             (*this)(i_d,j_d,k_d,x_d,y_d,z_d,i_c,j_c,k_c,x_c,y_c,z_c,i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b)

                                                = inteval.targets[0][iter]; 

                                             ++iter;

                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }

                  }

               }
               else{

                  if(l_a < l_b){

                     LibInt::prep_libint2(&inteval,l_c,input::gGaussInfo(i_c).galpha(j_c,k_c),input::gR_nc(i_c).gVector(), l_d, 

                           input::gGaussInfo(i_d).galpha(j_d,k_d),input::gR_nc(i_d).gVector(),l_b,input::gGaussInfo(i_b).galpha(j_b,k_b),

                           input::gR_nc(i_b).gVector(), l_a, input::gGaussInfo(i_a).galpha(j_a,k_a), input::gR_nc(i_a).gVector(), 0);


                     inteval.contrdepth = 1;

                     LIBINT2_PREFIXED_NAME (libint2_build_eri)[l_c][l_d][l_b][l_a](&inteval);

                     //extract the integrals!
                     int x_a,y_a,z_a;
                     int x_b,y_b,z_b;
                     int x_c,y_c,z_c;
                     int x_d,y_d,z_d;

                     int iter = 0;

                     //c
                     for(int I_c = 0;I_c <= l_c;++I_c){

                        x_c = l_c - I_c;

                        for(int J_c = 0;J_c <= I_c;++J_c){

                           y_c = I_c - J_c;
                           z_c = J_c;

                           //d
                           for(int I_d = 0;I_d <= l_d;++I_d){

                              x_d = l_d - I_d;

                              for(int J_d = 0;J_d <= I_d;++J_d){

                                 y_d = I_d - J_d;
                                 z_d = J_d;

                                 //b
                                 for(int I_b = 0;I_b <= l_b;++I_b){

                                    x_b = l_b - I_b;

                                    for(int J_b = 0;J_b <= I_b;++J_b){

                                       y_b = I_b - J_b;
                                       z_b = J_b;

                                       //a
                                       for(int I_a = 0;I_a <= l_a;++I_a){

                                          x_a = l_a - I_a;

                                          for(int J_a = 0;J_a <= I_a;++J_a){

                                             y_a  = I_a - J_a;
                                             z_a = J_a;

                                             (*this)(i_c,j_c,k_c,x_c,y_c,z_c,i_d,j_d,k_d,x_d,y_d,z_d,i_b,j_b,k_b,x_b,y_b,z_b,i_a,j_a,k_a,x_a,y_a,z_a)

                                                = inteval.targets[0][iter]; 

                                             ++iter;

                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }

                  }
                  else{

                     LibInt::prep_libint2(&inteval,l_c,input::gGaussInfo(i_c).galpha(j_c,k_c),input::gR_nc(i_c).gVector(), l_d, 

                           input::gGaussInfo(i_d).galpha(j_d,k_d),input::gR_nc(i_d).gVector(),l_a,input::gGaussInfo(i_a).galpha(j_a,k_a),

                           input::gR_nc(i_a).gVector(), l_b, input::gGaussInfo(i_b).galpha(j_b,k_b), input::gR_nc(i_b).gVector(), 0);


                     inteval.contrdepth = 1;

                     LIBINT2_PREFIXED_NAME (libint2_build_eri)[l_c][l_d][l_a][l_b](&inteval);

                     //extract the integrals!
                     int x_a,y_a,z_a;
                     int x_b,y_b,z_b;
                     int x_c,y_c,z_c;
                     int x_d,y_d,z_d;

                     int iter = 0;

                     //c
                     for(int I_c = 0;I_c <= l_c;++I_c){

                        x_c = l_c - I_c;

                        for(int J_c = 0;J_c <= I_c;++J_c){

                           y_c = I_c - J_c;
                           z_c = J_c;

                           //d
                           for(int I_d = 0;I_d <= l_d;++I_d){

                              x_d = l_d - I_d;

                              for(int J_d = 0;J_d <= I_d;++J_d){

                                 y_d = I_d - J_d;
                                 z_d = J_d;

                                 //a
                                 for(int I_a = 0;I_a <= l_a;++I_a){

                                    x_a = l_a - I_a;

                                    for(int J_a = 0;J_a <= I_a;++J_a){

                                       y_a = I_a - J_a;
                                       z_a = J_a;

                                       //b
                                       for(int I_b = 0;I_b <= l_b;++I_b){

                                          x_b = l_b - I_b;

                                          for(int J_b = 0;J_b <= I_b;++J_b){

                                             y_b  = I_b - J_b;
                                             z_b = J_b;

                                             (*this)(i_c,j_c,k_c,x_c,y_c,z_c,i_d,j_d,k_d,x_d,y_d,z_d,i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b)

                                                = inteval.targets[0][iter]; 

                                             ++iter;

                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }

                  }

               }

            }

         }

      }

   }

   LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval);

   this->symmetrize();

}

/**
 * access to the elements using quantumnumbers instead of matrixindex
 */
double CI_TPPM::operator()(int i_a,int j_a,int k_a,int x_a,int y_a,int z_a,int i_b,int j_b,int k_b,int x_b,int y_b,int z_b,

      int i_c,int j_c,int k_c,int x_c,int y_c,int z_c,int i_d,int j_d,int k_d,int x_d,int y_d,int z_d ) const {

   int omi_a = CI_SPPM::gijkxyz2omi(i_a,j_a,k_a,x_a,y_a,z_a);
   int omi_b = CI_SPPM::gijkxyz2omi(i_b,j_b,k_b,x_b,y_b,z_b);
   int omi_c = CI_SPPM::gijkxyz2omi(i_c,j_c,k_c,x_c,y_c,z_c);
   int omi_d = CI_SPPM::gijkxyz2omi(i_d,j_d,k_d,x_d,y_d,z_d);

   int I,J;

   if(omi_a > omi_b)
      I = omi2tomi[omi_b][omi_a];
   else
      I = omi2tomi[omi_a][omi_b];

   if(omi_c > omi_d)
      J = omi2tomi[omi_d][omi_c];
   else
      J = omi2tomi[omi_c][omi_d];

   return (*this)(I,J);

}

/**
 * access to the elements using quantumnumbers instead of matrixindex
 */
double &CI_TPPM::operator()(int i_a,int j_a,int k_a,int x_a,int y_a,int z_a,int i_b,int j_b,int k_b,int x_b,int y_b,int z_b,

      int i_c,int j_c,int k_c,int x_c,int y_c,int z_c,int i_d,int j_d,int k_d,int x_d,int y_d,int z_d ) {

   int omi_a = CI_SPPM::gijkxyz2omi(i_a,j_a,k_a,x_a,y_a,z_a);
   int omi_b = CI_SPPM::gijkxyz2omi(i_b,j_b,k_b,x_b,y_b,z_b);
   int omi_c = CI_SPPM::gijkxyz2omi(i_c,j_c,k_c,x_c,y_c,z_c);
   int omi_d = CI_SPPM::gijkxyz2omi(i_d,j_d,k_d,x_d,y_d,z_d);

   int I,J;

   if(omi_a > omi_b)
      I = omi2tomi[omi_b][omi_a];
   else
      I = omi2tomi[omi_a][omi_b];

   if(omi_c > omi_d)
      J = omi2tomi[omi_d][omi_c];
   else
      J = omi2tomi[omi_c][omi_d];

   if(I <= J)
      return (*this)(I,J);
   else
      return (*this)(J,I);

}
