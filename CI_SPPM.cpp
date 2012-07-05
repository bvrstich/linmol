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

vector< vector<int> > CI_SPPM::omi2ijkxyz;
int ******CI_SPPM::ijkxyz2omi;

int **CI_SPPM::ammax;

vector<int> CI_SPPM::amlist;

/** 
 * static function that reads in the input data and makes the corresponding lists and dimensions
 */
void CI_SPPM::init(){

   vector<int> v(6);

   ammax = new int * [input::gN_Z()];

   ijkxyz2omi = new int ***** [input::gN_Z()];

   int iter = 0;

   for(int i = 0;i < input::gN_Z();++i){

      v[0] = i;

      ammax[i] = new int [input::gGaussInfo(i).gNtypes()];

      ijkxyz2omi[i] = new int **** [input::gGaussInfo(i).gNtypes()];

      for(int j = 0;j < input::gGaussInfo(i).gNtypes();++j){

         v[1] = j;

         if(input::gGaussInfo(i).gtype(j) == 'S')
            ammax[i][j] = 0;
         else if(input::gGaussInfo(i).gtype(j) == 'P')
            ammax[i][j] = 1;
         else if(input::gGaussInfo(i).gtype(j) == 'D')
            ammax[i][j] = 2;
         else if(input::gGaussInfo(i).gtype(j) == 'F')
            ammax[i][j] = 3;
         else if(input::gGaussInfo(i).gtype(j) == 'G')
            ammax[i][j] = 4;
         else if(input::gGaussInfo(i).gtype(j) == 'H')
            ammax[i][j] = 5;
         else
            cout << "BASISSET TOO LARGE" << endl;

         ijkxyz2omi[i][j] = new int *** [input::gGaussInfo(i).gNcontr(j)];

         for(int k = 0;k < input::gGaussInfo(i).gNcontr(j);++k){

            //now allocate the rest of ijkxyz2omi
            ijkxyz2omi[i][j][k] = new int ** [ammax[i][j] + 1];

            for(int x = 0;x <= ammax[i][j];++x){

               ijkxyz2omi[i][j][k][x] = new int * [ammax[i][j] + 1];

               for(int y = 0;y <= ammax[i][j];++y)
                  ijkxyz2omi[i][j][k][x][y] = new int [ammax[i][j] + 1];

            }

            v[2] = k;

            //now loop over the different possibilities for x,y and z
            for(int am = 0;am <= ammax[i][j];++am){

               for(int x = am;x >= 0;x--)
                  for(int y = am - x;y >= 0;y--){

                     v[3] = x;
                     v[4] = y;
                     v[5] = am - x - y;

                     omi2ijkxyz.push_back(v);

                     amlist.push_back(am);

                     ijkxyz2omi[i][j][k][x][y][am - x - y] = iter;

                     ++iter;

                  }

            }

         }

      }

   }

}

/** 
 * function that deallocates the static members
 */
void CI_SPPM::clear(){

   for(int i = 0;i < input::gN_Z();++i){

      for(int j = 0;j < input::gGaussInfo(i).gNtypes();++j){

         for(int k = 0;k < input::gGaussInfo(i).gNcontr(j);++k){

            for(int x = 0;x <= ammax[i][j];++x){

               for(int y = 0;y <= ammax[i][j];++y)
                  delete [] ijkxyz2omi[i][j][k][x][y];

               delete [] ijkxyz2omi[i][j][k][x];

            }

            delete [] ijkxyz2omi[i][j][k];

         }

         delete [] ijkxyz2omi[i][j];

      }

      delete [] ijkxyz2omi[i];

   }

   delete [] ijkxyz2omi;

  for(int i = 0;i < input::gN_Z();++i)
     delete [] ammax[i];

   delete [] ammax;

}

/** 
 * Standard constructor: allocates the matrices and fills them with the correct elements
 */
CI_SPPM::CI_SPPM() : Matrix(omi2ijkxyz.size()) { }

/** 
 * copy constructor
 * @param ci_c CI_SPPM object to be copied in the newly constructed object
 */
CI_SPPM::CI_SPPM(const CI_SPPM &ci_c) : Matrix(ci_c){ }

/**
 * standard destructor
 */
CI_SPPM::~CI_SPPM(){ }

/**
 * access to the elements using quantumnumbers instead of matrixindex
 */
double CI_SPPM::operator()(int i_a,int j_a,int k_a,int x_a,int y_a,int z_a,int i_b,int j_b,int k_b,int x_b,int y_b,int z_b) const {

   return (*this)(ijkxyz2omi[i_a][j_a][k_a][x_a][y_a][z_a],ijkxyz2omi[i_b][j_b][k_b][x_b][y_b][z_b]);

}

/**
 * static function which allows access to the lists
 */
int CI_SPPM::gomi2ijkxyz(int omi,int option){

   return omi2ijkxyz[omi][option];

}

/**
 * static function which allows access to the lists
 */
int CI_SPPM::gijkxyz2omi(int i,int j,int k,int x,int y,int z){

   return ijkxyz2omi[i][j][k][x][y][z];

}

/**
 * construct the overlapmatrix of all the different primitives.
 */
void CI_SPPM::constr_overlap(){

   int i_a,j_a,k_a,x_a,y_a,z_a;
   int i_b,j_b,k_b,x_b,y_b,z_b;

   for(int omi_a = 0;omi_a < gn();++omi_a){

      i_a = omi2ijkxyz[omi_a][0];
      j_a = omi2ijkxyz[omi_a][1];
      k_a = omi2ijkxyz[omi_a][2];
      x_a = omi2ijkxyz[omi_a][3];
      y_a = omi2ijkxyz[omi_a][4];
      z_a = omi2ijkxyz[omi_a][5];

      for(int omi_b = omi_a;omi_b < gn();++omi_b){

         i_b = omi2ijkxyz[omi_b][0];
         j_b = omi2ijkxyz[omi_b][1];
         k_b = omi2ijkxyz[omi_b][2];
         x_b = omi2ijkxyz[omi_b][3];
         y_b = omi2ijkxyz[omi_b][4];
         z_b = omi2ijkxyz[omi_b][5];

         if(i_a == i_b){//on core

            if( (x_a + x_b)%2 != 0 || (y_a + y_b)%2 != 0 || (z_a + z_b)%2 != 0 )
               (*this)(omi_a,omi_b) = (*this)(omi_b,omi_a) = 0.0;
            else if(amlist[omi_a] == 0 && amlist[omi_b] == 0){//(0_a|0_b)

               (*this)(omi_a,omi_b) = pow(M_PI/Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b),1.5);

               (*this)(omi_b,omi_a) = (*this)(omi_a,omi_b);

            }
            else{//higher angular momenta

               //which has largest angular momentum, a or b?
               if(amlist[omi_b] >= amlist[omi_a]){

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

                        (*this)(omi_a,omi_b) = (x_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 2,y_b,z_b);
                     }
                     else
                        (*this)(omi_a,omi_b) = 0.0;

                     if(x_a > 0){

                        (*this)(omi_a,omi_b) += x_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b);

                     }

                  }
                  else if(flag == 1){//y_b is largest

                     if(y_b > 1){

                        (*this)(omi_a,omi_b) = (y_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 2,z_b);
                     }
                     else
                        (*this)(omi_a,omi_b) = 0.0;

                     if(y_a > 0){

                        (*this)(omi_a,omi_b) += y_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b);

                     }

                  }
                  else{//z_b is largest

                     if(z_b > 1){

                        (*this)(omi_a,omi_b) = (z_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 2);
                     }
                     else
                        (*this)(omi_a,omi_b) = 0.0;

                     if(z_a > 0){

                        (*this)(omi_a,omi_b) += z_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1);

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

                        (*this)(omi_a,omi_b) = (x_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a - 2,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);
                     }
                     else
                        (*this)(omi_a,omi_b) = 0.0;

                     if(x_b > 0){

                        (*this)(omi_a,omi_b) += x_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b);

                     }

                  }
                  else if(flag == 1){//y_a is largest

                     if(y_a > 1){

                        (*this)(omi_a,omi_b) = (y_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a - 2,z_a,i_b,j_b,k_b,x_b,y_b,z_b);
                     }
                     else
                        (*this)(omi_a,omi_b) = 0.0;

                     if(y_b > 0){

                        (*this)(omi_a,omi_b) += y_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b);

                     }

                  }
                  else{//z_a is largest

                     if(z_a > 1){

                        (*this)(omi_a,omi_b) = (z_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 2,i_b,j_b,k_b,x_b,y_b,z_b);
                     }
                     else
                        (*this)(omi_a,omi_b) = 0.0;

                     if(z_b > 0){

                        (*this)(omi_a,omi_b) += z_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1);

                     }

                  }

               }

               (*this)(omi_b,omi_a) = (*this)(omi_a,omi_b);

            }

         }
         else{//off core

            if(amlist[omi_a] == 0 && amlist[omi_b] == 0){//(0_A|0_B)

               (*this)(omi_a,omi_b) = pow(M_PI/Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b),1.5) * Tools::gprefac_overlap(i_a,j_a,i_b,j_b,k_a,k_b);

               (*this)(omi_b,omi_a) = (*this)(omi_a,omi_b);

            }
            else{

               //which has largest angular momentum, a or b?
               if(amlist[omi_b] >= amlist[omi_a]){

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

                     (*this)(omi_a,omi_b) = (Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[0] - input::gR(i_b)[0]) 

                        * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b);

                     if(x_b > 1){

                        (*this)(omi_a,omi_b) += (x_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 2,y_b,z_b);
                     }

                     if(x_a > 0){

                        (*this)(omi_a,omi_b) += x_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b);

                     }

                  }
                  else if(flag == 1){//y_b is largest

                     (*this)(omi_a,omi_b) = (Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[1] - input::gR(i_b)[1]) 

                        * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b);

                     if(y_b > 1){

                        (*this)(omi_a,omi_b) += (y_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 2,z_b);

                     }

                     if(y_a > 0){

                        (*this)(omi_a,omi_b) += y_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b);

                     }

                  }
                  else{//z_b is largest

                     (*this)(omi_a,omi_b) = (Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[2] - input::gR(i_b)[2]) 

                        * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 1);

                     if(z_b > 1){

                        (*this)(omi_a,omi_b) += (z_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 2);
                     }

                     if(z_a > 0){

                        (*this)(omi_a,omi_b) += z_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1);

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

                     (*this)(omi_a,omi_b) = (Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[0] - input::gR(i_a)[0]) 

                        * (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     if(x_a > 1){

                        (*this)(omi_a,omi_b) += (x_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a - 2,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);
                     }

                     if(x_b > 0){

                        (*this)(omi_a,omi_b) += x_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b);

                     }

                  }
                  else if(flag == 1){//y_a is largest

                     (*this)(omi_a,omi_b) = (Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[1] - input::gR(i_a)[1]) 

                        * (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     if(y_a > 1){

                        (*this)(omi_a,omi_b) += (y_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a - 2,z_a,i_b,j_b,k_b,x_b,y_b,z_b);
                     }

                     if(y_b > 0){

                        (*this)(omi_a,omi_b) += y_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b);

                     }

                  }
                  else{//z_a is largest

                     (*this)(omi_a,omi_b) = (Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[2] - input::gR(i_a)[2]) 

                        * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b);

                     if(z_a > 1){

                        (*this)(omi_a,omi_b) += (z_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 2,i_b,j_b,k_b,x_b,y_b,z_b);
                     }

                     if(z_b > 0){

                        (*this)(omi_a,omi_b) += z_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1);

                     }

                  }

               }

               (*this)(omi_b,omi_a) = (*this)(omi_a,omi_b);

            }

         }

      }
   }

}

/**
 * construct the kinetic energy matrix of all the different primitives.
 */
void CI_SPPM::constr_kinetic(){

   int i_a,j_a,k_a,x_a,y_a,z_a;
   int i_b,j_b,k_b,x_b,y_b,z_b;

   for(int omi_a = 0;omi_a < gn();++omi_a){

      i_a = omi2ijkxyz[omi_a][0];
      j_a = omi2ijkxyz[omi_a][1];
      k_a = omi2ijkxyz[omi_a][2];
      x_a = omi2ijkxyz[omi_a][3];
      y_a = omi2ijkxyz[omi_a][4];
      z_a = omi2ijkxyz[omi_a][5];

      for(int omi_b = omi_a;omi_b < gn();++omi_b){

         i_b = omi2ijkxyz[omi_b][0];
         j_b = omi2ijkxyz[omi_b][1];
         k_b = omi2ijkxyz[omi_b][2];
         x_b = omi2ijkxyz[omi_b][3];
         y_b = omi2ijkxyz[omi_b][4];
         z_b = omi2ijkxyz[omi_b][5];

         if(i_a == i_b){//on core

            if( (x_a + x_b)%2 != 0 || (y_a + y_b)%2 != 0 || (z_a + z_b)%2 != 0 )
               (*this)(omi_a,omi_b) = (*this)(omi_b,omi_a) = 0.0;
            else if(amlist[omi_a] == 0 && amlist[omi_b] == 0){//(0_a|T|0_b)

               double xi = input::gGaussInfo(i_a).galpha(j_a,k_a) * input::gGaussInfo(i_b).galpha(j_b,k_b)/Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b);

               (*this)(omi_a,omi_b) = 3.0 * xi * Tools::gS()(i_a,j_a,k_a,0,0,0,i_b,j_b,k_b,0,0,0);

               (*this)(omi_b,omi_a) = (*this)(omi_a,omi_b);

            }
            else{//higher am

               //which has largest angular momentum, a or b?
               if(amlist[omi_b] >= amlist[omi_a]){

                  double xi = input::gGaussInfo(i_a).galpha(j_a,k_a) * input::gGaussInfo(i_b).galpha(j_b,k_b)/Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b);

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

                     (*this)(omi_a,omi_b) = 2.0 * xi * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     if(x_b > 1){

                        (*this)(omi_a,omi_b) += (x_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 2,y_b,z_b);

                        (*this)(omi_a,omi_b) -= xi * (x_b - 1)/input::gGaussInfo(i_b).galpha(j_b,k_b)

                           * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 2,y_b,z_b);
                     }

                     if(x_a > 0){

                        (*this)(omi_a,omi_b) += x_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b);

                     }

                  }
                  else if(flag == 1){//y_b is largest

                     (*this)(omi_a,omi_b) = 2.0 * xi * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     if(y_b > 1){

                        (*this)(omi_a,omi_b) += (y_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 2,z_b);

                        (*this)(omi_a,omi_b) -= xi * (y_b - 1)/input::gGaussInfo(i_b).galpha(j_b,k_b)

                           * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 2,z_b);
                     }

                     if(y_a > 0){

                        (*this)(omi_a,omi_b) += y_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b);

                     }

                  }
                  else{//z_b is largest

                     (*this)(omi_a,omi_b) = 2.0 * xi * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     if(z_b > 1){

                        (*this)(omi_a,omi_b) += (z_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 2);

                        (*this)(omi_a,omi_b) -= xi * (z_b - 1)/input::gGaussInfo(i_b).galpha(j_b,k_b)

                           * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 2);
                     }

                     if(z_a > 0){

                        (*this)(omi_a,omi_b) += z_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1);

                     }

                  }

               }
               else{//am of a is larger than that of b

                  double xi = input::gGaussInfo(i_a).galpha(j_a,k_a) * input::gGaussInfo(i_b).galpha(j_b,k_b)/Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b);

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

                     (*this)(omi_a,omi_b) = 2.0 * xi * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     if(x_a > 1){

                        (*this)(omi_a,omi_b) += (x_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a - 2,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                        (*this)(omi_a,omi_b) -= xi * (x_a - 1)/input::gGaussInfo(i_a).galpha(j_a,k_a)

                           * Tools::gS()(i_a,j_a,k_a,x_a - 2,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);
                     }

                     if(x_b > 0){

                        (*this)(omi_a,omi_b) += x_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b);

                     }

                  }
                  else if(flag == 1){//y_a is largest

                     (*this)(omi_a,omi_b) = 2.0 * xi * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     if(y_a > 1){

                        (*this)(omi_a,omi_b) += (y_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a - 2,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                        (*this)(omi_a,omi_b) -= xi * (y_a - 1)/input::gGaussInfo(i_a).galpha(j_a,k_a)

                           * Tools::gS()(i_a,j_a,k_a,x_a,y_a - 2,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     }

                     if(y_b > 0){

                        (*this)(omi_a,omi_b) += y_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b);

                     }

                  }
                  else{//z_a is largest

                     (*this)(omi_a,omi_b) = 2.0 * xi * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     if(z_a > 1){

                        (*this)(omi_a,omi_b) += (z_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 2,i_b,j_b,k_b,x_b,y_b,z_b);

                        (*this)(omi_a,omi_b) -= xi * (z_a - 1)/input::gGaussInfo(i_a).galpha(j_a,k_a)

                           * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a - 2,i_b,j_b,k_b,x_b,y_b,z_b);

                     }

                     if(z_b > 0){

                        (*this)(omi_a,omi_b) += z_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1);

                     }

                  }

               }

               (*this)(omi_b,omi_a) = (*this)(omi_a,omi_b);

            }

         }
         else{//off-core

            if(amlist[omi_a] == 0 && amlist[omi_b] == 0){//(0_a|T|0_b)

               double xi = input::gGaussInfo(i_a).galpha(j_a,k_a) * input::gGaussInfo(i_b).galpha(j_b,k_b)/Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b);

               //distance between a and b
               R ab(input::gR(i_a));
               ab -= input::gR(i_b);

               (*this)(omi_a,omi_b) = xi * (3.0  - 2 * xi * ab.ddot(ab)) * Tools::gS()(i_a,j_a,k_a,0,0,0,i_b,j_b,k_b,0,0,0);

               (*this)(omi_b,omi_a) = (*this)(omi_a,omi_b);

            }
            else{

               //which has largest angular momentum, a or b?
               if(amlist[omi_b] >= amlist[omi_a]){

                  double xi = input::gGaussInfo(i_a).galpha(j_a,k_a) * input::gGaussInfo(i_b).galpha(j_b,k_b)/Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b);

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

                     (*this)(omi_a,omi_b) = 2.0 * xi * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     (*this)(omi_a,omi_b) += (Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[0] - input::gR(i_b)[0]) 

                        * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b);

                     if(x_b > 1){

                        (*this)(omi_a,omi_b) += (x_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 2,y_b,z_b);

                        (*this)(omi_a,omi_b) -= xi * (x_b - 1)/input::gGaussInfo(i_b).galpha(j_b,k_b)

                           * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b - 2,y_b,z_b);
                     }

                     if(x_a > 0){

                        (*this)(omi_a,omi_b) += x_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b);

                     }

                  }
                  else if(flag == 1){//y_b is largest

                     (*this)(omi_a,omi_b) = 2.0 * xi * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     (*this)(omi_a,omi_b) += (Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[1] - input::gR(i_b)[1]) 

                        * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b);

                     if(y_b > 1){

                        (*this)(omi_a,omi_b) += (y_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 2,z_b);

                        (*this)(omi_a,omi_b) -= xi * (y_b - 1)/input::gGaussInfo(i_b).galpha(j_b,k_b)

                           * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b - 2,z_b);
                     }

                     if(y_a > 0){

                        (*this)(omi_a,omi_b) += y_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b);

                     }

                  }
                  else{//z_b is largest

                     (*this)(omi_a,omi_b) = 2.0 * xi * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     (*this)(omi_a,omi_b) += (Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[2] - input::gR(i_b)[2]) 

                        * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 1);

                     if(z_b > 1){

                        (*this)(omi_a,omi_b) += (z_b - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 2);

                        (*this)(omi_a,omi_b) -= xi * (z_b - 1)/input::gGaussInfo(i_b).galpha(j_b,k_b)

                           * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b - 2);
                     }

                     if(z_a > 0){

                        (*this)(omi_a,omi_b) += z_a /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1);

                     }

                  }

               }
               else{//am of a is larger than that of b

                  double xi = input::gGaussInfo(i_a).galpha(j_a,k_a) * input::gGaussInfo(i_b).galpha(j_b,k_b)/Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b);

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

                     (*this)(omi_a,omi_b) = 2.0 * xi * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     (*this)(omi_a,omi_b) += (Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[0] - input::gR(i_a)[0]) 

                        * (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     if(x_a > 1){

                        (*this)(omi_a,omi_b) += (x_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a - 2,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                        (*this)(omi_a,omi_b) -= xi * (x_a - 1)/input::gGaussInfo(i_a).galpha(j_a,k_a)

                           * Tools::gS()(i_a,j_a,k_a,x_a - 2,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);
                     }

                     if(x_b > 0){

                        (*this)(omi_a,omi_b) += x_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a - 1,y_a,z_a,i_b,j_b,k_b,x_b - 1,y_b,z_b);

                     }

                  }
                  else if(flag == 1){//y_a is largest

                     (*this)(omi_a,omi_b) = 2.0 * xi * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     (*this)(omi_a,omi_b) += (Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[1] - input::gR(i_a)[1]) 

                        * (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     if(y_a > 1){

                        (*this)(omi_a,omi_b) += (y_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a - 2,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                        (*this)(omi_a,omi_b) -= xi * (y_a - 1)/input::gGaussInfo(i_a).galpha(j_a,k_a)

                           * Tools::gS()(i_a,j_a,k_a,x_a,y_a - 2,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     }

                     if(y_b > 0){

                        (*this)(omi_a,omi_b) += y_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a - 1,z_a,i_b,j_b,k_b,x_b,y_b - 1,z_b);

                     }

                  }
                  else{//z_a is largest

                     (*this)(omi_a,omi_b) = 2.0 * xi * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

                     (*this)(omi_a,omi_b) += (Tools::gP(i_a,j_a,i_b,j_b,k_a,k_b)[2] - input::gR(i_a)[2]) 

                        * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b);

                     if(z_a > 1){

                        (*this)(omi_a,omi_b) += (z_a - 1)/(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 2,i_b,j_b,k_b,x_b,y_b,z_b);

                        (*this)(omi_a,omi_b) -= xi * (z_a - 1)/input::gGaussInfo(i_a).galpha(j_a,k_a)

                           * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a - 2,i_b,j_b,k_b,x_b,y_b,z_b);

                     }

                     if(z_b > 0){

                        (*this)(omi_a,omi_b) += z_b /(2.0*Tools::ggamma(i_a,j_a,i_b,j_b,k_a,k_b))

                           * (*this)(i_a,j_a,k_a,x_a,y_a,z_a - 1,i_b,j_b,k_b,x_b,y_b,z_b - 1);

                     }

                  }

               }

               (*this)(omi_b,omi_a) = (*this)(omi_a,omi_b);

            }

         }

      }
   }

   this->symmetrize();

}

/**
 * @return the dimension of the matrix without the need for an actual object
 */
int CI_SPPM::gdim(){

   return omi2ijkxyz.size();

}

/**
 * access to some static lists
 * @param a the CI_SPPM matrix index
 * @return the angular momentum, i.e. x_a + y_a + z_a of this particular index
 */
int CI_SPPM::gamlist(int a){

   return amlist[a];

}

/**
 * access to some static lists
 * @param i the index of the core
 * @param j the gaussian index
 * @return the maximal angular momentum of this type
 */
int CI_SPPM::gammax(int i,int j){

   return ammax[i][j];

}
