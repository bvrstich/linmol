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

vector< vector<int> > CI_SPM::s2inlxyz;
int ******CI_SPM::inlxyz2s;

vector< vector<int> > CI_SPM::s2gauss;

int CI_SPM::l_max;
int CI_SPM::n_max;

int CI_SPM::dim;

/** 
 * static function that reads in the input data and makes the corresponding lists and dimensions
 */
void CI_SPM::init(){

   vector<int> v(6);
   vector<int> vg(2);

   int curtyp = 0;

   n_max = 0;
   l_max = 0;

   for(int i = 0;i < input::gN_Z();++i){

      v[0] = i;

      for(int j = 0;j < input::gGaussInfo(i).gNtypes();++j){

         if(input::gGaussInfo(i).gtype(j) == 'S')
            v[2] = 0;
         else if(input::gGaussInfo(i).gtype(j) == 'P')
            v[2] = 1;
         else if(input::gGaussInfo(i).gtype(j) == 'D')
            v[2] = 2;
         else if(input::gGaussInfo(i).gtype(j) == 'F')
            v[2] = 3;
         else if(input::gGaussInfo(i).gtype(j) == 'G')
            v[2] = 4;
         else if(input::gGaussInfo(i).gtype(j) == 'H')
            v[2] = 5;
         else
            cout << "BASISSET TOO LARGE" << endl;

         if(j == 0){

            v[1] = v[2] + 1;
            curtyp = v[2];

         }
         else{

            if(v[2] == curtyp)
               v[1]++;
            else{

               v[1] = v[2] + 1;
               curtyp = v[2];

            }

         }

         if(v[1] > n_max)
            n_max = v[1];

         if(v[2] > l_max)
            l_max = v[2];

         for(int x = v[2];x >= 0;x--)
            for(int y = v[2] - x;y >= 0;y--){

               v[3] = x;
               v[4] = y;
               v[5] = v[2] - x - y;

               s2inlxyz.push_back(v);

               vg[0] = i;
               vg[1] = j;

               s2gauss.push_back(vg);

            }

      }

   }
   
   //allocate the list
   inlxyz2s =  new int ***** [input::gN_Z()];

   for(int i = 0;i < input::gN_Z();++i){

      inlxyz2s[i] = new int **** [n_max];

      for(int n = 0;n < n_max;++n){

         inlxyz2s[i][n] =  new int *** [l_max + 1];

         for(int l = 0;l <= l_max;++l){

            inlxyz2s[i][n][l] =  new int ** [l + 1];

            for(int x = 0;x <= l;++x){

               inlxyz2s[i][n][l][x] =  new int * [l + 1];

               for(int y = 0;y <= l;++y)
                  inlxyz2s[i][n][l][x][y] =  new int [l + 1];

            }

         }

      }

   }

   //fill list using other list
   for(unsigned int s = 0;s < s2inlxyz.size();++s){

      v = s2inlxyz[s];

      inlxyz2s[v[0]][v[1] - v[2] - 1][v[2]][v[3]][v[4]][v[5]] = s;

   }

   dim = s2inlxyz.size();

}

/** 
 * function that deallocates the static members
 */
void CI_SPM::clear(){

   for(int i = 0;i < input::gN_Z();++i){

      for(int n = 0;n < n_max;++n){

         for(int l = 0;l <= l_max;++l){

            for(int x = 0;x <= l;++x){

               for(int y = 0;y <= l;++y)
                  delete [] inlxyz2s[i][n][l][x][y];

               delete [] inlxyz2s[i][n][l][x];

            }

            delete [] inlxyz2s[i][n][l];

         }

         delete [] inlxyz2s[i][n];

      }

      delete [] inlxyz2s[i];

   }

   delete [] inlxyz2s;

}

/** 
 * Standard constructor: allocates the matrices and fills them with the correct elements
 */
CI_SPM::CI_SPM() : Matrix(dim) { }

/** 
 * copy constructor
 * @param ci_c CI_SPM object to be copied in the newly constructed object
 */
CI_SPM::CI_SPM(const CI_SPM &ci_c) : Matrix(ci_c){ }

/**
 * standard destructor
 */
CI_SPM::~CI_SPM(){ }

/**
 * access the elements of the matrix using the physical quantum number i,n,l,x,y,z.
 * @return the number on place CI_SPM(i,j)
 */
double CI_SPM::operator()(int i_a,int n_a,int l_a,int x_a,int y_a,int z_a,int i_b,int n_b,int l_b,int x_b,int y_b,int z_b) const {

   return (*this)(ginlxyz2s(i_a,n_a,l_a,x_a,y_a,z_a),ginlxyz2s(i_b,n_b,l_b,x_b,y_b,z_b));

}

ostream &operator<<(ostream &output,CI_SPM &ci_p){

   for(int i = 0;i < ci_p.dim;++i)
      for(int j = 0;j < ci_p.dim;++j){

         output << ci_p.s2inlxyz[i][0] << "\t" << ci_p.s2inlxyz[i][1] << "\t" << ci_p.s2inlxyz[i][2]
         
            << "\t(" << ci_p.s2inlxyz[i][3] << "," << ci_p.s2inlxyz[i][4] << "," << ci_p.s2inlxyz[i][5] << ")\t|\t"

            << ci_p.s2inlxyz[j][0] << "\t" << ci_p.s2inlxyz[j][1] << "\t" << ci_p.s2inlxyz[j][2]
         
            << "\t(" << ci_p.s2inlxyz[j][3] << "," << ci_p.s2inlxyz[j][4] << "," << ci_p.s2inlxyz[j][5] << ")\t|\t" << ci_p(i,j) << endl;

      }

   return output;

}

/**
 * public static function to access the lists safely from other classes
 * @param s the single particle index
 * @param option indicates what quantumnumber will be returned
 * @return option == 0: i , option == 1 : n, option == 2 : l, option == 3 : x, option == 4 : y , option == 5 : z
 */
int CI_SPM::gs2inlxyz(int s,int option) {

   return s2inlxyz[s][option];

}

/**
 * public static function to access the lists safely from other classes
 * @param s the single particle index
 * @param option indicates what gaussian index will be returned
 * @return option == 0: i , option == 1 : j (the index of the Gaussian orbital on core i corresponding to s)
 */
int CI_SPM::gs2gauss(int s,int option) {

   return s2gauss[s][option];

}

/**
 * public static function to access the lists safely from other classes
 * @param i the i'th core
 * @param n the main quantumnumber
 * @param l the angular momentum
 * @param x the power of x in the Gaussian wavefucntion
 * @param y the power of y in the Gaussian wavefucntion
 * @param z the power of z in the Gaussian wavefucntion
 * @return s the single particle index corresponding to i n l x y z quantumnumbers
 */
int CI_SPM::ginlxyz2s(int i,int n,int l,int x,int y,int z) {

   return inlxyz2s[i][n - l - 1][l][x][y][z];

}

/**
 * static function that returns the dimension
 * @return the dimension of the basisset
 */
int CI_SPM::gdim(){

   return dim;

}

/**
 * static function
 * @return highest main quantumnumber in the basisset
 */
int CI_SPM::gn_max(){

   return n_max;

}

/**
 * static function
 * @return highest angular momentum quantumnumber in the basisset
 */
int CI_SPM::gl_max(){

   return l_max;

}

/**
 * Fill the CI_SPM object with the overlap matrix
 */
void CI_SPM::calc_overlap(){

   int i_a,x_a,y_a,z_a;
   int i_b,x_b,y_b,z_b;

   int j_a,j_b;//gaussian indices

   for(int a = 0;a < dim;++a){

      i_a = gs2inlxyz(a,0);
      j_a = gs2gauss(a,1);

      x_a = gs2inlxyz(a,3);
      y_a = gs2inlxyz(a,4);
      z_a = gs2inlxyz(a,5);

      for(int b = a;b < dim;++b){

         i_b = gs2inlxyz(b,0);
         j_b = gs2gauss(b,1);

         x_b = gs2inlxyz(b,3);
         y_b = gs2inlxyz(b,4);
         z_b = gs2inlxyz(b,5);

         (*this)(a,b) = 0.0;

         for(int k_a = 0;k_a < input::gGaussInfo(i_a).gNcontr(j_a);++k_a)
            for(int k_b = 0;k_b < input::gGaussInfo(i_b).gNcontr(j_b);++k_b){

               (*this)(a,b) += input::gGaussInfo(i_a).gprefactors(j_a,k_a) * input::gGaussInfo(i_b).gprefactors(j_b,k_b) 
               
                  * Tools::gnorm(a,k_a) * Tools::gnorm(b,k_b) * Tools::gS()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

            }

      }
   }

   this->symmetrize();

}

/**
 * Fill the CI_SPM object with the kinetic energy matrix
 */
void CI_SPM::calc_kinetic(){

   int i_a,x_a,y_a,z_a;
   int i_b,x_b,y_b,z_b;

   int j_a,j_b;//gaussian indices

   for(int a = 0;a < dim;++a){

      i_a = gs2inlxyz(a,0);
      j_a = gs2gauss(a,1);

      x_a = gs2inlxyz(a,3);
      y_a = gs2inlxyz(a,4);
      z_a = gs2inlxyz(a,5);

      for(int b = a;b < dim;++b){

         i_b = gs2inlxyz(b,0);
         j_b = gs2gauss(b,1);

         x_b = gs2inlxyz(b,3);
         y_b = gs2inlxyz(b,4);
         z_b = gs2inlxyz(b,5);

         (*this)(a,b) = 0.0;

         for(int k_a = 0;k_a < input::gGaussInfo(i_a).gNcontr(j_a);++k_a)
            for(int k_b = 0;k_b < input::gGaussInfo(i_b).gNcontr(j_b);++k_b){

               (*this)(a,b) += input::gGaussInfo(i_a).gprefactors(j_a,k_a) * input::gGaussInfo(i_b).gprefactors(j_b,k_b) 
               
                  * Tools::gnorm(a,k_a) * Tools::gnorm(b,k_b) * Tools::gT()(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b);

            }

      }
   }

   this->symmetrize();

}

/**
 * Fill the CI_SPM object with the nuclear attraction energy matrixelements to core "core"
 */
void CI_SPM::calc_nucl_attrac(int core){

   int i_a,x_a,y_a,z_a;
   int i_b,x_b,y_b,z_b;

   int j_a,j_b;//gaussian indices

   for(int a = 0;a < dim;++a){

      i_a = gs2inlxyz(a,0);
      j_a = gs2gauss(a,1);

      x_a = gs2inlxyz(a,3);
      y_a = gs2inlxyz(a,4);
      z_a = gs2inlxyz(a,5);

      for(int b = a;b < dim;++b){

         i_b = gs2inlxyz(b,0);
         j_b = gs2gauss(b,1);

         x_b = gs2inlxyz(b,3);
         y_b = gs2inlxyz(b,4);
         z_b = gs2inlxyz(b,5);

         (*this)(a,b) = 0.0;

         for(int k_a = 0;k_a < input::gGaussInfo(i_a).gNcontr(j_a);++k_a)
            for(int k_b = 0;k_b < input::gGaussInfo(i_b).gNcontr(j_b);++k_b){

               (*this)(a,b) += input::gGaussInfo(i_a).gprefactors(j_a,k_a) * input::gGaussInfo(i_b).gprefactors(j_b,k_b) 
               
                  * Tools::gnorm(a,k_a) * Tools::gnorm(b,k_b) * Tools::gU(core)(i_a,j_a,k_a,x_a,y_a,z_a,i_b,j_b,k_b,x_b,y_b,z_b,0);

            }

         (*this)(a,b) *= -input::gZ(core);

      }
   }

   this->symmetrize();

}

/* vim: set ts=3 sw=3 expandtab :*/
