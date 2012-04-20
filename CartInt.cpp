#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::vector;

#include "include.h"

input *CartInt::readin;
vector< vector<int> > CartInt::s2inlxyz;
int ******CartInt::inlxyz2s;

int CartInt::l_max;
int CartInt::n_max;
int CartInt::N_Z;

vector< vector<int> > CartInt::t2s;
int **CartInt::s2t;

int CartInt::dim;
int CartInt::N;

double CartInt::NucRepEn;

/** 
 * static function that reads in the input data and makes the matrix elements
 */
void CartInt::init(){

   readin = new input("start.stp");

   vector<int> v(6);

   int curtyp = 0;

   n_max = 0;
   l_max = 0;
   N_Z = readin->gNcores();

   N = readin->NumberOfElectrons();

   for(int i = 0;i < readin->gNcores();++i){

      v[0] = i;

      for(int j = 0;j < readin->gGaussInfo(i)->gNtypes();++j){

         if(readin->gGaussInfo(i)->gtype(j) == 'S')
            v[2] = 0;
         else if(readin->gGaussInfo(i)->gtype(j) == 'P')
            v[2] = 1;
         else if(readin->gGaussInfo(i)->gtype(j) == 'D')
            v[2] = 2;
         else if(readin->gGaussInfo(i)->gtype(j) == 'F')
            v[2] = 3;
         else if(readin->gGaussInfo(i)->gtype(j) == 'G')
            v[2] = 4;
         else if(readin->gGaussInfo(i)->gtype(j) == 'H')
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

            }

      }

   }
   
   //allocate the list
   inlxyz2s =  new int ***** [readin->gNcores()];

   for(int i = 0;i < readin->gNcores();++i){

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

   s2t = new int * [dim];

   for(int i = 0;i < dim;++i)
      s2t[i] = new int [dim];

   vector<int> vst(2);

   int iter = 0;

   for(int i = 0;i < dim;++i)
      for(int j = 0;j < dim;++j){

         vst[0] = i;
         vst[1] = j;

         t2s.push_back(vst);

         s2t[i][j] = iter;

         ++iter;

      }

   NucRepEn = readin->NucRepEn();

}

/** 
 * function that deallocates the static members
 */
void CartInt::clear(){

   for(int i = 0;i < readin->gNcores();++i){

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

   for(int i = 0;i < dim;++i)
      delete [] s2t[i];

   delete [] s2t;

   delete readin;

}

/** 
 * Standard constructor: allocates the matrices and fills them with the correct elements
 */
CartInt::CartInt(){ 

   S = new Matrix(dim);
   T = new Matrix(dim);
   U = new Matrix(dim);

   V = new Matrix(dim*dim);

   MxElem setup(*readin);
   setup.Init(*readin);

   for(int i = 0;i < dim;++i)
      for(int j = i;j < dim;++j){

         (*S)(i,j) = setup.gSoverlap(i,j);
         (*T)(i,j) = setup.gKEseparate(i,j);
         (*U)(i,j) = setup.gTelem(i,j) - setup.gKEseparate(i,j);

      }

   S->symmetrize();
   T->symmetrize();
   U->symmetrize();

   int a,b,c,d;

   for(int i = 0;i < dim*dim;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < dim*dim;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*V)(i,j) = setup.gVelem(a,b,c,d);

      }
   }

   V->symmetrize();
   
}

/** 
 * copy constructor
 * @param ci_c CartInt object to be copied in the newly constructed object
 */
CartInt::CartInt(const CartInt &ci_c){ 

   S = new Matrix(ci_c.gS());
   T = new Matrix(ci_c.gT());
   U = new Matrix(ci_c.gU());
   
   V = new Matrix(ci_c.gV());

}

/**
 * standard destructor
 */
CartInt::~CartInt(){ 

   delete S;
   delete T;
   delete U;

   delete V;
   
}

/** 
 * @return the overlapmatrix, const version
 */
const Matrix &CartInt::gS() const { 

   return *S;

}

/** 
 * @return the overlapmatrix
 */
Matrix &CartInt::gS() { 

   return *S;

}

/** 
 * @return the kinetic energy matrix, const version
 */
const Matrix &CartInt::gT() const { 

   return *T; 
}

/** 
 * @return the kinetic energy matrix
 */
Matrix &CartInt::gT() { 

   return *T;

}

/** 
 * @return the nuclear attraction matrix, const version
 */
const Matrix &CartInt::gU() const { 

   return *U; 
}

/** 
 * @return the nuclear attraction matrix
 */
Matrix &CartInt::gU() { 

   return *U;

}

/** 
 * @return the electronic repulsion matrix
 */
const Matrix &CartInt::gV() const { 

   return *V; 
}

/** 
 * @return the electronic repulsion matrix
 */
Matrix &CartInt::gV() { 

   return *V;

}

/**
 * normalize the cartesian wavefunctions, needed for transformation to spherical
 */
void CartInt::norm() {

   double norm[dim];

   for(int i = 0;i < dim;++i)
      norm[i] = sqrt((*S)(i,i));

   //first the one body matrices
   for(int i = 0;i < dim;++i)
      for(int j = i;j < dim;++j){

         (*S)(i,j) /= norm[i] * norm[j];
         (*T)(i,j) /= norm[i] * norm[j];
         (*U)(i,j) /= norm[i] * norm[j];

      }

   int a,b,c,d;

   //then the two body matrices
   for(int i = 0;i < dim*dim;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < dim*dim;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*V)(i,j) /= norm[a] * norm[b] * norm[c] * norm[d];

      }

   }

   S->symmetrize();
   T->symmetrize();
   U->symmetrize();
   V->symmetrize();

}

ostream &operator<<(ostream &output,CartInt &ci_p){

   output << endl;
   output << "Overlap Matrix:" << endl;
   output << endl;

   //first overlap
   for(int s_i = 0;s_i < ci_p.dim;++s_i)
      for(int s_j = 0;s_j < ci_p.dim;++s_j){

         output << ci_p.s2inlxyz[s_i][0] << "\t" << ci_p.s2inlxyz[s_i][1] << "\t" << ci_p.s2inlxyz[s_i][2]
         
            << "\t(" << ci_p.s2inlxyz[s_i][3] << "," << ci_p.s2inlxyz[s_i][4] << "," << ci_p.s2inlxyz[s_i][5] << ")\t|\t"

            << ci_p.s2inlxyz[s_j][0] << "\t" << ci_p.s2inlxyz[s_j][1] << "\t" << ci_p.s2inlxyz[s_j][2]
         
            << "\t(" << ci_p.s2inlxyz[s_j][3] << "," << ci_p.s2inlxyz[s_j][4] << "," << ci_p.s2inlxyz[s_j][5] << ")\t|\t" << (ci_p.gS())(s_i,s_j) << endl;

      }

   output << endl;
   output << "Kinetic energy:" << endl;
   output << endl;

   //kinetic energy
   for(int s_i = 0;s_i < ci_p.dim;++s_i)
      for(int s_j = 0;s_j < ci_p.dim;++s_j){

         output << ci_p.s2inlxyz[s_i][0] << "\t" << ci_p.s2inlxyz[s_i][1] << "\t" << ci_p.s2inlxyz[s_i][2]
         
            << "\t(" << ci_p.s2inlxyz[s_i][3] << "," << ci_p.s2inlxyz[s_i][4] << "," << ci_p.s2inlxyz[s_i][5] << ")\t|\t"

            << ci_p.s2inlxyz[s_j][0] << "\t" << ci_p.s2inlxyz[s_j][1] << "\t" << ci_p.s2inlxyz[s_j][2]
         
            << "\t(" << ci_p.s2inlxyz[s_j][3] << "," << ci_p.s2inlxyz[s_j][4] << "," << ci_p.s2inlxyz[s_j][5] << ")\t|\t" << (ci_p.gT())(s_i,s_j) << endl;

      }

   output << endl;
   output << "Nuclear attraction energy:" << endl;
   output << endl;

   //kinetic energy
   for(int s_i = 0;s_i < ci_p.dim;++s_i)
      for(int s_j = 0;s_j < ci_p.dim;++s_j){

         output << ci_p.s2inlxyz[s_i][0] << "\t" << ci_p.s2inlxyz[s_i][1] << "\t" << ci_p.s2inlxyz[s_i][2]

            << "\t(" << ci_p.s2inlxyz[s_i][3] << "," << ci_p.s2inlxyz[s_i][4] << "," << ci_p.s2inlxyz[s_i][5] << ")\t|\t"

            << ci_p.s2inlxyz[s_j][0] << "\t" << ci_p.s2inlxyz[s_j][1] << "\t" << ci_p.s2inlxyz[s_j][2]

            << "\t(" << ci_p.s2inlxyz[s_j][3] << "," << ci_p.s2inlxyz[s_j][4] << "," << ci_p.s2inlxyz[s_j][5] << ")\t|\t" << (ci_p.gU())(s_i,s_j) << endl;

      }

   output << endl;
   output << "Electronic repulsion energy:" << endl;
   output << endl;

   int s_i,s_j,s_k,s_l;

   for(int t_i = 0;t_i < ci_p.dim*ci_p.dim;++t_i){

      s_i = ci_p.t2s[t_i][0];
      s_j = ci_p.t2s[t_i][1];

      for(int t_j = 0;t_j < ci_p.dim*ci_p.dim;++t_j){

         s_k = ci_p.t2s[t_j][0];
         s_l = ci_p.t2s[t_j][1];

         output << "[\t" << ci_p.s2inlxyz[s_i][0] << "\t" << ci_p.s2inlxyz[s_i][1] << "\t" << ci_p.s2inlxyz[s_i][2]

            << "\t(" << ci_p.s2inlxyz[s_i][3] << "," << ci_p.s2inlxyz[s_i][4] << "," << ci_p.s2inlxyz[s_i][5] << ")\t|\t"

            << ci_p.s2inlxyz[s_j][0] << "\t" << ci_p.s2inlxyz[s_j][1] << "\t" << ci_p.s2inlxyz[s_j][2]

            << "\t(" << ci_p.s2inlxyz[s_j][3] << "," << ci_p.s2inlxyz[s_j][4] << "," << ci_p.s2inlxyz[s_j][5] << ")\t]\t||\t[" 

            << ci_p.s2inlxyz[s_k][0] << "\t" << ci_p.s2inlxyz[s_k][1] << "\t" << ci_p.s2inlxyz[s_k][2]

            << "\t(" << ci_p.s2inlxyz[s_k][3] << "," << ci_p.s2inlxyz[s_k][4] << "," << ci_p.s2inlxyz[s_k][5] << ")\t|\t"

            << ci_p.s2inlxyz[s_l][0] << "\t" << ci_p.s2inlxyz[s_l][1] << "\t" << ci_p.s2inlxyz[s_l][2]

            << "\t(" << ci_p.s2inlxyz[s_l][3] << "," << ci_p.s2inlxyz[s_l][4] << "," << ci_p.s2inlxyz[s_l][5] << ")\t]\t" <<(ci_p.gV())(t_i,t_j) << endl;

      }
   }

   return output;

}


/**
 * public static function to access the lists safely from other classes
 * @param s the single particle index
 * @param option indicates what quantumnumber will be returned
 * @return option == 0: i , option == 1 : n, option == 2 : l, option == 3 : x, option == 4 : y , option == 5 : z
 */
int CartInt::gs2inlxyz(int s,int option) {

   return s2inlxyz[s][option];

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
int CartInt::ginlxyz2s(int i,int n,int l,int x,int y,int z) {

   return inlxyz2s[i][n - l - 1][l][x][y][z];

}

/**
 * static function that returns the dimension
 * @return the dimension of the basisset
 */
int CartInt::gdim(){

   return dim;

}

/**
 * static function that returns the dimension
 * @return nr of cores N_Z
 */
int CartInt::gN_Z(){

   return N_Z;

}

/**
 * static function
 * @return highest main quantumnumber in the basisset
 */
int CartInt::gn_max(){

   return n_max;

}

/**
 * static function
 * @return highest angular momentum quantumnumber in the basisset
 */
int CartInt::gl_max(){

   return l_max;

}

/**
 * static function
 * @return nr of electrons
 */
int CartInt::gN(){

   return N;

}

/**
 * orthogonalizes the basis: inverse sqrt of S
 */
void CartInt::orthogonalize() {

   //first inverse sqrt of S
   S->sqrt(-1);

   Matrix T_copy(dim);
   Matrix U_copy(dim);

   T_copy = 0.0;
   U_copy = 0.0;

   //transform T
   for(int i = 0;i < dim;++i)
      for(int j = 0;j < dim;++j){

         for(int k = 0;k < dim;++k){

            T_copy(i,j) += (*S)(i,k) * (*T)(k,j);
            U_copy(i,j) += (*S)(i,k) * (*U)(k,j);

         }

      }

   *T = 0.0;
   *U = 0.0;

   for(int i = 0;i < dim;++i)
      for(int j = i;j < dim;++j){

         for(int k = 0;k < dim;++k){

            (*T)(i,j) += T_copy(i,k) * (*S)(k,j);
            (*U)(i,j) += U_copy(i,k) * (*S)(k,j);

         }

      }

   Matrix V_copy(dim*dim);

   V_copy = 0.0;

   int a,b,c,d;

   //contract a
   for(int i = 0;i < dim*dim;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = 0;j < dim*dim;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         for(int a_ = 0;a_ < dim;++a_)
            V_copy(i,j) += (*S)(a,a_) * (*V)(s2t[a_][b],j); 

      }
   }

   *V = 0.0;

   //contract b
   for(int i = 0;i < dim*dim;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = 0;j < dim*dim;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         for(int b_ = 0;b_ < dim;++b_)
            (*V)(i,j) += (*S)(b,b_) * V_copy(s2t[a][b_],j); 

      }
   }

   V_copy = 0.0;

   //contract c
   for(int i = 0;i < dim*dim;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = 0;j < dim*dim;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         for(int c_ = 0;c_ < dim;++c_)
            V_copy(i,j) += (*V)(i,s2t[c_][d]) * (*S)(c_,c); 

      }
   }

   *V = 0.0;

   //contract d
   for(int i = 0;i < dim*dim;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = 0;j < dim*dim;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         for(int d_ = 0;d_ < dim;++d_)
            (*V)(i,j) += V_copy(i,s2t[c][d_]) * (*S)(d_,d); 

      }
   }

   T->symmetrize();
   U->symmetrize();

}

/**
 * access to the individual elements of the matrices
 */
double CartInt::gS(int i,int j) const {

   return (*S)(i,j);

}

/**
 * access to the individual elements of the matrices
 */
double CartInt::gT(int i,int j) const {

   return (*T)(i,j);

}

/**
 * access to the individual elements of the matrices
 */
double CartInt::gU(int i,int j) const {

   return (*U)(i,j);

}

/**
 * access to the individual elements of the matrices
 */
double CartInt::gV(int a,int b,int c,int d) const {

   return (*V)(s2t[a][b],s2t[c][d]);

}

/**
 * static function that returns the nuclear repulsion energy for this configuration
 */
double CartInt::gNucRepEn() {

   return NucRepEn;

}
