/*
 
   Program information
   
      ThING Is Not Gaussian is a rudimentary program to calculate
      matrix elements of gaussian contractions using recursion
      relations.
      
   Copyright notice

      Copyright 2011, 2011 Sebastian Wouters
      <sebastianwouters@gmail.com>
      
   Copyright permission
   
      This file is part of ThING Is Not Gaussian.

      ThING Is Not Gaussian is free software: you can redistribute
      it and/or modify it under the terms of the GNU General 
      Public License as published by the Free Software Foundation,
      either version 3 of the License, or (at your option) any
      later version.

      ThING Is Not Gaussian is distributed in the hope that it will
      be useful, but WITHOUT ANY WARRANTY; without even the implied
      warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
      See the GNU General Public License for more details.

      You should have received a copy of the GNU General Public
      License along with ThING Is Not Gaussian. If not, see 
      <http://www.gnu.org/licenses/>.
      
   File information MxElem.h and MxElem.cpp
   
      Class to fetch & store the matrix elements of a chemical
      problem. It contains:
         - the overlap matrix elements
         - the one body matrix elements
         - the kinetic energy matrix elements (so that nuclear charge
           rescaling can be used in HF if desired)
         - the teo body matrix elements

      Basic operations concerning the overlap matrix (Lodwin tfo,
      canonical tfo, inverse overlap matrix) are available
      "to a certain extent".
    
*/

#include "preamble.h"
#include "input.h"
#include "R.h"
#include "MxElem.h"
#include "MxElemFiller.h"


/**
 * Constructor for the MxElem class
 * @param readin the problem to be solved
 */
MxElem::MxElem(input & readin){

   this->N_Z = readin.gNcores();
   this->dim = CalcTotalNumberOfOrbitals(readin);

   //allocate the different matrixelements
   Telem = new double[(dim*(dim+1))/2];
   T = new double[(dim*(dim+1))/2];
   S = new double[(dim*(dim+1))/2];

   U = new double * [N_Z];

   for(int i = 0;i < N_Z;++i)
      U[i] = new double [dim*(dim + 1)/2];

   Velem = new double***[dim];

   for (int i=0; i<dim; i++){
      Velem[i] = new double**[dim-i+1];
      for (int x=0; x<dim-i+1; x++){
         Velem[i][x] = new double*[dim-i+1];
         for (int y=0; y<dim-i+1; y++){
            Velem[i][x][y] = new double[dim-i-x+1];
         }
      }
   }

   this->Init(readin);

}


/**
 * Copy constructor for the MxElem class
 * @param tocopy the MxElem object to be copied
 */
MxElem::MxElem(MxElem & tocopy){

   dim = tocopy.gdim();
   N_Z = tocopy.gN_Z();

   //allocate
   Telem = new double[(dim*(dim+1))/2];
   T = new double[(dim*(dim+1))/2];
   S = new double[(dim*(dim+1))/2];

   U = new double * [N_Z];

   for(int i = 0;i < N_Z;++i)
      U[i] = new double [dim*(dim + 1)/2];

   Velem = new double***[dim];

   for (int i=0; i<dim; i++){
      Velem[i] = new double**[dim-i+1];
      for (int x=0; x<dim-i+1; x++){
         Velem[i][x] = new double*[dim-i+1];
         for (int y=0; y<dim-i+1; y++){
            Velem[i][x][y] = new double[dim-i-x+1];
         }
      }
   }

   for (int i=0; i<dim; i++){
      for (int j=i; j<dim; j++){
         setTelem(i,j,tocopy.gTelem(i,j));
         setT(i,j,tocopy.gT(i,j));
         setS(i,j,tocopy.gS(i,j));
         for (int k=i; k<dim; k++)
            for (int l=j; l<dim; l++)
               setVelem(i,j,k,l,tocopy.gVelem(i,j,k,l));
      }
   }

}

/**
 * Standard destructor
 */
MxElem::~MxElem(){

   delete [] Telem;
   delete [] T;
   delete [] S;

   for(int i = 0;i < N_Z;++i)
      delete [] U[i];

   delete [] U;

   for (int i=0; i<dim; i++){
      for (int x=0; x<dim-i+1; x++){
         for (int y=0; y<dim-i+1; y++){
            delete [] Velem[i][x][y];
         }
         delete [] Velem[i][x];
      }
      delete [] Velem[i];
   }
   delete [] Velem;

}


/**
 * Function to find the orbital angular momentum corresponding to a subshell
 * @param type the type of the subshell: s, p, d, f, g, h
 * @return L the angular momentum
 */
int MxElem::GetLofType(char type){

   int L = -1;

   switch (type){
      case 's':
      case 'S':
         L = 0;
         break;
      case 'p':
      case 'P':
         L = 1;
         break;
      case 'd':
      case 'D':
         L = 2;
         break;
      case 'f':
      case 'F':
         L = 3;
         break;
      case 'g':
      case 'G':
         L = 4;
         break;
      case 'h':
      case 'H':
         L = 5;
         break;
      case 'i':
      case 'I':
         L = 6;
         break;
      case 'j':
      case 'J':
         L = 7;
         break;
      default:
         cout << "Subshell orbitals don't go higher than J in this program. Given : " << type << "." << endl;
         break;
   }

   assert(L>-1);
   return L;

}


/**
 * Function to find the total number of orbitals corresponding to a problem
 * @param readin the problem to be solved
 * @return the number of different orbitals
 */
int MxElem::CalcTotalNumberOfOrbitals(input & readin){

   int counter = 0;
   int Ncores = readin.gNcores();

   for (int cnt=0; cnt<Ncores; cnt++){

      Gauss * atom = readin.gGaussInfo(cnt);
      int Ntypes = atom->gNtypes();

      for (int cnt2=0; cnt2<Ntypes; cnt2++){

         char type = atom->gtype(cnt2);
         int L = GetLofType(type);
         counter += ((L+1)*(L+2))/2;

      }

   }

   return counter;

}


/**
 * Getter of the total number of orbitals
 * @return the total number of orbitals
 */
int MxElem::gdim() const {

   return dim;

}


/**
 * Getter of the one body matrix element (i|O|j)
 * @param i the first orbital
 * @param j the second orbital
 * @return (i|O|j)
 */
double MxElem::gTelem(int i, int j){

   if (i>j)
      return Telem[j + (i*(i+1))/2];

   return Telem[i + (j*(j+1))/2];

}


/**
 * Setter of the one body matrix element (i|O|j)
 * @param i the first orbital
 * @param j the second orbital
 * @param v the new value of the mx element
 */
void MxElem::setTelem(int i, int j, double v){

   if (i>j)
      Telem[j + (i*(i+1))/2] = v;
   else
      Telem[i + (j*(j+1))/2] = v;

}


/**
 * Getter of the kinetic energy matrix element (i|T|j)
 * @param i the first orbital
 * @param j the second orbital
 * @return (i|T|j)
 */
double MxElem::gT(int i, int j){

   if (i>j)
      return T[j + (i*(i+1))/2];

   return T[i + (j*(j+1))/2];

}


/**
 * Setter of the kinetic energy matrix element (i|T|j)
 * @param i the first orbital
 * @param j the second orbital
 * @param v the new value of the mx element
 */
void MxElem::setT(int i, int j, double v){

   if (i>j)
      T[j + (i*(i+1))/2] = v;
   else
      T[i + (j*(j+1))/2] = v;

}


/**
 * Getter of the overlap matrix element (i|j)
 * @param i the first orbital
 * @param j the second orbital
 * @return (i|j)
 */
double MxElem::gS(int i, int j){

   if (i>j)
      return S[j + (i*(i+1))/2];

   return S[i + (j*(j+1))/2];

}

/**
 * Getter of the nuclear attraction integral
 * @param core the index of the core
 * @param i the first orbital
 * @param j the second orbital
 * @return (i|j)
 */
double MxElem::gU(int core,int i, int j) {

   if (i>j)
      return U[core][j + (i*(i+1))/2];

   return U[core][i + (j*(j+1))/2];

}

/**
 * Setter of the overlap matrix element (i|j)
 * @param i the first orbital
 * @param j the second orbital
 * @param v the new value of the mx element
 */
void MxElem::setS(int i, int j, double v){

   if (i>j)
      S[j + (i*(i+1))/2] = v;
   else
      S[i + (j*(j+1))/2] = v;

}

/**
 * Setter of the nuclear attraction element
 @ param core index of the nuclear core
 * @param i the first orbital
 * @param j the second orbital
 * @param v the new value of the mx element
 */
void MxElem::sU(int core,int i, int j, double v){

   if (i>j)
      U[core][j + (i*(i+1))/2] = v;
   else
      U[core][i + (j*(j+1))/2] = v;

}

/**
 * Getter of the two body matrix element (ij|V|kl)
 * @param i the first orbital
 * @param j the second orbital
 * @param k the third orbital
 * @param l the fourth orbital
 * @return (ij|V|kl)
 */
double MxElem::gVelem(int i, int j, int k, int l){

   // i smallest
   if ((i<=j) && (i<=k) && (i<=l)){
      if (j<=l)
         return gVelemOK(i,j,k,l);
      else
         return gVelemOK(i,l,k,j);
   }
   //j smallest
   if ((j<=k) && (j<=l)){
      if (i<=k)
         return gVelemOK(j,i,l,k);
      else
         return gVelemOK(j,k,l,i);
   }
   //k smallest
   if (k<=l){
      if (j<=l)
         return gVelemOK(k,j,i,l);
      else
         return gVelemOK(k,l,i,j);
   }
   // l smallest
   if (i<=k)
      return gVelemOK(l,i,j,k);

   return gVelemOK(l,k,j,i);

}


/**
 * Setter of the two body matrix element (ij|V|kl)
 * @param i the first orbital
 * @param j the second orbital
 * @param k the third orbital
 * @param l the fourth orbital
 * @param v the new value of the mx element
 */
void MxElem::setVelem(int i, int j, int k, int l, double v){

   // i smallest
   if ((i<=j) && (i<=k) && (i<=l)){
      if (j<=l)
         setVelemOK(i,j,k,l,v);
      else
         setVelemOK(i,l,k,j,v);
   } else {
      //j smallest
      if ((j<=k) && (j<=l)){
         if (i<=k)
            setVelemOK(j,i,l,k,v);
         else
            setVelemOK(j,k,l,i,v);
      } else {
         //k smallest
         if (k<=l){
            if (j<=l)
               setVelemOK(k,j,i,l,v);
            else
               setVelemOK(k,l,i,j,v);
         } else {
            // l smallest
            if (i<=k)
               setVelemOK(l,i,j,k,v);
            else
               setVelemOK(l,k,j,i,v);
         }
      }
   }

}


/**
 * Getter of the two body matrix element (ij|V|kl) if i<j<l and i<k
 * @param i the first orbital
 * @param j the second orbital
 * @param k the third orbital
 * @param l the fourth orbital
 * @return (ij|V|kl)
 */
double MxElem::gVelemOK(int i, int j, int k, int l){

   return Velem[i][j-i][k-i][l-j];

}

/**
 * Setter of the two body matrix element (ij|V|kl) if i<j<l and i<k
 * @param i the first orbital
 * @param j the second orbital
 * @param k the third orbital
 * @param l the fourth orbital
 * @param v the new value of the mx element
 */
void MxElem::setVelemOK(int i, int j, int k, int l, double v){

   //cout << "dim = " << dim << " and (" << i << "," << j << "," << k << "," << l << ")" << endl;
   Velem[i][j-i][k-i][l-j] = v;

}

/**
 * Calculate the nuclear-nuclear potential energy of a problem
 * @param problem the problem to be solved
 */
double MxElem::NuclPotEn(input & problem){

   double energy = 0.0;
   int Ncores = problem.gNcores();

   for (int i=0; i<Ncores; i++){

      R Ri(*(problem.gvector(i)));
      int Zi = problem.gcore(i);

      for (int j=i+1; j<Ncores; j++){

         energy += Zi * problem.gcore(j) / sqrt(Ri.DistanceSquared(*(problem.gvector(j))));
      }
   }

   return energy;

}

/**
 * Initialise the matrix elements
 * @param readin the problem to be solved
 */
void MxElem::Init(input & readin){

   MxElemFiller filler(readin);

   double centerx = 0.0;
   double centery = 0.0;
   double centerz = 0.0;

   for (int i=0; i<readin.gNcores(); i++){
      R Ratom(*(readin.gvector(i)));
      centerx += Ratom.gxco();
      centery += Ratom.gyco();
      centerz += Ratom.gzco();
   }
   centerx = centerx/readin.gNcores();
   centery = centery/readin.gNcores();
   centerz = centerz/readin.gNcores();

   R Center(centerx,centery,centerz);

   int cnt1 = -1;

   int * CoreNumber = new int[dim];
   int * Type = new int[dim];
   int * n1xvalues = new int[dim];
   int * n1yvalues = new int[dim];
   int * n1zvalues = new int[dim];

   for (int i=0; i<readin.gNcores(); i++){
      Gauss * first = readin.gGaussInfo(i);

      for (int k=0; k<(first->gNtypes()); k++){
         int L1 = GetLofType(first->gtype(k));

         // loop over different [n1x,n1y,n1z] contributions with n1x+n1y+n1z=L1;
         for (int n1x=L1; n1x>=0; n1x--){
            for (int n1y=L1-n1x; n1y>=0; n1y--){
               int n1z = L1-n1x-n1y;
               cnt1++;

               CoreNumber[cnt1] = i;
               Type[cnt1] = k;
               n1xvalues[cnt1] = n1x;
               n1yvalues[cnt1] = n1y;
               n1zvalues[cnt1] = n1z;

            }
         }
      }
   }

   for (int count1 = 0; count1<dim; count1++){

      int i = CoreNumber[count1];
      int k = Type[count1];
      int n1x = n1xvalues[count1];
      int n1y = n1yvalues[count1];
      int n1z = n1zvalues[count1];

      double OneBodyElement;
      int count2, count3, count4, start, start2, start3;

      count2 = -1;

      for (int j=i; j<readin.gNcores(); j++){
         Gauss * second = readin.gGaussInfo(j);
         start = 0;
         if (i==j) start=k;

         for (int l=start; l<(second->gNtypes()); l++){
            int L2 = GetLofType(second->gtype(l));

            // loop over different [n2x,n2y,n2z] contributions with n2x+n2y+n2z=L2;
            for (int n2x=L2; n2x>=0; n2x--){
               for (int n2y=L2-n2x; n2y>=0; n2y--){
                  int n2z = L2-n2x-n2y;

                  if ((i==j)&&(l==k)){
                     if (n2x>n1x){
                        n2x=n1x;
                        n2y=n1y;
                        n2z=n1z;
                     } else {
                        if (n2x==n1x){
                           if (n2y>n1y){
                              n2y=n1y;
                              n2z=n1z;
                           }
                        }
                     }
                  }

                  count2++;
                  //If they're centered on the same atom and n1i+n2i odd for i=x,y or z -> Overlap & T 0.0
                  if ((i==j)&&((((n1x+n2x)%2)!=0)||(((n1y+n2y)%2)!=0)||(((n1z+n2z)%2)!=0))){
                     setS(count1,count1+count2,0.0);
                     setT(count1,count1+count2,0.0);
                     OneBodyElement = 0.0;
                  } else {
                     setS(count1,count1+count2,filler.Overlap(i,k,n1x,n1y,n1z,j,l,n2x,n2y,n2z));
                     OneBodyElement = filler.KE(i,k,n1x,n1y,n1z,j,l,n2x,n2y,n2z);
                     setT(count1,count1+count2,OneBodyElement);
                  }

                  for (int Ncore=0; Ncore<readin.gNcores(); Ncore++){

                     double ward = filler.ElNucl(i,k,n1x,n1y,n1z,j,l,n2x,n2y,n2z,Ncore);

                     sU(Ncore,count1,count1+count2,ward);

                     OneBodyElement += ward;

                  }

                  setTelem(count1,count1+count2,OneBodyElement);

                  count3 = -1;

                  for (int m=i; m<readin.gNcores(); m++){
                     Gauss * third = readin.gGaussInfo(m);
                     start2 = 0;
                     if (i==m) start2=k;

                     for (int n=start2; n<(third->gNtypes()); n++){
                        int L3 = GetLofType(third->gtype(n));
                        for (int n3x=L3; n3x>=0; n3x--){
                           for (int n3y=L3-n3x; n3y>=0; n3y--){
                              int n3z = L3-n3x-n3y;

                              if ((i==m)&&(n==k)){
                                 if (n3x>n1x){
                                    n3x=n1x;
                                    n3y=n1y;
                                    n3z=n1z;
                                 } else {
                                    if (n3x==n1x){
                                       if (n3y>n1y){
                                          n3y=n1y;
                                          n3z=n1z;
                                       }
                                    }
                                 }
                              }

                              count3++;

                              count4 = -1;
                              for (int o=j; o<readin.gNcores(); o++){
                                 Gauss * fourth = readin.gGaussInfo(o);
                                 start3 = 0;
                                 if (j==o) start3=l;

                                 for (int p=start3; p<(fourth->gNtypes()); p++){
                                    int L4 = GetLofType(fourth->gtype(p));
                                    for (int n4x=L4; n4x>=0; n4x--){
                                       for (int n4y=L4-n4x; n4y>=0; n4y--){
                                          int n4z = L4-n4x-n4y;

                                          if ((j==o)&&(p==l)){
                                             if (n4x>n2x){
                                                n4x=n2x;
                                                n4y=n2y;
                                                n4z=n2z;
                                             } else {
                                                if (n4x==n2x){
                                                   if (n4y>n2y){
                                                      n4y=n2y;
                                                      n4z=n2z;
                                                   }
                                                }
                                             }
                                          }

                                          count4++;

                                          //alpha=(i,k,1) beta=(j,l,2) gamma=(m,n,3) delta=(o,p,4)
                                          //(alpha beta|V|gamma delta) = (alpha gamma beta delta) in MxElemFiller!!
                                          setVelem(count1, count1+count2, count1+count3, count1+count2+count4, filler.ElEl(i,k,n1x,n1y,n1z,m,n,n3x,n3y,n3z,j,l,n2x,n2y,n2z,o,p,n4x,n4y,n4z));
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

   delete [] CoreNumber;
   delete [] Type;
   delete [] n1xvalues;
   delete [] n1yvalues;
   delete [] n1zvalues;

}

/**
 * @return the number of cores
 */
int MxElem::gN_Z() const {

   return N_Z;

}
