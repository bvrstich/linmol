#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::complex;

#include "include.h"

/** 
 * Standard constructor: 
 * fills the coefficients and index arrays with the correct CartInt indices corresponding to the spherical i,n,l and m quantumnumbers
 * @param i the index of the core
 * @param n the main quantumnumber
 * @param l the angular momentum
 * @param m the angular momentum projection
 */
Transform::Transform(int i,int n,int l,int m){ 

   this->i = i;
   this->n = n;
   this->l = l;
   this->m = m;

   if(l == 0){

      dim = 1;

      coef = new complex<double> * [dim];
      ind = new int [dim];

      coef[0] = new complex<double>(1.0,0.0);
      ind[0] = CartInt::ginlxyz2s(i,n,l,0,0,0);

   }
   else if(l == 1){

      if(m == 1){

         dim = 2;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(1.0/sqrt(2.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,1,0,0);

         coef[1] = new complex<double>(0.0,1.0/sqrt(2.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,1,0);

      }
      else if(m == 0){

         dim = 1;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         coef[0] = new complex<double>(1.0,0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,0,0,1);

      }
      else if(m == -1){

         dim = 2;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(1.0/sqrt(2.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,1,0,0);

         coef[1] = new complex<double>(0.0,-1.0/sqrt(2.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,1,0);

      }

   }
   else if(l == 2){

      if(m == 2){

         dim = 3;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(3.0/8.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,2,0,0);

         coef[1] = new complex<double>(-sqrt(3.0/8.0),0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,2,0);

         coef[2] = new complex<double>(0.0,1.0/sqrt(2.0));
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,1,0);

      }
      else if(m == 1){

         dim = 2;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(1.0/sqrt(2.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,1,0,1);

         coef[1] = new complex<double>(0.0,1.0/sqrt(2.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,1,1);

      }
      else if(m == 0){

         dim = 3;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(1.0,0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,0,0,2);

         coef[1] = new complex<double>(-0.5,0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,2,0,0);

         coef[2] = new complex<double>(-0.5,0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,0,2,0);

      }
      else if(m == -1){

         dim = 2;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(1.0/sqrt(2.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,1,0,1);

         coef[1] = new complex<double>(0.0,-1.0/sqrt(2.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,1,1);

      }
      else if(m == -2){

         dim = 3;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(3.0/8.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,2,0,0);

         coef[1] = new complex<double>(-sqrt(3.0/8.0),0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,2,0);

         coef[2] = new complex<double>(0.0,-1.0/sqrt(2.0));
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,1,0);

      }

   }
   else if(l == 3){

      if(m == 3){

         dim = 4;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(5.0)/4.0,0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,3,0,0);

         coef[1] = new complex<double>(0.0,-sqrt(5.0)/4.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,3,0);

         coef[2] = new complex<double>(-0.75,0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,2,0);

         coef[3] = new complex<double>(0.0,0.75);
         ind[3] = CartInt::ginlxyz2s(i,n,l,2,1,0);

      }
      else if(m == 2){

         dim = 3;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         coef[0] = new complex<double>(sqrt(3.0/8.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,2,0,1);

         coef[1] = new complex<double>(-sqrt(3.0/8.0),0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,2,1);

         coef[2] = new complex<double>(0.0,1.0/sqrt(2.0));
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,1,1);

      }
      else if(m == 1){

         dim = 6;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         coef[0] = new complex<double>(sqrt(3.0/5.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,1,0,2);

         coef[1] = new complex<double>(0.0,sqrt(3.0/5.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,1,2);

         coef[2] = new complex<double>(-sqrt(3.0)/4.0,0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,3,0,0);

         coef[3] = new complex<double>(0.0,-sqrt(3.0)/4.0);
         ind[3] = CartInt::ginlxyz2s(i,n,l,0,3,0);

         coef[4] = new complex<double>(-sqrt(3.0/5.0)*0.25,0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,1,2,0);

         coef[5] = new complex<double>(0.0,-sqrt(3.0/5.0)*0.25);
         ind[5] = CartInt::ginlxyz2s(i,n,l,2,1,0);

      }
      else if(m == 0){

         dim = 3;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         coef[0] = new complex<double>(1.0,0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,0,0,3);

         coef[1] = new complex<double>(-1.5/sqrt(5.0),0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,2,0,1);

         coef[2] = new complex<double>(-1.5/sqrt(5.0),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,0,2,1);

      }
      else if(m == -1){

         dim = 6;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         coef[0] = new complex<double>(sqrt(3.0/5.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,1,0,2);

         coef[1] = new complex<double>(0.0,-sqrt(3.0/5.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,1,2);

         coef[2] = new complex<double>(-sqrt(3.0)/4.0,0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,3,0,0);

         coef[3] = new complex<double>(0.0,sqrt(3.0)/4.0);
         ind[3] = CartInt::ginlxyz2s(i,n,l,0,3,0);

         coef[4] = new complex<double>(-sqrt(3.0/5.0)*0.25,0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,1,2,0);

         coef[5] = new complex<double>(0.0,sqrt(3.0/5.0)*0.25);
         ind[5] = CartInt::ginlxyz2s(i,n,l,2,1,0);

      }
      else if(m == -2){

         dim = 3;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         coef[0] = new complex<double>(sqrt(3.0/8.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,2,0,1);

         coef[1] = new complex<double>(-sqrt(3.0/8.0),0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,2,1);

         coef[2] = new complex<double>(0.0,-1.0/sqrt(2.0));
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,1,1);

      }
      else if(m == -3){

         dim = 4;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(5.0)/4.0,0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,3,0,0);

         coef[1] = new complex<double>(0.0,sqrt(5.0)/4.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,3,0);

         coef[2] = new complex<double>(-0.75,0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,2,0);

         coef[3] = new complex<double>(0.0,-0.75);
         ind[3] = CartInt::ginlxyz2s(i,n,l,2,1,0);

      }

   }
   else if(l == 4){

      if(m == 4){

         int dim = 5;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(35.0/2.0)/8.0,0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,4,0,0);

         coef[1] = new complex<double>(sqrt(35.0/2.0)/8.0,0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,4,0);

         coef[2] = new complex<double>(-0.75*sqrt(1.5),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,2,2,0);

         coef[3] = new complex<double>(0.0,sqrt(5.0/8.0));
         ind[3] = CartInt::ginlxyz2s(i,n,l,3,1,0);

         coef[4] = new complex<double>(0.0,-sqrt(5.0/8.0));
         ind[4] = CartInt::ginlxyz2s(i,n,l,1,3,0);

      }
      else if(m == 3){

         int dim = 4;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(5.0)/4.0,0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,3,0,1);

         coef[1] = new complex<double>(0.0,-sqrt(5.0)/4.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,3,1);

         coef[2] = new complex<double>(-0.75,0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,2,1);

         coef[3] = new complex<double>(0.0,0.75);
         ind[3] = CartInt::ginlxyz2s(i,n,l,2,1,1);

      }
      else if(m == 2){

         int dim = 7;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(1.5*sqrt(3.0/14),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,2,0,2);

         coef[1] = new complex<double>(-1.5*sqrt(3.0/14),0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,2,2);

         coef[2] = new complex<double>(0.0,3.0/sqrt(14));
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,1,2);

         coef[3] = new complex<double>(0.25*sqrt(2.5),0.0);
         ind[3] = CartInt::ginlxyz2s(i,n,l,4,0,0);

         coef[4] = new complex<double>(-0.25*sqrt(2.5),0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,0,4,0);

         coef[5] = new complex<double>(0.0,-0.5*sqrt(5.0/14.0));
         ind[5] = CartInt::ginlxyz2s(i,n,l,3,1,0);

         coef[6] = new complex<double>(0.0,-0.5*sqrt(5.0/14.0));
         ind[6] = CartInt::ginlxyz2s(i,n,l,1,3,0);

      }
      else if(m == 1){

         int dim = 6;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(5.0/7.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,1,0,3);

         coef[1] = new complex<double>(0.0,sqrt(5.0/7.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,1,3);

         coef[2] = new complex<double>(-0.75*sqrt(5.0/7.0),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,3,0,1);

         coef[3] = new complex<double>(0.0,-0.75*sqrt(5.0/7.0));
         ind[3] = CartInt::ginlxyz2s(i,n,l,0,3,1);

         coef[4] = new complex<double>(-0.75/sqrt(7.0),0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,1,2,1);

         coef[5] = new complex<double>(0.0,-0.75/sqrt(7.0));
         ind[5] = CartInt::ginlxyz2s(i,n,l,2,1,1);

      }
      else if(m == 0){

         int dim = 6;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(1.0,0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,0,0,4);

         coef[1] = new complex<double>(3.0/8.0,0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,4,0,0);

         coef[2] = new complex<double>(3.0/8.0,0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,0,4,0);

         coef[3] = new complex<double>(-3.0*sqrt(3.0/35.0),0.0);
         ind[3] = CartInt::ginlxyz2s(i,n,l,2,0,2);

         coef[4] = new complex<double>(-3.0*sqrt(3.0/35.0),0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,0,2,2);

         coef[5] = new complex<double>(0.75*sqrt(3.0/35.0),0.0);
         ind[5] = CartInt::ginlxyz2s(i,n,l,2,2,0);

      }
      else if(m == -1){

         int dim = 6;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(5.0/7.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,1,0,3);

         coef[1] = new complex<double>(0.0,-sqrt(5.0/7.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,1,3);

         coef[2] = new complex<double>(-0.75*sqrt(5.0/7.0),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,3,0,1);

         coef[3] = new complex<double>(0.0,0.75*sqrt(5.0/7.0));
         ind[3] = CartInt::ginlxyz2s(i,n,l,0,3,1);

         coef[4] = new complex<double>(-0.75/sqrt(7.0),0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,1,2,1);

         coef[5] = new complex<double>(0.0,0.75/sqrt(7.0));
         ind[5] = CartInt::ginlxyz2s(i,n,l,2,1,1);

      }
      else if(m == -2){

         int dim = 7;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(1.5*sqrt(3.0/14),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,2,0,2);

         coef[1] = new complex<double>(-1.5*sqrt(3.0/14),0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,2,2);

         coef[2] = new complex<double>(0.0,-3.0/sqrt(14));
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,1,2);

         coef[3] = new complex<double>(0.25*sqrt(2.5),0.0);
         ind[3] = CartInt::ginlxyz2s(i,n,l,4,0,0);

         coef[4] = new complex<double>(-0.25*sqrt(2.5),0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,0,4,0);

         coef[5] = new complex<double>(0.0,0.5*sqrt(5.0/14.0));
         ind[5] = CartInt::ginlxyz2s(i,n,l,3,1,0);

         coef[6] = new complex<double>(0.0,0.5*sqrt(5.0/14.0));
         ind[6] = CartInt::ginlxyz2s(i,n,l,1,3,0);

      }
      else if(m == -3){

         int dim = 4;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(5.0)/4.0,0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,3,0,1);

         coef[1] = new complex<double>(0.0,sqrt(5.0)/4.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,3,1);

         coef[2] = new complex<double>(-0.75,0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,2,1);

         coef[3] = new complex<double>(0.0,-0.75);
         ind[3] = CartInt::ginlxyz2s(i,n,l,2,1,1);


      }
      else if(m == -4){

         int dim = 5;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(35.0/2.0)/8.0,0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,4,0,0);

         coef[1] = new complex<double>(sqrt(35.0/2.0)/8.0,0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,4,0);

         coef[2] = new complex<double>(-0.75*sqrt(1.5),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,2,2,0);

         coef[3] = new complex<double>(0.0,-sqrt(5.0/8.0));
         ind[3] = CartInt::ginlxyz2s(i,n,l,3,1,0);

         coef[4] = new complex<double>(0.0,sqrt(5.0/8.0));
         ind[4] = CartInt::ginlxyz2s(i,n,l,1,3,0);

      }

   }
   else if(l == 5){

      if(m == 5){

         int dim = 6;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(3.0/16.0*sqrt(7.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,5,0,0);

         coef[1] = new complex<double>(0.0,3.0/16.0*sqrt(7.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,5,0);

         coef[2] = new complex<double>(5.0/16.0*sqrt(7.0),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,4,0);

         coef[3] = new complex<double>(0.0,5.0/16.0*sqrt(7.0));
         ind[3] = CartInt::ginlxyz2s(i,n,l,4,1,0);

         coef[4] = new complex<double>(-5.0/8.0*sqrt(3.0),0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,3,2,0);

         coef[5] = new complex<double>(0.0,-5.0/8.0*sqrt(3.0));
         ind[5] = CartInt::ginlxyz2s(i,n,l,2,3,0);

      }
      else if(m == 4){

         int dim = 5;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(17.5)/8.0,0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,4,0,1);

         coef[1] = new complex<double>(sqrt(17.5)/8.0,0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,4,1);

         coef[2] = new complex<double>(-0.75*sqrt(1.5),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,2,2,1);

         coef[3] = new complex<double>(0.0,0.5*sqrt(2.5));
         ind[3] = CartInt::ginlxyz2s(i,n,l,3,1,1);

         coef[4] = new complex<double>(0.0,-0.5*sqrt(2.5));
         ind[4] = CartInt::ginlxyz2s(i,n,l,1,3,1);

      }
      else if(m == 3){

         int dim = 10;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(0.5*sqrt(5.0/3.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,3,0,2);

         coef[1] = new complex<double>(0.0,-0.5*sqrt(5.0/3.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,3,2);

         coef[2] = new complex<double>(-0.5*sqrt(3.0),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,2,2);

         coef[3] = new complex<double>(0.0,0.5*sqrt(3.0));
         ind[3] = CartInt::ginlxyz2s(i,n,l,2,1,2);

         coef[4] = new complex<double>(-sqrt(35.0)/16.0,0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,5,0,0);

         coef[5] = new complex<double>(0.0,sqrt(35.0)/16.0);
         ind[5] = CartInt::ginlxyz2s(i,n,l,0,5,0);

         coef[6] = new complex<double>(sqrt(35.0)/16.0,0.0);
         ind[6] = CartInt::ginlxyz2s(i,n,l,1,4,0);

         coef[7] = new complex<double>(0.0,-sqrt(35.0)/16.0);
         ind[7] = CartInt::ginlxyz2s(i,n,l,4,1,0);

         coef[8] = new complex<double>(sqrt(5.0/3.0)/8.0,0.0);
         ind[8] = CartInt::ginlxyz2s(i,n,l,3,2,0);

         coef[9] = new complex<double>(0.0,-sqrt(5.0/3.0)/8.0);
         ind[9] = CartInt::ginlxyz2s(i,n,l,2,3,0);

      }
      else if(m == 2){

         int dim = 7;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(5.0/8.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,2,0,3);

         coef[1] = new complex<double>(-sqrt(5.0/8.0),0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,2,3);

         coef[2] = new complex<double>(0.0,sqrt(5.0/6.0));
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,1,3);

         coef[3] = new complex<double>(-0.25*sqrt(35.0/6.0),0.0);
         ind[3] = CartInt::ginlxyz2s(i,n,l,4,0,1);

         coef[4] = new complex<double>(0.25*sqrt(35.0/6.0),0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,0,4,1);

         coef[5] = new complex<double>(0.0,-0.5*sqrt(5.0/6.0));
         ind[5] = CartInt::ginlxyz2s(i,n,l,3,1,1);

         coef[6] = new complex<double>(0.0,-0.5*sqrt(5.0/6.0));
         ind[6] = CartInt::ginlxyz2s(i,n,l,1,3,1);

      }
      else if(m == 1){

         int dim = 12;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(5.0/6.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,1,0,4);

         coef[1] = new complex<double>(0.0,sqrt(5.0/6.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,1,4);

         coef[2] = new complex<double>(-1.5*sqrt(5.0/14.0),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,3,0,2);

         coef[3] = new complex<double>(0.0,-1.5*sqrt(5.0/14.0));
         ind[3] = CartInt::ginlxyz2s(i,n,l,0,3,2);

         coef[4] = new complex<double>(-1.5/sqrt(14.0),0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,1,2,2);

         coef[5] = new complex<double>(0.0,-1.5/sqrt(14.0));
         ind[5] = CartInt::ginlxyz2s(i,n,l,2,1,2);

         coef[6] = new complex<double>(sqrt(7.5)/8.0,0.0);
         ind[6] = CartInt::ginlxyz2s(i,n,l,5,0,0);

         coef[7] = new complex<double>(0.0,sqrt(7.5)/8.0);
         ind[7] = CartInt::ginlxyz2s(i,n,l,0,5,0);

         coef[8] = new complex<double>(sqrt(5.0/6.0)/8.0,0.0);
         ind[8] = CartInt::ginlxyz2s(i,n,l,1,4,0);

         coef[9] = new complex<double>(0.0,sqrt(5.0/6.0)/8.0);
         ind[9] = CartInt::ginlxyz2s(i,n,l,4,1,0);

         coef[10] = new complex<double>(0.25*sqrt(5.0/14.0),0.0);
         ind[10] = CartInt::ginlxyz2s(i,n,l,3,2,0);

         coef[11] = new complex<double>(0.0,0.25*sqrt(5.0/14.0));
         ind[11] = CartInt::ginlxyz2s(i,n,l,2,3,0);

      }
      else if(m == 0){

         int dim = 6;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(1.0,0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,0,0,5);

         coef[1] = new complex<double>(-5.0/sqrt(21.0),0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,2,0,3);

         coef[2] = new complex<double>(-5.0/sqrt(21.0),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,0,2,3);

         coef[3] = new complex<double>(5.0/8.0,0.0);
         ind[3] = CartInt::ginlxyz2s(i,n,l,4,0,1);

         coef[4] = new complex<double>(5.0/8.0,0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,0,4,1);

         coef[5] = new complex<double>(0.25*sqrt(15.0/7.0),0.0);
         ind[5] = CartInt::ginlxyz2s(i,n,l,2,2,1);

      }
      else if(m == -1){

         int dim = 12;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(5.0/6.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,1,0,4);

         coef[1] = new complex<double>(0.0,-sqrt(5.0/6.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,1,4);

         coef[2] = new complex<double>(-1.5*sqrt(5.0/14.0),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,3,0,2);

         coef[3] = new complex<double>(0.0,1.5*sqrt(5.0/14.0));
         ind[3] = CartInt::ginlxyz2s(i,n,l,0,3,2);

         coef[4] = new complex<double>(-1.5/sqrt(14.0),0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,1,2,2);

         coef[5] = new complex<double>(0.0,1.5/sqrt(14.0));
         ind[5] = CartInt::ginlxyz2s(i,n,l,2,1,2);

         coef[6] = new complex<double>(sqrt(7.5)/8.0,0.0);
         ind[6] = CartInt::ginlxyz2s(i,n,l,5,0,0);

         coef[7] = new complex<double>(0.0,-sqrt(7.5)/8.0);
         ind[7] = CartInt::ginlxyz2s(i,n,l,0,5,0);

         coef[8] = new complex<double>(sqrt(5.0/6.0)/8.0,0.0);
         ind[8] = CartInt::ginlxyz2s(i,n,l,1,4,0);

         coef[9] = new complex<double>(0.0,-sqrt(5.0/6.0)/8.0);
         ind[9] = CartInt::ginlxyz2s(i,n,l,4,1,0);

         coef[10] = new complex<double>(0.25*sqrt(5.0/14.0),0.0);
         ind[10] = CartInt::ginlxyz2s(i,n,l,3,2,0);

         coef[11] = new complex<double>(0.0,-0.25*sqrt(5.0/14.0));
         ind[11] = CartInt::ginlxyz2s(i,n,l,2,3,0);

      }
      else if(m == -2){

         int dim = 7;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(5.0/8.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,2,0,3);

         coef[1] = new complex<double>(-sqrt(5.0/8.0),0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,2,3);

         coef[2] = new complex<double>(0.0,-sqrt(5.0/6.0));
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,1,3);

         coef[3] = new complex<double>(-0.25*sqrt(35.0/6.0),0.0);
         ind[3] = CartInt::ginlxyz2s(i,n,l,4,0,1);

         coef[4] = new complex<double>(0.25*sqrt(35.0/6.0),0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,0,4,1);

         coef[5] = new complex<double>(0.0,0.5*sqrt(5.0/6.0));
         ind[5] = CartInt::ginlxyz2s(i,n,l,3,1,1);

         coef[6] = new complex<double>(0.0,0.5*sqrt(5.0/6.0));
         ind[6] = CartInt::ginlxyz2s(i,n,l,1,3,1);

      }
      else if(m == -3){

         int dim = 10;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(0.5*sqrt(5.0/3.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,3,0,2);

         coef[1] = new complex<double>(0.0,0.5*sqrt(5.0/3.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,3,2);

         coef[2] = new complex<double>(-0.5*sqrt(3.0),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,2,2);

         coef[3] = new complex<double>(0.0,-0.5*sqrt(3.0));
         ind[3] = CartInt::ginlxyz2s(i,n,l,2,1,2);

         coef[4] = new complex<double>(-sqrt(35.0)/16.0,0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,5,0,0);

         coef[5] = new complex<double>(0.0,-sqrt(35.0)/16.0);
         ind[5] = CartInt::ginlxyz2s(i,n,l,0,5,0);

         coef[6] = new complex<double>(sqrt(35.0)/16.0,0.0);
         ind[6] = CartInt::ginlxyz2s(i,n,l,1,4,0);

         coef[7] = new complex<double>(0.0,sqrt(35.0)/16.0);
         ind[7] = CartInt::ginlxyz2s(i,n,l,4,1,0);

         coef[8] = new complex<double>(sqrt(5.0/3.0)/8.0,0.0);
         ind[8] = CartInt::ginlxyz2s(i,n,l,3,2,0);

         coef[9] = new complex<double>(0.0,sqrt(5.0/3.0)/8.0);
         ind[9] = CartInt::ginlxyz2s(i,n,l,2,3,0);

      }
      else if(m == -4){

         int dim = 5;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(sqrt(17.5)/8.0,0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,4,0,1);

         coef[1] = new complex<double>(sqrt(17.5)/8.0,0.0);
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,4,1);

         coef[2] = new complex<double>(-0.75*sqrt(1.5),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,2,2,1);

         coef[3] = new complex<double>(0.0,-0.5*sqrt(2.5));
         ind[3] = CartInt::ginlxyz2s(i,n,l,3,1,1);

         coef[4] = new complex<double>(0.0,0.5*sqrt(2.5));
         ind[4] = CartInt::ginlxyz2s(i,n,l,1,3,1);

      }
      else if(m == -5){

         int dim = 6;

         coef = new complex<double> * [dim];
         ind = new int [dim];

         //start the transformation
         coef[0] = new complex<double>(3.0/16.0*sqrt(7.0),0.0);
         ind[0] = CartInt::ginlxyz2s(i,n,l,5,0,0);

         coef[1] = new complex<double>(0.0,-3.0/16.0*sqrt(7.0));
         ind[1] = CartInt::ginlxyz2s(i,n,l,0,5,0);

         coef[2] = new complex<double>(5.0/16.0*sqrt(7.0),0.0);
         ind[2] = CartInt::ginlxyz2s(i,n,l,1,4,0);

         coef[3] = new complex<double>(0.0,-5.0/16.0*sqrt(7.0));
         ind[3] = CartInt::ginlxyz2s(i,n,l,4,1,0);

         coef[4] = new complex<double>(-5.0/8.0*sqrt(3.0),0.0);
         ind[4] = CartInt::ginlxyz2s(i,n,l,3,2,0);

         coef[5] = new complex<double>(0.0,5.0/8.0*sqrt(3.0));
         ind[5] = CartInt::ginlxyz2s(i,n,l,2,3,0);

      }

   }
   else
      cout << "No more transformations" << endl;

}

/** 
 * copy constructor
 * @param tf_c Transform object to be copied in the newly constructed object
 */
Transform::Transform(const Transform &tf_c){ 

   dim = tf_c.gdim();

   coef = new complex<double> * [dim];
   ind = new int [dim];

   for(int i = 0;i < dim;++i){

      coef[i] = new complex<double>(tf_c.gcoef(i));
      ind[i] = tf_c.gind(i);

   }

}

/**
 * standard destructor
 */
Transform::~Transform(){

   for(int i = 0;i < dim;++i)
      delete coef[i];

   delete [] coef;

   delete [] ind;

}

/**
 * @return the dimension, i.e. number of terms in the transformation
 */
int Transform::gdim() const {

   return dim;

}

/**
 * @param i the "i'th" term in the transformation
 * @return the coefficient corresponding to the "i'th" term in the transformation
 */
complex<double> Transform::gcoef(int i) const {

   return *coef[i];

}

/**
 * @param i the "i'th" term in the transformation
 * @return the cartesian sp index corresponding to the "i'th" term in the transformation
 */
int Transform::gind(int i) const {

   return ind[i];

}
