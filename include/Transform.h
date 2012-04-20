#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <iostream>
#include <fstream>
#include <complex>

using std::ostream;
using std::complex;

/**
 * @author Brecht Verstichel
 * @date 16-04-2012\n\n
 * This class generates the transformation between spherical and cartesian gaussian wavefunctions
 */
class Transform {

   /**
    * Output stream operator overloaded
    * @param tf the Transform you want to print
    */ friend ostream &operator<<(ostream &output,Transform &tf);

   public:
      
      //constructor
      Transform(int,int,int,int);

      //copy constructor
      Transform(const Transform &);

      //destructor
      virtual ~Transform();

      int gdim() const;

      complex<double> gcoef(int) const;

      int gind(int) const;

   private:

      //!spherical indices
      int i,n,l,m;
      
      //!array of complex doubles containing the coefficients of the transformation
      complex<double> **coef;

      //!array of ints containing the indices of the tranformation
      int *ind;

      //!dimension of the transformation, i.e. the number of different terms in the transformation
      int dim;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
