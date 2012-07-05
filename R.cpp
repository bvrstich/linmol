#include "preamble.h"
#include "R.h"

/**
 * emptyp constructor
 */
R::R() : Vector(3) { }

/** 
 * Constructor for the R class
 * @param xco the x coordinate
 * @param yco the y coordinate
 * @param zco the z coordinate
 */
R::R(double xco, double yco, double zco) : Vector(3) {

   (*this)[0] = xco;
   (*this)[1] = yco;
   (*this)[2] = zco;

}

/** 
 * Copy contructor for the R class
 * @param r_c the R to be copied
 */
R::R(const R& r_c) : Vector(r_c){ }


/**
 * Standard destructor
 */
R::~R(){ }


ostream &operator<<(ostream &output,const R &r_p){

   output << "[ " << r_p[0] << " ; " << r_p[1] << " ; " << r_p[2] << " ]";
   return output;

}

/**
 * Setter function for the coordinates
 * @param xco x coordinate
 * @param yco y coordinate
 * @param zco z coordinate
 */
void R::set(double xco, double yco, double zco){

   (*this)[0] = xco;
   (*this)[1] = yco;
   (*this)[2] = zco;

}

/* vim: set ts=3 sw=3 expandtab :*/
