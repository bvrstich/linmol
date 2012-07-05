#ifndef R_H
#define R_H

#include "preamble.h"
#include "Vector.h"

/**
 * @author Brecht Verstichel
 * @date 15-06-2012\n\n
 * This is a class written for real space vectors, it enherits from the Vector class, it's special
 * because it only has three entries, and the set function makes use of that, and it's too much of 
 * a hassle to remove it from Sebastian's original ThING.
 */
class R : public Vector {

   /**
    * Ostream overloader
    * @param output the ostream to write the info to
    * @param r_p the R of which the coordinates need to be written
    */
   friend ostream &operator<<(ostream &output,const R &r_p);

   public:

      //Constructor
      R();

      R(double, double, double);

      //Copy constructor
      R(const R &);

      //Destructor
      virtual ~R();

      //Setter
      void set(double, double, double);

   private:


};

#endif


/* vim: set ts=3 sw=3 expandtab :*/
