#ifndef CI_SPPM_m_H
#define CI_SPPM_m_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"
#include "input.h"

/**
 * @author Brecht Verstichel
 * @date 22-04-2012\n\n
 * This class is specifically written for the storage of the nuclear attraction matrix elements of contracted gaussians.
 * It's a tensor of rank three, just a basic CI_SPPM matrix with a layers of (m) added on top, depending on how many
 * orders I need to calculate the full (0) order nuclear attraction
 */
class CI_SPPM_m {

   /**
    * Output stream operator overloaded
    * @param output the output stream e.g. cout
    * @param ci_p the CI_SPPM_m you want to print
    */
   friend ostream &operator<<(ostream &output,CI_SPPM_m &ci_p);

   public:
      
      //constructor
      CI_SPPM_m();

      //copy constructor
      CI_SPPM_m(const CI_SPPM_m &);

      double operator()(int,int,int,int,int,int,int,int,int,int,int,int,int) const;

      double operator()(int,int,int) const;

      double &operator()(int,int,int);

      void constr_nucl_attrac(int core);

      //destructor
      virtual ~CI_SPPM_m();

      static void init();

      static void clear();

   private:

      //!matrix containing the nuclear attraction elements
      double ***U;

      //!list returning the heigth of the stack on element (a,b)
      static int **ab2m;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
