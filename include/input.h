#ifndef INPUT_H
#define INPUT_H

#include "preamble.h"
#include "Gauss.h"
#include "R.h"

/**
 * @author Sebastian Wouters, Brecht Verstichel
 * @date 11-06-2012\n\n
 * This is a class that handles the ugly input of the program. Reads in the necessary files for further processing.
 */
class input{

   public:
      
      static void init(string);

      static void clear();

      static int gcharge();

      static int gN_Z();

      static string gbasisset();

      static int gZ(int);

      static R &gR_nc(int);

      static const R &gR(int);

      static const Gauss &gGaussInfo(int);

      static int NumberOfElectrons();

      static double gNucRepEn();

   private:

      //!To know the proton number of an element: elements[Z-1] = "name".
      static string *elements;

      //!charge = number of protons - number of electrons
      static int charge;
      
      //!Basisset
      static string basisset;

      //!Number of cores
      static int N_Z;

      //!Cores (atomic number)
      static int *Z;

      //!Positions of the cores
      static R **r;

      //!Per core: a Gauss object that stores the basic info to construct the desired Gaussian contractions
      static Gauss **GaussInfo;

      //!nr of elements?
      static int MENDELJEV;

      //!the nuclear repuslion energy
      static double NucRepEn;

      //Helper functions
      static void initelements();
      static void readinsetupfile(string);
      static void readinsetupfile(istream&);
      static void fillgaussinfo();
      static int getZ(string);

};

#endif


/* vim: set ts=3 sw=3 expandtab :*/
