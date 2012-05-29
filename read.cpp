/**
 * @author Ward Poelmans, Brecht Verstichel
 * @date 09-05-2012
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

//includes all important headers and defines which conditions are
//going to be used:
#include "include.h"

/**
 * This program reads in a TPM object from a HDF5 file.
 * It first gets the setupdata needed for the init's and
 * then get the actual TPM data.
 */
int main(int argc, char **argv)
{
   cout.precision(10);

   char *filename = 0;

   struct option long_options[] =
   {
      {"file",  required_argument, 0, 'f'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;
   while( (j = getopt_long (argc, argv, "hf:", long_options, &i)) != -1)
      switch(j)
      {
         case 'h':
         case '?':
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
               "\n"
               "    -f, --file=file.h5           file to read (in HDF5 format)\n"
               "    -h, --help                   Display this help\n"
               "\n";
            return 0;
            break;
         case 'f':
            filename = optarg;
            if( strlen(optarg) <= 0 )
            {
               std::cerr << "Couldn't find a filename. Please specify the file to read." << endl;
               return -1;
            }
            break;
      }

   if(!filename)
   {
      cout << "Usage: " << argv[0] << " [OPTIONS]\n"
         "\n"
         "    -f, --file=file.h5           file to read (in HDF5 format)\n"
         "    -h, --help                   Display this help\n"
         "\n";
      return 0;
   }

   cout << "Reading: " << filename << endl;

   string setupdata;

   TPM::ReadInitfromFile(filename,setupdata);

   CartInt::init(setupdata,false);
   SphInt::init();

   int M = 2*SphInt::gdim();//dim sp hilbert space
   int N = SphInt::gN();//nr of particles

   Tools::init(M,N);

   SPM::init(M,N);
   TPM::init(M,N);
   PHM::init(M,N);
   DPM::init(M,N);
   PPHM::init(M,N);

   SUP::init(M,N);
   EIG::init(M,N);

   CartInt ci;
   ci.norm();

   SphInt si(ci);
   si.orthogonalize();

   LinCon::init(M,N);
   LinIneq::init(M,N,si);

   //hamiltoniaan
   TPM ham;
   ham.molecule(si);

   TPM rdm;
   rdm.ReadfromFile(filename);

   cout << "Energy = " << rdm.ddot(ham) << endl;

   LinIneq::clear();

   PPHM::clear();
   DPM::clear();
   PHM::clear();
   TPM::clear();
   SPM::clear();

   Tools::clear();

   SphInt::clear();
   CartInt::clear();

   return 0;
}

/* vim: set ts=3 sw=3 expandtab :*/
