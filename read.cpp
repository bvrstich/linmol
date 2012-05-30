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
   char *setupfile = 0;

   struct option long_options[] =
   {
      {"file",  required_argument, 0, 'f'},
      {"setup",  required_argument, 0, 's'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;
   while( (j = getopt_long (argc, argv, "hf:s:", long_options, &i)) != -1)
      switch(j)
      {
         case 'h':
         case '?':
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
               "\n"
               "    -f, --file=file.h5           file to read (in HDF5 format)\n"
               "    -s  --setup=file.h5          use file to read in the hamiltonian matrix elements\n"
               "    -h, --help                   Display this help\n"
               "\n"
               "By default, the program rebuilds the hamiltonian matrix unless -s is given. The file\n"
               "and setupfile can be the same.\n"
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
         case 's':
            setupfile = optarg;
            if( strlen(optarg) <= 0 )
            {
               std::cerr << "Couldn't find a setupfile. Please specify the file to read." << endl;
               return -1;
            }
            break;
      }

   if(!filename)
   {
      cout << "Usage: " << argv[0] << " [OPTIONS]\n"
         "\n"
         "    -f, --file=file.h5           file to read (in HDF5 format)\n"
         "    -s  --setup=file.h5          use file to read in the hamiltonian matrix elements\n"
         "    -h, --help                   Display this help\n"
         "\n"
         "By default, the program rebuilds the hamiltonian matrix unless -s is given. The file\n"
         "and setupfile can be the same.\n"
         "\n";
      return 0;
   }

   cout << "Reading: " << filename << endl;

   if(setupfile)
      CartInt::initfromfile(setupfile);
   else
   {
      string setupdata;

      TPM::ReadInitfromFile(filename,setupdata);

      CartInt::init(setupdata,false);
   }

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

   CartInt *ci;

   if(setupfile)
      ci = new CartInt(setupfile);
   else
      ci = new CartInt;

   ci->norm();

   SphInt si(*ci);
   si.orthogonalize();

   LinCon::init(M,N);
   LinIneq::init(M,N,si);

   //hamiltoniaan
   TPM ham;
   ham.molecule(si);

   TPM rdm;
   rdm.ReadfromFile(filename);

   cout << "Energy = " << rdm.ddot(ham) + CartInt::gNucRepEn() << endl;

   LinIneq::clear();

   PPHM::clear();
   DPM::clear();
   PHM::clear();
   TPM::clear();
   SPM::clear();

   Tools::clear();

   SphInt::clear();
   CartInt::clear();

   delete ci;

   return 0;
}

/* vim: set ts=3 sw=3 expandtab :*/
