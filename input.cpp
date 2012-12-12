#include "preamble.h"
#include "input.h"
#include "Gauss.h"
#include "R.h"

int input::MENDELJEV = 118;

string *input::elements;

int input::charge;

string input::basisset;

int input::N_Z;

int *input::Z;

R **input::r;

Gauss **input::GaussInfo;

double input::NucRepEn;

/**
 * initializes the static objects, reads in the file
 * @param setupfile filename of the setupfile of the job
 */
void input::init(string setupfile){

   initelements();

   readinsetupfile(setupfile);

   fillgaussinfo();

   //set the Nuclear repulsion energy
   NucRepEn = 0.0;

   for (int i = 0;i < N_Z;i++)
      for (int j = i + 1;j < N_Z;j++)
         NucRepEn += Z[i] * Z[j] / sqrt(r[i]->dist_sqrd(*r[j]));

}

/**
 * deallocate the memory
 */
void input::clear(){

   delete [] elements;
   delete [] Z;

   for (int cnt=0; cnt < N_Z; cnt++){

      delete r[cnt];
      delete GaussInfo[cnt];

   }

   delete [] r;
   delete [] GaussInfo;

}


/**
 * Getter function for the global charge
 * @return the charge
 */
int input::gcharge(){

   return charge;

}


/**
 * Getter function for the number of Z in the problem
 * @return the number of cores in the problem
 */
int input::gN_Z(){

   return N_Z;

}


/**
 * Getter function for the basis set name
 * @return the basis set name
 */
string input::gbasisset(){

   return basisset;

}


/**
 * Getter function for charge of the i-th core
 * @param i the number of the core
 * @return the charge of the i-th core
 */
int input::gZ(int i){

   return Z[i];

}
/**
 * Getter function for the position vector of the i-th core
 * @param i the number of the core
 * @return the position vector of the i-th core
 */
R &input::gR_nc(int i){

   return *r[i];

}

/**
 * Getter function for the position vector of the i-th core
 * @param i the number of the core
 * @return the position vector of the i-th core
 */
const R &input::gR(int i){

   return *r[i];

}

/**
 * Getter function for the Gauss information object of the i-th core
 * @param cnt the number of the core
 * @return the Gauss information object of the i-th core
 */
const Gauss &input::gGaussInfo(int cnt){

   return *GaussInfo[cnt];

}

/**
 * Function to read in the elements mapper : elements[Z-1] = "name"
 */
void input::initelements(){

   elements = new string[MENDELJEV];
   string temp;
   int counter = 0;

   ifstream elem("basissets/elements.bs", ios::in);

   while(getline(elem,temp)) {

      int place = temp.find("\t") + 1;
      elements[counter] = temp.substr(place,temp.size()-place);
      counter++;

   }

   elem.close();

}


/**
 * Function that establishes the reverse elements mapper : getZ("name")=Z
 */
int input::getZ(string elem){

   for (int cnt=0; cnt<MENDELJEV; cnt++)
      if (!elem.compare(elements[cnt]))
         return cnt+1;

   return -1;

}

int input::NumberOfElectrons(){

   int count = 0;

   for (int i=0; i< N_Z; i++)
      count += Z[i];

   count -= charge;

   return count;

}

/*
 * Function to read in the setup file
 * @param setupfile filename of the setup file
 */
void input::readinsetupfile(string setupfile)
{
   ifstream inputfile(setupfile.c_str());

   readinsetupfile(inputfile);

   inputfile.close();
}

/**
 * Function to read in the setup file
 * @param inputfile stream which containts setup data
 */
void input::readinsetupfile(istream &inputfile){

   string temp, temp2, temp3;
   int pos1, pos2, pos3, length;

   //1: the basisset
   getline(inputfile,temp);
   pos1 = temp.find("=")+1;
   length = temp.size()-pos1;
   basisset = temp.substr(pos1,length);

   //2: the charge
   getline(inputfile,temp);
   pos1 = temp.find("=")+1;
   length = temp.size()-pos1;
   charge = atoi((temp.substr(pos1,length)).c_str());

   //4: white line
   getline(inputfile,temp);

   //5: predefined distances
   getline(inputfile,temp);
   pos1 = temp.find("=")+1;
   length = temp.size()-pos1;
   int numberofpred = atoi((temp.substr(pos1,length)).c_str());
   string * variables = new string[numberofpred];
   double * values = new double[numberofpred];

   //6: fill the arrays
   for (int cnt=0; cnt<numberofpred; cnt++){

      getline(inputfile,temp);
      pos1 = temp.find("=");
      variables[cnt] = temp.substr(0,pos1);
      pos1++;
      length = temp.size()-pos1;
      values[cnt] = atof((temp.substr(pos1,length)).c_str());

   }

   //7: white line
   getline(inputfile,temp);

   //8: number of cores;
   getline(inputfile,temp);
   pos1 = temp.find("=")+1;
   length = temp.size()-pos1;
   N_Z = atoi((temp.substr(pos1,length)).c_str());
   Z = new int[N_Z];
   r = new R * [N_Z];

   double * cotemp = new double[3];

   //9: read in cores
   for (int cnt=0; cnt<N_Z; cnt++){
      getline(inputfile,temp);

      //9.1 Core name -> Z
      pos1 = temp.find(" ")+1;
      pos2 = temp.find(" ",pos1);
      temp2 = temp.substr(pos1,pos2-pos1);
      Z[cnt] = getZ(temp2);

      //9.2 Core coordinates -> R
      for (int cnt2=0; cnt2<3; cnt2++){
         pos1 = pos2+1;
         if (cnt2<2) 
            pos2 = temp.find(" ",pos1);
         else
            pos2 = temp.size();
         temp2 = temp.substr(pos1,pos2-pos1);

         //Predefined distance used?
         pos3 = temp2.find("*");
         if (pos3==(int)string::npos){
            cotemp[cnt2] = atof(temp2.c_str());
         } else {
            temp3 = temp2.substr(pos3+1,temp2.size()-pos3-1);
            bool found = false;
            int counter = 0;
            while ((!found)&&(counter<numberofpred)){
               if (!temp3.compare(variables[counter]))
                  found = true;
               counter++;
            }
            cotemp[cnt2] = values[counter-1] * atof((temp2.substr(0,pos3)).c_str());
         }
      }

      r[cnt] = new R(cotemp[0],cotemp[1],cotemp[2]);
   }

   delete [] cotemp;
   delete [] variables;
   delete [] values;
}


/**
 * Function to create the GaussInfo array and fill it.
 */
void input::fillgaussinfo(){

   GaussInfo = new Gauss*[N_Z];

   string marge = basisset;
   while (marge.find("*")!=string::npos){
      marge.replace(marge.find("*"),1,"star");
   }

   marge = "basissets/" + marge + ".bs";
   ifstream basiss(marge.c_str(), ios::in);
   int maggy = 0;
   int numberofgauss, Z_i, cnt, count;

   bool itchy = true;

   //jump to start of the input of the first atom
   do {
      getline(basiss,marge);
      if (marge.find("***")!=string::npos)
         itchy = false;
   } while (itchy);

   bool scratchy = false;

   //loop over the atoms in the basisset and stop if (i) maggy equals N_Z or (ii) the end of the file has been reached
   do {

      //End of file? Scratchy kills the loop.
      if (!getline(basiss,marge)) scratchy = true;
      else {
         //which atom?
         Z_i = getZ(marge.substr(0,marge.find(" ")));

         //atom required? if yes: cnt contains the index of the core in our problem
         cnt = 0;
         itchy = true;
         while ((itchy)&&(cnt<N_Z)) {
            if (Z_i==Z[cnt]) itchy = false;
            else cnt++;
         }

         //if required: readin
         if (!itchy){

            //dump everything in an ugly vector and count i.t.m. the number of different contractions
            vector<string> milhouse;
            itchy = true;
            count = 0;
            do {
               getline(basiss,marge);
               if (marge.find("***")!=string::npos) itchy = false;
               else {
                  milhouse.push_back(marge);
                  if ((marge.substr(0,2)).compare("SP")==0) count += 2;
                  else
                     if ((marge.substr(0,1)).compare(" ")) count++;
               }
            } while (itchy);

            //make a Gauss object to store the contractions
            GaussInfo[cnt] = new Gauss(count);
            maggy++;

            //loop over the vector and store everything in the Gauss object
            count = 0;
            for (int cnt2=0; cnt2<(int)milhouse.size(); cnt2++){
               marge = milhouse[cnt2];
               if ((marge.substr(0,2)).compare("SP")==0) {
                  numberofgauss = atoi((marge.substr(4,3)).c_str());
                  double * alphas = new double[numberofgauss];
                  double * fors = new double[numberofgauss];
                  double * forp = new double[numberofgauss];

                  for (int cnt3=0; cnt3<numberofgauss; cnt3++){
                     marge = milhouse[cnt2+1+cnt3];

                     string n1 = marge.substr(0,27);
                     if (n1.find("D")!=string::npos) n1.replace(n1.find("D"),1,"E");
                     alphas[cnt3] = atof(n1.c_str());

                     n1 = marge.substr(27,23);
                     if (n1.find("D")!=string::npos) n1.replace(n1.find("D"),1,"E");
                     fors[cnt3] = atof(n1.c_str());

                     n1 = marge.substr(50,marge.size()-50);
                     if (n1.find("D")!=string::npos) n1.replace(n1.find("D"),1,"E");
                     forp[cnt3] = atof(n1.c_str());
                  }

                  GaussInfo[cnt]->set(count,numberofgauss,alphas,fors,'S');
                  GaussInfo[cnt]->set(count+1,numberofgauss,alphas,forp,'P');
                  count+=2;

                  cnt2+=numberofgauss;

                  delete [] alphas;
                  delete [] fors;
                  delete [] forp;

               } else {
                  numberofgauss = atoi((marge.substr(4,3)).c_str());
                  double * alphas = new double[numberofgauss];
                  double * coeff = new double[numberofgauss];

                  for (int cnt3=0; cnt3<numberofgauss; cnt3++){
                     marge = milhouse[cnt2+1+cnt3];

                     string n1 = marge.substr(0,27);
                     if (n1.find("D")!=string::npos) n1.replace(n1.find("D"),1,"E");
                     alphas[cnt3] = atof(n1.c_str());

                     n1 = marge.substr(27,marge.size()-50);
                     if (n1.find("D")!=string::npos) n1.replace(n1.find("D"),1,"E");
                     coeff[cnt3] = atof(n1.c_str());
                  }

                  GaussInfo[cnt]->set(count,numberofgauss,alphas,coeff,milhouse[cnt2].at(0));
                  count++;

                  cnt2+=numberofgauss;

                  delete [] alphas;
                  delete [] coeff;
               }

            }

            //check if other cores in the problem have the same number of protons. if yes: copy the gauss object.
            count = cnt;
            cnt++;
            while (cnt<N_Z) {
               if (Z_i ==Z[cnt]){
                  GaussInfo[cnt] = new Gauss(*(GaussInfo[count]));
                  maggy++;
               }
               cnt++;
            }

            //if the atom is not required move on to the next atom.
         } else {
            itchy = true;
            do {
               if (marge.find("***")!=string::npos) itchy = false;
            } while ((itchy)&&(getline(basiss,marge)));
         }
      }

   } while ((maggy<N_Z)&&(!scratchy));

   //terminate the program if the basisset doesn't contain all the required atoms
   if (maggy<N_Z){
      cout << "Basisset doesn't contain all the required atoms!" << endl;
      assert(maggy>=N_Z);
   }

   basiss.close();

}

/**
 * @return the nuclear-nuclear repulsion energy of a problem
 */
double input::gNucRepEn(){

   return NucRepEn;

}

/* vim: set ts=3 sw=3 expandtab :*/
