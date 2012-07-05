#ifndef LIBINT_H
#define LIBINT_H

#include <iostream>
#include <cstdlib>
#include <libint2.h>

class LibInt{

   public:
      
      static void init();

      static void clear();

      static int gMAXFAC();

      static int gEPS();

      static int gL_MAX();

      static double gF(int);

      static double gf(int,int,int,int,int,int);

      static void calc_f(int,double);

      template<typename LibintEval>
         static void prep_libint2(LibintEval* erieval, unsigned int am1, double alpha1,
               double A[3], unsigned int am2, double alpha2, double B[3],
               unsigned int am3, double alpha3, double C[3],
               unsigned int am4, double alpha4, double D[3], int norm_flag) {

            const unsigned int am = am1 + am2 + am3 + am4;

            const double gammap = alpha1 + alpha2;

            const double Px = (alpha1 * A[0] + alpha2 * B[0]) / gammap;
            const double Py = (alpha1 * A[1] + alpha2 * B[1]) / gammap;
            const double Pz = (alpha1 * A[2] + alpha2 * B[2]) / gammap;
            const double PAx = Px - A[0];
            const double PAy = Py - A[1];
            const double PAz = Pz - A[2];
            const double PBx = Px - B[0];
            const double PBy = Py - B[1];
            const double PBz = Pz - B[2];
            const double AB_x = A[0] - B[0];
            const double AB_y = A[1] - B[1];
            const double AB_z = A[2] - B[2];
            const double AB2 = AB_x * AB_x + AB_y * AB_y + AB_z * AB_z;

#if LIBINT2_DEFINED(eri,PA_x)
            erieval->PA_x[0] = PAx;
#endif
#if LIBINT2_DEFINED(eri,PA_y)
            erieval->PA_y[0] = PAy;
#endif
#if LIBINT2_DEFINED(eri,PA_z)
            erieval->PA_z[0] = PAz;
#endif
#if LIBINT2_DEFINED(eri,PB_x)
            erieval->PB_x[0] = PBx;
#endif
#if LIBINT2_DEFINED(eri,PB_y)
            erieval->PB_y[0] = PBy;
#endif
#if LIBINT2_DEFINED(eri,PB_z)
            erieval->PB_z[0] = PBz;
#endif

#if LIBINT2_DEFINED(eri,AB_x)
            erieval->AB_x[0] = AB_x;
#endif
#if LIBINT2_DEFINED(eri,AB_y)
            erieval->AB_y[0] = AB_y;
#endif
#if LIBINT2_DEFINED(eri,AB_z)
            erieval->AB_z[0] = AB_z;
#endif
#if LIBINT2_DEFINED(eri,oo2z)
            erieval->oo2z[0] = 0.5/gammap;
#endif

            const double gammaq = alpha3 + alpha4;
            const double gammapq = gammap * gammaq / (gammap + gammaq);
            const double Qx = (alpha3 * C[0] + alpha4 * D[0]) / gammaq;
            const double Qy = (alpha3 * C[1] + alpha4 * D[1]) / gammaq;
            const double Qz = (alpha3 * C[2] + alpha4 * D[2]) / gammaq;
            const double QCx = Qx - C[0];
            const double QCy = Qy - C[1];
            const double QCz = Qz - C[2];
            const double QDx = Qx - D[0];
            const double QDy = Qy - D[1];
            const double QDz = Qz - D[2];
            const double CD_x = C[0] - D[0];
            const double CD_y = C[1] - D[1];
            const double CD_z = C[2] - D[2];
            const double CD2 = CD_x * CD_x + CD_y * CD_y + CD_z * CD_z;

#if LIBINT2_DEFINED(eri,QC_x)
            erieval->QC_x[0] = QCx;
#endif
#if LIBINT2_DEFINED(eri,QC_y)
            erieval->QC_y[0] = QCy;
#endif
#if LIBINT2_DEFINED(eri,QC_z)
            erieval->QC_z[0] = QCz;
#endif
#if LIBINT2_DEFINED(eri,QD_x)
            erieval->QD_x[0] = QDx;
#endif
#if LIBINT2_DEFINED(eri,QD_y)
            erieval->QD_y[0] = QDy;
#endif
#if LIBINT2_DEFINED(eri,QD_z)
            erieval->QD_z[0] = QDz;
#endif

#if LIBINT2_DEFINED(eri,CD_x)
            erieval->CD_x[0] = CD_x;
#endif
#if LIBINT2_DEFINED(eri,CD_y)
            erieval->CD_y[0] = CD_y;
#endif
#if LIBINT2_DEFINED(eri,CD_z)
            erieval->CD_z[0] = CD_z;
#endif
#if LIBINT2_DEFINED(eri,oo2e)
            erieval->oo2e[0] = 0.5/gammaq;
#endif

            // Prefactors for interelectron transfer relation
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_x)
            erieval->TwoPRepITR_pfac0_0_x[0] = - (alpha2*AB_x + alpha4*CD_x)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_y)
            erieval->TwoPRepITR_pfac0_0_y[0] = - (alpha2*AB_y + alpha4*CD_y)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_z)
            erieval->TwoPRepITR_pfac0_0_z[0] = - (alpha2*AB_z + alpha4*CD_z)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_x)
            erieval->TwoPRepITR_pfac0_1_x[0] = - (alpha2*AB_x + alpha4*CD_x)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_y)
            erieval->TwoPRepITR_pfac0_1_y[0] = - (alpha2*AB_y + alpha4*CD_y)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_z)
            erieval->TwoPRepITR_pfac0_1_z[0] = - (alpha2*AB_z + alpha4*CD_z)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac1_0)
            erieval->TwoPRepITR_pfac1_0[0] = -gammaq/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac1_1)
            erieval->TwoPRepITR_pfac1_1[0] = -gammap/gammaq;
#endif

            const double PQx = Px - Qx;
            const double PQy = Py - Qy;
            const double PQz = Pz - Qz;
            const double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;
            const double Wx = (gammap * Px + gammaq * Qx) / (gammap + gammaq);
            const double Wy = (gammap * Py + gammaq * Qy) / (gammap + gammaq);
            const double Wz = (gammap * Pz + gammaq * Qz) / (gammap + gammaq);

#if LIBINT2_DEFINED(eri,WP_x)
            erieval->WP_x[0] = Wx - Px;
#endif
#if LIBINT2_DEFINED(eri,WP_y)
            erieval->WP_y[0] = Wy - Py;
#endif
#if LIBINT2_DEFINED(eri,WP_z)
            erieval->WP_z[0] = Wz - Pz;
#endif
#if LIBINT2_DEFINED(eri,WQ_x)
            erieval->WQ_x[0] = Wx - Qx;
#endif
#if LIBINT2_DEFINED(eri,WQ_y)
            erieval->WQ_y[0] = Wy - Qy;
#endif
#if LIBINT2_DEFINED(eri,WQ_z)
            erieval->WQ_z[0] = Wz - Qz;
#endif
#if LIBINT2_DEFINED(eri,oo2ze)
            erieval->oo2ze[0] = 0.5/(gammap+gammaq);
#endif
#if LIBINT2_DEFINED(eri,roz)
            erieval->roz[0] = gammapq/gammap;
#endif
#if LIBINT2_DEFINED(eri,roe)
            erieval->roe[0] = gammapq/gammaq;
#endif

            double K1 = exp(-alpha1 * alpha2 * AB2 / gammap);
            double K2 = exp(-alpha3 * alpha4 * CD2 / gammaq);
            double pfac = 2 * pow(M_PI, 2.5) * K1 * K2 / (gammap * gammaq * sqrt(gammap + gammaq));

            LibInt::calc_f(am, PQ2 * gammapq);

            // using dangerous macros from libint2.h
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
            erieval->LIBINT_T_SS_EREP_SS(0)[0] = pfac*LibInt::gF(0);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
            erieval->LIBINT_T_SS_EREP_SS(1)[0] = pfac*LibInt::gF(1);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
            erieval->LIBINT_T_SS_EREP_SS(2)[0] = pfac*LibInt::gF(2);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
            erieval->LIBINT_T_SS_EREP_SS(3)[0] = pfac*LibInt::gF(3);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
            erieval->LIBINT_T_SS_EREP_SS(4)[0] = pfac*LibInt::gF(4);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
            erieval->LIBINT_T_SS_EREP_SS(5)[0] = pfac*LibInt::gF(5);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
            erieval->LIBINT_T_SS_EREP_SS(6)[0] = pfac*LibInt::gF(6);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
            erieval->LIBINT_T_SS_EREP_SS(7)[0] = pfac*LibInt::gF(7);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
            erieval->LIBINT_T_SS_EREP_SS(8)[0] = pfac*LibInt::gF(8);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
            erieval->LIBINT_T_SS_EREP_SS(9)[0] = pfac*LibInt::gF(9);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
            erieval->LIBINT_T_SS_EREP_SS(10)[0] = pfac*LibInt::gF(10);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
            erieval->LIBINT_T_SS_EREP_SS(11)[0] = pfac*LibInt::gF(11);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
            erieval->LIBINT_T_SS_EREP_SS(12)[0] = pfac*LibInt::gF(12);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
            erieval->LIBINT_T_SS_EREP_SS(13)[0] = pfac*LibInt::gF(13);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
            erieval->LIBINT_T_SS_EREP_SS(14)[0] = pfac*LibInt::gF(14);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
            erieval->LIBINT_T_SS_EREP_SS(15)[0] = pfac*LibInt::gF(15);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
            erieval->LIBINT_T_SS_EREP_SS(16)[0] = pfac*LibInt::gF(16);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
            erieval->LIBINT_T_SS_EREP_SS(17)[0] = pfac*LibInt::gF(17);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
            erieval->LIBINT_T_SS_EREP_SS(18)[0] = pfac*LibInt::gF(18);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
            erieval->LIBINT_T_SS_EREP_SS(19)[0] = pfac*LibInt::gF(19);
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
            erieval->LIBINT_T_SS_EREP_SS(20)[0] = pfac*LibInt::gF(20);
#endif

         }


   private:

      //!array needed for the calculation of the Boys function.
      static double *F;

      //!array containing the f functions needed for the overlap calculation
      static double *****f;

      //!definition of the max angular momentum value allowed by the program
      static const int L_MAX;

      //!constant defining the size of the array containing the double faculties
      static const int MAXFAC;

      //!constant defining the precision with which to calculate the boys function F
      static const double EPS;

};

#endif
