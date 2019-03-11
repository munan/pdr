/* Date: May 2, 2015, Author: Munan Gong
 * Shielding functions*/
#ifndef SHIELDING_H_
#define SHIELDING_H_

#include <math.h> /*a^x = pow(a,x)*/

class Shielding {
  public:
		/*van Dishoeck & Black's shielding function.
		 * NCO, NH2: column densith of CO and H2 in cm^-2.*/
    Shielding();
    ~Shielding();
		static double fShield_CO_vDB(const double NCO, const double NH2);
    /*CO sheilding from Visser+2009, Table 5*/
		static double fShield_CO_V09(const double NCO, const double NH2);
    /*H2 self shielding from Draine+Bertoldi1996*/
		static double fShield_H2(const double NH2, const double bH2);
    /*CI self shielding.*/
		static double fShield_C(const double NC, const double NH2);
    /*CI shielding of CO*/
		static double fShield_CO_C(const double NC);
  private:
		/*CO column density for DB table*/
		static const int len_NCO_DB_ = 8;
		static const double logNCOvDB_[len_NCO_DB_];
		/* H2 column densities for DB table*/
		static const int len_NH2_DB_ = 6;
		static const double logNH2vDB_[len_NH2_DB_];
		/* Tabulated shielding factors */
		static const double ThetavDB_[len_NH2_DB_][len_NCO_DB_];

    /*Visser+ 2009 Table 5, b(CO)=0.3, Tex(CO)=5K*/
    static const int len_NCO_V09_ = 47;
    static const double logNCOV09_[len_NCO_V09_];
    static const int len_NH2_V09_ = 42;
    static const double logNH2V09_[len_NH2_V09_];
		static const double ThetaV09_[len_NH2_V09_][len_NCO_V09_];

		/* Find index of linear interpretation, return the first index i for i,
		 * i+1.*/
		static int LinearInterpIndex(const int len, const double xarr[], 
                                 const double x){
      int i = 0;
      if ( x < xarr[0]) {
        return 0;
      } else if ( x > xarr[len-1]) {
        return len-2;
      } else {
        for (i=0; x>xarr[i]; i++) {}
        return i-1;
      }
    }

		static double LinearInterp(const double x0, const double x1,
                               const double y0, const double y1,
                               const double x){
      return y0 + ( (y1-y0)/(x1-x0) ) * (x-x0);
    }
};

#endif /*SHIELDING_H_*/
