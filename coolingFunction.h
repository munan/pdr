/* Date Sept 5, 2018  Author: Munan Gong
 * Class cooling function. Read-in the abundances from Slab class and calculate
 * cooling. Prototype for cooling function in TIGRESS.
 */
#ifndef COOLINGFUNCTION_H_
#define COOLINGFUNCTION_H_

#include "slab.h"
#include "NL99p.h"

class CoolingFunction {
  public:
    CoolingFunction(Slab &myslab);
    ~CoolingFunction();
    void ComputeAbundances();
		void WriteAbundances(FILE *pf);
    void ComputeThermoRates();
		void WriteThermoRates(FILE *pf);
  private:
    Slab &myslab_;
    //input parameters
		const int nE_; /*number of heating and cooling processes*/
    double kcr_;
    double nH_;
    double Zd_;
    double xCtot_;
    double xOtot_;
    double xHetot_;
    long int ngrid_;
    double *T_;
    double *GPE_;
    double *GCI_;
    double *GCO_;

    //output paramters
    double *fe_;
    double *fHplus_;
    double *fH2_;
    double *fHI_;
    double *fCplus_;
    double *fCO_;
    double *fCI_;
    double *fOI_;
		double **yE_; /*heating and cooling processes*/

    //functions to compute abundances
    void get_fH2_();
    void get_fCO_(); //Its necessary to calculate fCplus_ and fH2_ first
    double fHplus_e_(double x_e, double x_H2, double temp, double G_PE);
    double CII_rec_rate_(const double temp);
    double fCplus_e_(double x_e, double x_H2, double temp, double G_PE, double G_CI);
    double fe_e_(double x_e, double x_H2, double temp, double G_PE, double G_CI);
    void get_fe_(); //need fH2_ first, solve the electron abundance with iteration
    void get_fHplus_();//need to calculate fe_ and fH2_ first
    void get_fCplus_();//need to calculate fe_ and fH2_ first
};

#endif /*COOLINGFUNCTION_H_*/
