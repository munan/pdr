#include <stdio.h>
#include "cool_tigress.h"

int main() {
  //output directory TODO: mkdir
	const char dir[] 
	 =
   "/Users/munangong/chemistry_Athena/PDR_cvode/tigress/out_sn/";
  //physical parameters TODO: check this, Setdvdr(true) SetNCOeffGlobal(false)
  const Real dvdr = 3 * 3e-14;
  const Real Z = 1;
  const Real xi_CR = 2e-16;
  const Real G0 = 1.;

  long int inH = 0;
  const Real G_PE = G0;
  const Real G_CI = G0;
  const Real G_CO = G0;
  const Real G_H2 = 0.;
	const Real nH_first = 1e-6;
	const Real nH_last = 1.0e4;
	const Real nH_fac = 10;
	const Real T_first = 10;
	const Real T_last = 1.0e9;
	const Real T_fac = 1.1;
	char fn_nH[100];  
  sprintf(fn_nH, "%snH_arr.dat", dir);
	char fn_T[100];  
  sprintf(fn_T, "%sT_arr.dat", dir);
  char file_out[100];
  FILE *pf_out = NULL;
  Real nH, T;
  Real x_e, x_HI, x_H2, x_Cplus, x_CI, x_CO, x_OI;
  Real h_CR, h_PE, h_H2pump, h_tot;
  Real c_Lya, c_OI, c_CII, c_CI, c_CO, c_Rec, c_tot;

	for (nH=nH_first; nH<nH_last; nH*=nH_fac) {
		sprintf(file_out, "%s%06ld.dat", dir, inH);
		pf_out = fopen(file_out, "w+");
    for (T=T_first; T<T_last; T*=T_fac) {
      get_abundances(nH, T, dvdr, Z, xi_CR, G_PE, G_CI, G_CO, G_H2,
                     &x_e, &x_HI, &x_H2, &x_Cplus, &x_CI, &x_CO, &x_OI);
      h_tot = heating(x_e, x_HI, x_H2, nH, T, Z, xi_CR, G_PE, G_H2);
      c_tot = cooling(x_e, x_HI, x_H2, x_Cplus, x_CI, x_CO, x_OI, 
                      nH, T, dvdr, Z, G_PE);
      h_PE = heatingPE(x_e, nH, T, Z, G_PE);
      h_CR = heatingCR(x_e, x_HI, x_H2, nH, xi_CR);
      h_H2pump =  heatingH2pump(x_HI, x_H2, nH, T, G_H2);
      c_Lya = coolingLya(x_e, x_HI, nH, T);
      c_OI = coolingOI(x_e, x_OI, x_HI, x_H2, nH, T);
      c_CII = coolingCII(x_e, x_Cplus, x_HI, x_H2, nH, T);
      c_CI = coolingCI(x_e, x_CI, x_HI, x_H2, nH, T);
      c_CO = coolingCO(x_e, x_CO, x_HI, x_H2, nH, T, dvdr);
      c_Rec = coolingRec(x_e, nH, T, Z, G_PE);
			fprintf(pf_out, "%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e %12.4e %12.4e\n",
             x_e, x_HI, x_H2, x_Cplus, x_CI, x_CO, x_OI, h_CR, h_PE, h_H2pump, h_tot,
             c_Lya, c_OI, c_CII, c_CI, c_CO, c_Rec, c_tot);
    }
    fclose(pf_out);
    inH++;
  }

	/*write range of nH in file*/
	FILE *pf_nH = fopen(fn_nH, "w+");
	for (nH=nH_first; nH<nH_last; nH*=nH_fac) {
		fprintf(pf_nH, "%12.4e     ", nH);
	}
	fclose(pf_nH);

	/*write range of T in file*/
	FILE *pf_T = fopen(fn_T, "w+");
	for (T=T_first; T<T_last; T*=T_fac) {
		fprintf(pf_T, "%12.4e     ", T);
	}
	fclose(pf_T);

  return 0;
}
