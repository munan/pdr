/*Date: June 9, 2015, Author: Munan Gong
 * CNM slab models for fixed nH/G0 = 23 and vary nH.
 * Use to make plot for chemical species abundances at grids of nH and NH.
 */

#include "gow17.h"
#include "cvodeDense.h"
#include "slab.h"
#include "math.h"
#include <stdio.h>

int main() {
  gow17 odegow17;
	const double t0 = 0;
	const int dim = odegow17.Dimen();
	double y0[dim]; 
  double abstol[dim];
  const double reltol = 1.0e-2;//relative tolerance of the ode solver 
	const double abstol0 = 1.0e-9;//absolute tolerance of the ode solver
	const int mxsteps = 5000000; //maximum number of steps in ode solver
	const double temp = 100; //initial temperature
	double nH = 0; //density variable
	double G0 = 0; //radiation field variable

	const long int ngrid = 2e3; //number of grids
	const double tolfac = 1e-1 / reltol; // tolerance for equilibrium
  //maximum time in seconds for the chemical evolution
	const double tmax = 2000 * 1e6 * 3.16e7; 
  //minimum time in seconds for the chemical evolution
	const double tmin = 0.1 * 1e6 * 3.16e7; 
	const double verbose = false; //print detailed information
  //flag for dust shielding
	const bool isdust = true;
  //flag for H2 self-shielding
	const bool isfsH2 = true;
  //flag for CO self-shielding
	const bool isfsCO = true;
  //flag for CI self-shielding
	const bool isfsC = true;
  //dust and gas metallicity
  const double Zdg = 1.;
  //minimum column of the PDR
	const double NH_min = 1.0e17/Zdg;
  //maximum column of the PDR
	const double NH_total = 1.0e22/Zdg;
  //logarithmically or linearly spaced grid
  const bool islogNH = true;
  //radiation field: 0: beamed/plane-parallel,
  //1: Isotropic approximation/plane-parallel with 60 degrees incident angle
  //2: Isotropic
  const int field_geo = 0; 

  //set up densities of a series of PDR
	const double nH_first = 10;
	const double nH_last = 1.0e3;
	const double nH_fac = 1.05;
	//const double nH2G0 = 94. / ( 1. + 3.1*pow(Zdg, 0.365) ); /*nH/G0 constant in CNM*/
	long int slab_id = 0;
  //output directory. Note this should be the same as that in examples.in
	const char dir[] = "../out_example_large/";
	char fn_nH[100];  
  sprintf(fn_nH, "%snH_arr.dat", dir);
	char fn_colH[100];  
  sprintf(fn_colH, "%scolH_arr.dat", dir);
	char fn_para[100];  
  sprintf(fn_para, "%spara.dat", dir);
	char fn_chemnet[100];  
  sprintf(fn_chemnet, "%schemnet.dat", dir);
	FILE *pf_chemnet = fopen(fn_chemnet, "w+");
  odegow17.PrintChemNet(pf_chemnet);
  fclose(pf_chemnet);

	/*initialize y0 and absolute tolerance*/
	for (int i=0; i<dim; i++) {
		abstol[i] = abstol0;
		y0[i] = 0.;
	}
	y0[odegow17.id("He+")] = 1.450654e-08;
	y0[odegow17.id("H3+")] = 2.681411e-07;
	y0[odegow17.id("C+")] = 0.0001;
	y0[odegow17.id("S+")] = 2.7e-6 * 0;
	y0[odegow17.id("Si+")] = 1.e-6 * 0;
	y0[odegow17.id("CO")] = 1.0e-7;
	y0[odegow17.id("H2")] = 0.1;
	y0[odegow17.id("E")] = odegow17.GetE(temp, 0.1, 0.);

	abstol[odegow17.id("He+")] = 1e-15;
	abstol[odegow17.id("OHx")] = 1.0e-15;
	abstol[odegow17.id("CHx")] = 1.0e-15;
	abstol[odegow17.id("CO")] = 1.0e-15;
	abstol[odegow17.id("C+")] = 1.0e-15;
	abstol[odegow17.id("HCO+")] = 1.0e-30;
	abstol[odegow17.id("H2")] = 1.0e-8;
	abstol[odegow17.id("H+")] = 1.0e-15;
	abstol[odegow17.id("H3+")] = 1.0e-15;
	abstol[odegow17.id("H2+")] = 1.0e-15;
	abstol[odegow17.id("E")] = odegow17.GetE(1., 0.1, 0);

  CvodeDense cvodeDense(odegow17, reltol, abstol);
	cvodeDense.SetMxsteps(mxsteps);
	cvodeDense.SetMaxOrd(3);

	//odegow17.SetConstTemp(temp); //set temperature to be constant
  //parameters for CO cooling
  odegow17.SetGradv(3 * 3e-14);
  odegow17.SetNCOeffGlobal(true);
  //odegow17.SetbCO(1 * 1e5);
  odegow17.SetbCOL(true);
  //CR ionisation rate
  odegow17.SetIonRate(2e-16);
  //dust and gas metallicity
  odegow17.SetZg(Zdg);

  //Set grain reaction rates
  odegow17.SetfH2gr(10./sqrt(temp));
  odegow17.SetfHplusgr(0.6);
  odegow17.SetfHeplusgr(0.6);
  odegow17.SetfCplusgr(0.6);
  odegow17.SetfSplusgr(0.6);
  odegow17.SetfSiplusgr(0.6);
  
	bool isWriteNH = true;
	for (nH=nH_first; nH<nH_last; nH*=nH_fac) {
		char *file_slab = new char[100];
		char *file_rates = new char[100];
		sprintf(file_slab, "%sslab%06ld.dat", dir, slab_id);
		sprintf(file_rates, "%srates%06ld.dat", dir, slab_id);
		FILE *pf = fopen(file_slab, "w+");
		FILE *pf_rates = fopen(file_rates, "w+");
		//G0 = nH / nH2G0;
    //odegow17.SetIonRate(2e-16 * sqrt(G0));
    G0 = 2.; //incident radiation field strength
		odegow17.SetInit(t0, y0);
		odegow17.SetnH(nH);
		Slab *myslab = new Slab(odegow17, cvodeDense, ngrid, NH_total, G0, Zdg,
                            islogNH, NH_min);
    myslab->SetFieldGeo(field_geo);
		myslab->IsH2MolSheilding(isfsH2);
		myslab->IsCOMolSheilding(isfsCO);
		myslab->IsCselfSheilding(isfsC);
		myslab->IsDustSheilding(isdust);
		myslab->SolveEq(tolfac, tmin, tmax, verbose, pf_rates);
		myslab->WriteAbd(pf);
		fclose(pf);
		fclose(pf_rates);
		delete [] file_slab;
		delete [] file_rates;
		/*write profiling parameters*/
		char *file_prof = new char[100];
		sprintf(file_prof, "%sprof%06ld.dat", dir, slab_id);
		FILE *pf_prof = fopen(file_prof, "w+");
		myslab->WriteProf(pf_prof);
		fclose(pf_prof);
		delete [] file_prof;
		/*write rates of thermo processes*/
		char *file_thermo = new char[100];
		sprintf(file_thermo, "%sthermo%06ld.dat", dir, slab_id);
		FILE *pf_thermo = fopen(file_thermo, "w+");
		myslab->WriteThermoRates(pf_thermo);
		fclose(pf_thermo);
		delete [] file_thermo;

    if (isWriteNH) {
      FILE *pf_colH = fopen(fn_colH, "w+");
      myslab->WriteNH(pf_colH);
      fclose(pf_colH);
      isWriteNH = false;
    }

		delete myslab;
		slab_id++;
	}

	/*write range of nH in file*/
	FILE *pf_nH = fopen(fn_nH, "w+");
	for (nH=nH_first; nH<nH_last; nH*=nH_fac) {
		fprintf(pf_nH, "%12.4e     ", nH);
	}
	fclose(pf_nH);


	/*write parameters of the last slab to file*/
	FILE *pf_para = fopen(fn_para, "w+");
	odegow17.WriteParam(pf_para);
	fprintf(pf_para, "reltol = %0.4e\n", reltol);
	fprintf(pf_para, "tolfac = %0.4e\n", tolfac);
	fprintf(pf_para, "tmin = %0.4e\n", tmin);
	fclose(pf_para);

	return 0;
}
