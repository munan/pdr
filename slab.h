/*Date: June 2, 2015, Author: Munan Gong
 * Grid of a plane-parallel slab.
 * TODO: implement both beamed and isotropic radiation.
 */

#ifndef SLAB_H_
#define SLAB_H_

#include "NL99p.h"
#include "cvodeDense.h"
#include <stdio.h>

class Slab {
  friend class CoolingFunction;
	public:
		Slab(NL99p &ode, CvodeDense &solver,
				 const long int ngrid, const double NH_total, const double G0,
         const double Zd, const bool logNH=false, const double NH_min=1.0e18);
		~Slab();
		/*Solve the chemistry to equalibrium*/
		void SolveEq(const double tolfac, const double tmin, const double tmax, 
								 const bool verbose, FILE *pf_rates=NULL);
		/*turn on/off dust sheilding*/
		void IsDustSheilding(const bool isdust);
		/*turn on/off H2 self sheilding*/
		void IsH2MolSheilding(const bool isH2mol);
		/*turn on/off CO sheilding by H2 and self*/
		void IsCOMolSheilding(const bool isCOmol);
		void IsCselfSheilding(const bool isCfs);
		/*Output abundances to file, ngrid line, dimen column*/
		void WriteAbd(FILE *pf);
		/*Output array of H2 self sheilding factor*/
		void WritefShieldH2mol(FILE *pf);
		/*Output array of CO sheilding factor by CO and H2*/
		void WritefShieldCOmol(FILE *pf);
    /* Output array of NH*/
    void WriteNH(FILE *pf);
		/*Output parameter array to file*/
		void WriteProf(FILE *pf);
		/*Output rates of thermo processes to file*/
		void WriteThermoRates(FILE *pf);
    /*Set radiation field geometery*/
    void SetFieldGeo(int field_geo);

	private:
		/*TODO: think about new data structure of only have one ode system,
		 * but store chemical species in an array that belongs to the grid*/
		NL99p &ode_;
		CvodeDense &solver_;
    RadField *prad_;
		const int dimen_; /*dimention of ode*/
		const int nE_; /*number of heating and cooling processes*/
		/*slab structure*/
		const long int ngrid_; /*number of grids in the slab*/
		const double NH_total_; /*Total column density in the slab*/
		/*chemical status*/
		double t_;
		double **y_;
		double **yE_; /*heating and cooling processes*/
		double *fShieldH2mol_;
		double *fShieldCOmol_;
		double *hLast_;
		double *tSolve_;
		double *nstepLast_;
		double *NH_arr_;
		/*calculate with/without dust/self sheilding*/
    const double G0_;
    /*whether we use logarithm spacing for NH and what is the minimum if so*/
    const bool logNH_;
    const double NH_min_;
    /*field geometry: 
     * 0: beamed, 1: IsotropicApprox, 2: Isotropic*/
    int field_geo_;
};

#endif /*SLAB_H_*/
