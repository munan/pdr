#include "slab.h"

Slab::Slab(gow17 &ode, CvodeDense &solver,
		       const long int ngrid, const double NH_total,
           const double G0, const double Zd,
           const bool logNH, const double NH_min)
	:ode_(ode),
	 solver_(solver),
	 dimen_(ode.Dimen()),
	 nE_(ode.GetnE()),
	 ngrid_(ngrid),
	 NH_total_(NH_total),
   G0_(G0),
   logNH_(logNH),
   NH_min_(NH_min),
   field_geo_(0)
{
  prad_ = new RadField(ngrid, G0_, Zd);
	y_ = new double* [ngrid];
	yE_ = new double* [ngrid];
	for (int i=0; i<ngrid_; i++) {
		y_[i] = new double [dimen_];
		yE_[i] = new double [nE_];
	}
	/*initialize t and y to be the initial value in the ode*/
	t_ = ode.GetTime();
	fShieldH2mol_ = new double [ngrid_];
	fShieldCOmol_ = new double [ngrid_];
	hLast_ = new double [ngrid_];
	tSolve_ = new double [ngrid_];
	nstepLast_ = new double [ngrid_];
	NH_arr_ = new double [ngrid_];
	for (int i=0; i<ngrid_; i++) {
		ode.CopyAbd(y_[i]);
		fShieldH2mol_[i] = 0.;
		fShieldCOmol_[i] = 0.;
	}
	/*initialize yE_*/
	for (int i=0; i<ngrid_; i++) {
		for (int j=0; j<nE_; j++) {
			yE_[i][j] = 0.;
		}
	}
}

Slab::~Slab() {
	for (int i=0; i<ngrid_; i++) {
		delete [] y_[i];
		delete [] yE_[i];
	}
	delete [] y_;
	delete [] yE_;
	delete [] fShieldH2mol_;
	delete [] fShieldCOmol_;
	delete [] hLast_;
	delete [] tSolve_;
	delete [] nstepLast_;
	delete [] NH_arr_;
  delete prad_;
}

void Slab::WriteAbd(FILE *pf) {
	for (int i=0; i<ngrid_; i++) {
		for (int j=0; j<dimen_; j++) {
			fprintf(pf, "%12.4e  ", y_[i][j]);
		}
		fprintf(pf, "\n");
	}
}

void Slab::WriteThermoRates(FILE *pf) {
	for (int i=0; i<ngrid_; i++) {
		for (int j=0; j<nE_; j++) {
			fprintf(pf, "%12.4e  ", yE_[i][j]);
		}
		fprintf(pf, "\n");
	}
}


void Slab::SolveEq(const double tolfac, const double tmin, 
		               const double tmax, 
		               const bool verbose, FILE *pf_rates) {
  double dNH = 0;
  double NH = 0.;
	double NH2 = 0.;
	double NCO = 0.;
	double NC = 0.;
  double xCI = 0.;
  double fac = 0.;
  double xCtot = ode_.GetxCtot();
  /* create array of NH */
  if (logNH_) {
    NH_arr_[0] = NH_min_;
    fac = pow( 10, ( log10(NH_total_) - log10(NH_min_) )/ngrid_ );
    for (int i=1; i<ngrid_; i++) {
      NH_arr_[i] = NH_arr_[i-1]*fac;
    }
  } else {
    fac = NH_total_/ngrid_;
    for (int i=0; i<ngrid_; i++) {
      NH_arr_[i] = fac * (i + 1);
    }
  }

	for (int i=0; i<ngrid_; i++) {
		NH = NH_arr_[i];
    if (i == ngrid_ - 1) {
      dNH = 0.;
    } else {
      dNH = NH_arr_[i+1] - NH;
    }
		ode_.SetNCO(NCO);
		ode_.SetNH(NH);
    /*solve radiation field*/
    if (field_geo_ == 0) {
      prad_->Beamed(i, NH, NH2, NCO, NC);
    } else if (field_geo_ == 1) {
      prad_->IsotropicApprox(i, NH, NH2, NCO, NC);
    } else if (field_geo_ == 2) {
      prad_->Isotropic(i, NH, NH2, NCO, NC, 20);
    } else {
      printf("field_geo_ = %d, not recogonized, use beamed geometery.\n",
             field_geo_);
    }
    /*assign radiation field to chemistry*/
    ode_.SetRadField(prad_->GPE + i, *(prad_->Gph + i), prad_->GISRF + i);
    /*get sheilding factor*/
		fShieldH2mol_[i] = prad_->GetfShieldH2mol();
		fShieldCOmol_[i] = prad_->GetfShieldCOmol();
    /*solve to equalibrium*/
		solver_.SolveEq(tolfac, tmax, verbose, tmin);
		ode_.CopyAbd(y_[i]);
		ode_.CopyThermoRates(yE_[i]);
		NH2 += y_[i][ode_.id("H2")] * dNH;
		NCO += y_[i][ode_.id("CO")] * dNH;
    xCI =  xCtot - y_[i][ode_.id("HCO+")] - y_[i][ode_.id("CHx")] 
                 - y_[i][ode_.id("CO")]- y_[i][ode_.id("C+")];
		NC += xCI * dNH;
		hLast_[i] = solver_.GethLast();
		tSolve_[i] = solver_.GettSolve();
		nstepLast_[i] = solver_.GetnstepLast();
		if (pf_rates != NULL) {
			ode_.WriteRates(pf_rates);
		}
	}
}

void Slab::WritefShieldH2mol(FILE *pf) {
	for (int i=0; i<ngrid_; i++) {
		fprintf(pf, "%12.4e  ", fShieldH2mol_[i]);
	}
	fprintf(pf, "\n");
}

void Slab::WritefShieldCOmol(FILE *pf) {
	for (int i=0; i<ngrid_; i++) {
		fprintf(pf, "%12.4e  ", fShieldCOmol_[i]);
	}
	fprintf(pf, "\n");
}

void Slab::WriteNH(FILE *pf) {
	for (int i=0; i<ngrid_; i++) {
		fprintf(pf, "%12.4e  ", NH_arr_[i]);
	}
	fprintf(pf, "\n");
}

void Slab::WriteProf(FILE *pf) {
	fprintf(pf, "tSolve = ");
	for (int i=0; i<ngrid_; i++) {
		fprintf(pf, "%12.4e  ", tSolve_[i]);
	}
	fprintf(pf, "\n");

	fprintf(pf, "nstepLast = ");
	for (int i=0; i<ngrid_; i++) {
		fprintf(pf, "%12.4e  ", nstepLast_[i]);
	}
	fprintf(pf, "\n");

	fprintf(pf, "hLast = ");
	for (int i=0; i<ngrid_; i++) {
		fprintf(pf, "%12.4e  ", hLast_[i]);
	}
	fprintf(pf, "\n");
}

void Slab::IsDustSheilding(const bool isdust) {
	prad_->IsdustSheild(isdust);
	return;
}

void Slab::IsH2MolSheilding(const bool isH2mol) {
	prad_->IsMolSheildH2(isH2mol);
	return;
}

void Slab::IsCOMolSheilding(const bool isCOmol) {
	prad_->IsMolSheildCO(isCOmol);
	return;
}

void Slab::IsCselfSheilding(const bool isCfs) {
	prad_->IsSelfSheildC(isCfs);
	return;
}

void Slab::SetFieldGeo(int field_geo) {
  field_geo_ = field_geo;
  return;
}
