/* Date: Oct 2, 2015, Author: Munan Gong
 * Radiation field for gow17 network
 */
#ifndef RADFIELD_H_
#define RADFIELD_H_

#include <stdio.h>
#include "thermo.h"
#include "gow17.h"
#include <math.h>

const double pi = 3.1415926535897;

/*ISRF strengh, normalized to Solar neighborhood value.
 *In Draine 1978 field, G0 = 1.*/
class RadField {
  friend class gow17;
  public:
    RadField(const long int ngrid, const double G0, const double Zd);
    ~RadField();
    /*photo-electric heating on dust and grain assisted recombination of H+*/
    double *GPE; 
    /*ISRF*/
    double *GISRF; 
    /*photo reactions of gow17 network, including dust and self shielding for
     * CO and H2, Gph[igrid][ispec]*/
    double **Gph;
    /*calculate radiation field, and store in igrid*/
    void Beamed(const long int igrid, const double NH, 
                const double NH2, const double NCO, const double NC);
    /* approximation of isotropic field in WHM2010*/
    void IsotropicApprox(const long int igrid, const double NH, 
                         const double NH2, const double NCO, const double NC);
    /* isotropic aprroximation in angles*/
    void Isotropic(const long int igrid, const double NH, 
                   const double NH2, const double NCO, const double NC,
                   const long int nangle);
		/*turn on/off self shielding (and sheilding by H2) of H2 and CO*/
		void IsMolSheildH2(const bool isfsH2);
		void IsMolSheildCO(const bool isfsCO);
		void IsSelfSheildC(const bool isfsC);
		void IsdustSheild(const bool isdust);
		void SetbH2(const double bH2);
    /*return self-sheilding factors*/
    double GetfShieldH2mol();
    double GetfShieldCOmol();
  private:
    const double ngrid_;
    const double G0_; /* radiation field at boundary*/
		/* turn on/off self shielding (and sheilding by H2) of H2 and CO*/
		bool isfsH2_;
		bool isfsCO_;
		bool isfsC_;
    bool isdust_;
		/*Velocity dispersion of H2 used to calculate self-shielding of H2*/
		double bH2_;
    /* self sheilding (and sheilding by H2) factor of H2, CO, and C*/
    double fs_CO_;
    double fs_H2_;
    double fs_C_;
    double Zd_;
};

#endif /*RADFIELD_H_*/
