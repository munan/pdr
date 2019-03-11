#include "cvodeDense.h"

CvodeDense::CvodeDense(Ode &ode, 
                       const double reltol, const double *abstol,
											 const bool userJac)
    :ode_(ode),
     kDimen(ode.Dimen()),
     reltol_(reltol),
     abstol_(NULL),
     cvode_mem_(NULL),
     dense_matrix_(NULL),
     dense_ls_(NULL),
     t_solve_(0.),
     t_solve_eq_(0.),
     nst_last_(0),
     nst_eq_(0) {
  int flag;
  /* -----------Initialize absolute value vector----------*/
  abstol_ = N_VNew_Serial(kDimen);
  CheckFlag((void *)abstol_, "N_VNew_Serial", 0);
  for (int i=0; i<kDimen; i++) {
    NV_Ith_S(abstol_, i) = abstol[i];
  }

  /* ----------------Intialize solver---------------------*/
  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem_ = CVodeCreate(CV_BDF);
  CheckFlag((void *)cvode_mem_, "CVodeCreate", 0);

  /* Set the user data pointer to ode*/
  flag = CVodeSetUserData(cvode_mem_, &ode_);
  CheckFlag(&flag, "CVodeSetUserData", 1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem_, ode_.wrapRHS, 
                   ode_.t_, ode_.y_);
  CheckFlag(&flag, "CVodeInit", 1);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem_, reltol_, abstol_);
  CheckFlag(&flag, "CVodeSVtolerances", 1);

  /* Call CVDense to specify the CVDENSE dense linear solver */
  //flag = CVDense(cvode_mem_, kDimen);
  //CheckFlag(&flag, "CVDense", 1);

  /* Create dense SUNMatrix for use in linear solves */
  dense_matrix_ = SUNDenseMatrix(kDimen, kDimen);
  CheckFlag((void *)dense_matrix_, "SUNDenseMatrix", 0);

  /* Create dense SUNLinearSolver object for use by CVode */
  dense_ls_ = SUNDenseLinearSolver(ode_.y_, dense_matrix_);
  CheckFlag((void *)dense_ls_, "SUNDenseLinearSolver", 0);

  /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
  flag = CVDlsSetLinearSolver(cvode_mem_, dense_ls_, dense_matrix_);
  CheckFlag(&flag, "CVDlsSetLinearSolver", 1);

  /* Set the Jacobian routine to Jac (user-supplied) */
	if (userJac) {
		flag = CVDlsSetJacFn(cvode_mem_, ode_.wrapJac);
		CheckFlag(&flag, "CVDlsSetJacFn", 1);
	}
}

void CvodeDense::SetStabLimDet(const bool stldet) {
  int flag;
  flag = CVodeSetStabLimDet(cvode_mem_, stldet);
  CheckFlag(&flag, "CVodeSetStabLimDet", 1);
  return;
}

void CvodeDense::SetMaxOrd(const int maxord) {
  int flag;
  flag = CVodeSetMaxOrd(cvode_mem_, maxord);
  CheckFlag(&flag, "CVodeSetMaxOrd", 1);
  return;
}

void CvodeDense::SetMxsteps(const int mxsteps) {
  int flag;
  flag = CVodeSetMaxNumSteps(cvode_mem_, mxsteps);
  CheckFlag(&flag, "CVodeSetMaxNumSteps", 1);
  return;
}

void CvodeDense::SetInitStep(const double hinit) {
  int flag;
  flag = CVodeSetInitStep(cvode_mem_, hinit);
  CheckFlag(&flag, "CVodeSetInitStep", 1);
  return;
}


CvodeDense::~CvodeDense() {
  /*Free abstol_ vector*/
  N_VDestroy_Serial(abstol_);
  /* Free integrator memory */
  CVodeFree(&cvode_mem_);
}

void CvodeDense::Solve(const double tfinal) {
  int flag;
  clock_t t = clock();
  long int nst = GetNstep();
  flag = CVode(cvode_mem_, tfinal, ode_.y_, &ode_.t_,
               CV_NORMAL);
  CheckFlag(&flag, "CVode", 3);
  t = clock() - t;
  t_solve_ = ((float)t)/CLOCKS_PER_SEC; /*record the run time*/ 
  nst_last_ = GetNstep() - nst;
  /*floor small negative abundances to avoid numerical problems*/
  /*TODO: adjust species abundances to conserve the elements*/
  for (int i=0; i<kDimen; i++) {
    if ( NV_Ith_S(ode_.y_, i) < 0. ) {
      printf("Warning: correcting negative aboundances of species ");
      //printf("%d from %0.4e to abstol_[i] = %0.4e.\n", i,  
      //       NV_Ith_S(ode_.y_, i), NV_Ith_S(abstol_, i));
      printf("%d from %0.4e to zero.\n", i,  
             NV_Ith_S(ode_.y_, i));
      //NV_Ith_S(ode_.y_, i) = NV_Ith_S(abstol_, i);
      NV_Ith_S(ode_.y_, i) = 0.;
    }
  }
  return;
}
void CvodeDense::SolveEq(const double tolfac, const double tmax,
												 const bool verbose, const double tmin) {
  /* Sanity check*/
  if (tolfac <= 1.) {
    throw std::logic_error("ERROR: CvodeDense::SoveEq - tolfac <= 1. ");
  }
  if (tmax <= 0.) {
    throw std::logic_error("ERROR: CvodeDense::SoveEq - tmax <= 0. ");
  }

  /* Guess the time toward equalibrium.
   *  To make a guess, we compute dx/dt at the initial conditions and at a point 
   *  slightly perturbed from them, and then use the most restrictive timestep we find.
   *  We also require the minmum time of 10yr=3.16e8s and maximum time of
   *  tmax/10 */
  /* Compute current time derivatives */
	int flag;
  const double tevol_min = tmin;
  const double tevol_max = tmax / 16.;
  const double small__ = 1.e-50;
  double dt1 = tevol_max;
  double dt2 = tevol_max;
  double tevol = 0.;
  N_Vector ydot = N_VNew_Serial(kDimen);
  N_Vector ypert = N_VNew_Serial(kDimen);
  /* get the minimum component of dt1, restrict to be smaller than tevol_max*/
  ode_.RHS(ode_.t_, ode_.y_, ydot);
  for (int i=0; i<kDimen; i++) {
    dt1 = std::min( dt1, fabs( 
               ( NV_Ith_S(ode_.y_, i) + NV_Ith_S(abstol_, i) )  
                / ( NV_Ith_S(ydot, i) + small__ ) )
                  );
  }
  /*set the perturbed y*/
  for (int i=0; i<kDimen; i++) {
    NV_Ith_S(ypert, i) = NV_Ith_S(ode_.y_, i) + dt1*NV_Ith_S(ydot, i);
  }
  /*calculate the minimum component of the perturbed y, restrict to be smaller
   * than tevol_max*/
  ode_.RHS(ode_.t_+dt1, ypert, ydot);
  for (int i=0; i<kDimen; i++) {
    dt2 = std::min( dt2, fabs( 
               ( NV_Ith_S(ypert, i) + NV_Ith_S(abstol_, i) )  
                / ( NV_Ith_S(ydot, i) + small__ ) )
                  );
  }
  /*Set the guess for equalibrium timescale*/
  tevol = std::max(dt1, dt2)/10.; /*conservative evaluation*/
  /*make sure tevol is above tevol_min*/
  tevol = std::max(tevol, tevol_min);
  if (verbose) {
    printf("SolveEq: estimated equilibration timescale = %0.4e seconds.\n", tevol);
  }

  /*destroy created vectors*/
  N_VDestroy_Serial(ydot);
  N_VDestroy_Serial(ypert);

  /*find equalibrium.*/
  N_Vector yprev = N_VNew_Serial(kDimen);
  double err[kDimen];
  double maxerr;
  clock_t t = clock();
  long int nst = GetNstep();
  while (true) {
		/*reinitialize ode*/
		flag = CVodeReInit(cvode_mem_, ode_.t_, ode_.y_);
		CheckFlag(&flag, "CVodeReInit", 1);
		
    /*copy y to yprev*/
    for (int i=0; i<kDimen; i++) {
      NV_Ith_S(yprev, i) = NV_Ith_S(ode_.y_, i);
    }
    /*solve y to t+tevol*/
    Solve(ode_.t_ + tevol);
    /*Calculate the relative differences of y from last timestep*/
    for (int i=0; i<kDimen; i++) {
      err[i] = fabs(  NV_Ith_S(ode_.y_, i) / ( NV_Ith_S(yprev, i) + small__ ) 
          - 1.  );
      /*not consider aboundances below 10 times absolute tolerance*/
      if ( NV_Ith_S(ode_.y_, i) < 10*NV_Ith_S(abstol_, i) ) {
        err[i] = 0.;
      }
    }
    /*test equalibrium statndard*/
    maxerr = *std::max_element(err, err+kDimen);
    if (verbose) {
      printf("SolveEq: evolved from t = %0.4e to %0.4e seconds, residual = %0.4e.\n",
              ode_.t_ - tevol, ode_.t_, maxerr);
    }
    if ( maxerr < reltol_*tolfac ) {
      N_VDestroy_Serial(yprev);
      t = clock() - t;
      t_solve_eq_ = ((float)t)/CLOCKS_PER_SEC; /*record the run time*/ 
      nst_eq_ = GetNstep() - nst;
      return;
    }
    if (tevol > tevol_max) {
      N_VDestroy_Serial(yprev);
      t = clock() - t;
      t_solve_eq_ = ((float)t)/CLOCKS_PER_SEC; /*record the run time*/ 
      nst_eq_ = GetNstep() - nst;
      printf(
			 "Warning: not reaching equalibrium in tevol_max, residual=%0.4e.\n", 
			 maxerr);
      return;
    }
    tevol *= 2.;
  }

}
  
void CvodeDense::PrintStats() const {
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem_, &nst);
  CheckFlag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem_, &nfe);
  CheckFlag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem_, &nsetups);
  CheckFlag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem_, &netf);
  CheckFlag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem_, &nni);
  CheckFlag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem_, &ncfn);
  CheckFlag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem_, &nje);
  CheckFlag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem_, &nfeLS);
  CheckFlag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem_, &nge);
  CheckFlag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nCvodeDens: Final Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
  return;
}

void CvodeDense::WriteRunTime(FILE *pf) const {
  fprintf(pf, "tSolve = %0.4e\n", t_solve_);
  fprintf(pf, "tSolveEq = %0.4e\n", t_solve_eq_);
  return;
}
void CvodeDense::WriteNstep(FILE *pf) const {
  fprintf(pf, "nstepLast = %ld\n", nst_last_);
  fprintf(pf, "nstepEq = %ld\n", nst_eq_);
  return;
}
void CvodeDense::WriteFirstLastStep(FILE *pf) const {
  double hfirst, hlast;
  int flag;
  flag = CVodeGetActualInitStep(cvode_mem_, &hfirst);
  CheckFlag(&flag, "CVodeGetActualInitStep", 1);
  flag = CVodeGetLastStep(cvode_mem_, &hlast);
  CheckFlag(&flag, "CVodeGetLastStep", 1);
  fprintf(pf, "hFirst = %0.4e\n", hfirst);
  fprintf(pf, "hLast = %0.4e\n", hlast);
  return;
}

int CvodeDense::GetNstep() const {
  int flag;
  long int nst = 0;
  flag = CVodeGetNumSteps(cvode_mem_, &nst);
  CheckFlag(&flag, "CVodeGetNumSteps", 1);
  return nst;
}

void CvodeDense::ReInit() {
	int flag;
	flag = CVodeReInit(cvode_mem_, ode_.t_, ode_.y_);
	CheckFlag(&flag, "CVodeReInit", 1);
}

double CvodeDense::GethLast() const {
	double hlast;
	int flag;
  flag = CVodeGetLastStep(cvode_mem_, &hlast);
  CheckFlag(&flag, "CVodeGetLastStep", 1);
	return hlast;
}

double CvodeDense::GettSolve() const {
	return t_solve_;
}

double CvodeDense::GetnstepLast() const {
	return nst_last_;
}
