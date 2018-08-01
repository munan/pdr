/* Date: Mar 17, 2015 Author: Munan Gong
 * Implimentation of base class Ode for ode solver
 */
#include "ode.h"

Ode::Ode() {
  y_ = NULL;
}

Ode::~Ode() {}

void Ode::SetInit(const double t0, const double *y0) {
  int dimen = Dimen();
  t_ = t0;
  for (int i=0; i<dimen; i++) {
    NV_Ith_S(y_, i) = y0[i];
  }
  return;
}

void Ode::PrintStat () const {
  int dimen = Dimen();
  realtype y1, y2, y3;
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =", t_);
  for (int i=0; i<dimen; i++) {
    printf("%14.6Le  ", NV_Ith_S(y_,i));
  }
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4le      y =", t_);
  for (int i=0; i<dimen; i++) {
    printf("%14.6le  ", NV_Ith_S(y_,i));
  }
#else
  printf("At t = %0.4e      y =", t_);
  for (int i=0; i<dimen; i++) {
    printf("%14.6e  ", NV_Ith_S(y_,i));
  }
#endif
  printf("\n");
  return;
}

void Ode::WriteStat (FILE *pf) const {
  int dimen = Dimen();
  realtype y1, y2, y3;
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(pf, "%0.4Le      ", t_);
  for (int i=0; i<dimen; i++) {
    fprintf(pf, "%14.6Le  ", NV_Ith_S(y_,i));
  }
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  fprintf(pf, "%0.4le      ", t_);
  for (int i=0; i<dimen; i++) {
    fprintf(pf, "%14.6le  ", NV_Ith_S(y_,i));
  }
#else
  fprintf(pf, "%0.4e      ", t_);
  for (int i=0; i<dimen; i++) {
    fprintf(pf, "%14.6e  ", NV_Ith_S(y_,i));
  }
#endif
  fprintf(pf, "\n");
  return;
}

int Ode::wrapJac(const realtype t,
                 const N_Vector y, const N_Vector fy, 
                 SUNMatrix J, void *user_data,
                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  Ode *meptr = (Ode*) user_data;
  meptr->Jac(t, y, fy, J, tmp1, tmp2, tmp3);
  return 0;
}
int Ode::wrapRHS(const realtype t, const N_Vector y,
                        N_Vector ydot, void *user_data) {
  Ode *meptr = (Ode*) user_data;
  meptr->RHS(t, y, ydot);
  return 0;
}
