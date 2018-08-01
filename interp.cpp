#include "interp.h"

Interp::Interp() {}
Interp::~Interp() {}

double Interp::LP1D(const int len, const double *xarr, const double *data,
                    const double x) {
  const int ix = LinearInterpIndex_(len, xarr, x);
  return LP1Di_(xarr, data, ix, x);

}

double Interp::LP2D(const int lenx, const double *xarr, 
                    const int leny, const double *yarr,
                    const double *data, const double x, const double y) {
  const int ix = LinearInterpIndex_(lenx, xarr, x);
  const int iy = LinearInterpIndex_(leny, yarr, y);
  return LP2Di_(xarr, yarr, lenx, ix, iy, data, x, y);
}

double Interp::LP1Di_(const double *xarr, const double *data, const int ix,
                     const double x) {
  return LinearInterp_(xarr[ix], xarr[ix+1], data[ix], data[ix+1], x);
}

double Interp::LP2Di_(const double *xarr, const double *yarr,
                      const int lenx, const int ix, const int iy,
                      const double *data, const double x, const double y) {
  double fl1, fl2;
  const double x0 = xarr[ix];
  const double x1 = xarr[ix+1];
  fl1 = LinearInterp_(x0, x1, data[iy*lenx + ix], data[iy*lenx + ix+1], x);
  fl2 = LinearInterp_(x0, x1, data[(iy+1)*lenx + ix], data[(iy+1)*lenx + ix+1], x);
  return LinearInterp_(yarr[iy], yarr[iy+1], fl1, fl2, y);
}
