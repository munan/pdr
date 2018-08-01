/*Date: July 22, 2015. Author: Munan Gong
 * Linear and bi-linear interpolations*/

#ifndef INTERP_H_
#define INTERP_H_
#include <math.h>

class Interp {
  friend class Thermo;
  public:
    Interp();
    ~Interp();
    /* ID arrray linear interpolation*/
    static double LP1D(const int len, const double *xarr, const double *data,
                       const double x);
    /* 2D array bi-linear interpolation*/
    static double LP2D(const int lenx, const double *xarr, 
                       const int leny, const double *yarr,
                       const double *data, const double x, const double y);
  private:
		static int LinearInterpIndex_(const int len, const double xarr[], 
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

		static double LinearInterp_(const double x0, const double x1,
                               const double y0, const double y1,
                               const double x){
      return y0 + ( (y1-y0)/(x1-x0) ) * (x-x0);
    }

    /* Interpolation with index provided.
     * ix, iy: return from LinearInterpIndex */
    static double LP1Di_(const double *xarr, const double *data, const int ix,
                        const double x);
    /* 2D array bi-linear interpolation*/
    static double LP2Di_(const double *xarr, const double *yarr,
                         const int lenx, const int ix, const int iy,
                         const double *data, const double x, const double y);
};

#endif /*INTERP_H_*/
