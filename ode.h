/* Date Mar 14, 2015  Author: Munan Gong
 * Base class Ode for ode solver
 */
#ifndef ODE_H_
#define ODE_H_

#include <sundials/sundials_types.h> /* realtype type*/
#include <nvector/nvector_serial.h> /* N_Vector type*/
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix  */
#include <stdio.h>
#include "sundial.h" /*Ith IJth macro and CheckFlag function*/
class Ode {
  friend class CvodeDense;
  public:
    /*constructor: set y_ to NULL pointer*/
    Ode();
    /*destructor: do nothing*/
    virtual ~Ode();
    /*set initial time and vector*/
    virtual void SetInit(const double t0, const double *y0);
    /*print status of t_ and y_*/
    virtual void PrintStat () const;
    /*Wrapper function for Jac and RHS to pass the member function pointer.
     *This is going to directly communiate with CVODE package.
     *Note that user_data is always set to pointer to Ode class.*/
    virtual void WriteStat (FILE *pf) const;
    static int wrapJac(const realtype t,
                       const N_Vector y, const N_Vector fy, 
                       SUNMatrix J, void *user_data,
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    static int wrapRHS(const realtype t, const N_Vector y,
                       N_Vector ydot, void *user_data);

    /*------------All functions below has to be overload------------*/
    /*return dimention*/
    virtual int Dimen() const = 0;
    /* Note that the RHS and Jac does NOT have user_data. All parameters should
     * be passed to the class as private variables.*/
    /* right hand side of ode */
    virtual int RHS(const realtype t, const N_Vector y,
                   N_Vector ydot) = 0;
    /*Jacobian*/
    virtual int Jac(const realtype t,
                    const N_Vector y, const N_Vector fy, 
                    SUNMatrix J, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) = 0;
  protected:
    /*status of ODE. Note the the ODE solver will change this. */
    realtype t_;
    N_Vector y_;
};

#endif /*ODE_H_*/

