/* Date: Mar 4, 2015, author: Munan Gong
 * sundial.h
 * -------------------------------------------------------------------
 * Macros and functions that are useful in Sundial package.
 * Adopted from examples in cvode/cvRoberts_dns.*/

#ifndef SUNDIAL_H_
#define SUNDIAL_H_

#include <stdio.h>
#include <nvector/nvector_serial.h> /* N_Vector type*/
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <cvode/cvode.h> /* CV_SUCCESS */
#include <stdexcept> /*throw exceptions*/

/*Access or assign N_Vector or DlsMat elements*/
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

/*check flags function for Sundial*/
void CheckFlag(const void *flagvalue, const char *funcname, 
               const int opt);

#endif /*SUNDIAL_H_*/
