/******************************************************************************
 *                      Mex Gateway routine for cartProd                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexCartprod.c math_utils.c -output cartprod -DMEXCOMPILE
*/

#ifndef MEXCARTPROD_H
#define MEXCARTPROD_H

/* Include */
#include "mex.h"
#include "mathUtils.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXCARTPROD_H */
