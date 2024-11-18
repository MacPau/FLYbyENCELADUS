/******************************************************************************
 *                    Mex Gateway routine for ephNEO                          *
 *                                                                            *
 *                                                      Matteo Ceriotti, 2009 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexEphNEO.c ephNEO.c mathUtils.c anomConv.c -output ephNEO -DMEXCOMPILE
*/

#ifndef MEXEPHNEO_H
#define MEXEPHNEO_H

/* Include */
#include "mex.h"
#include "matrix.h"
#include "ephNEO.h"
#include "string.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXEPHNEO_H */
