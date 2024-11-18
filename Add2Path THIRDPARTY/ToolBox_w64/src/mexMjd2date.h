/******************************************************************************
 *                      Mex Gateway routine for mjd2date                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexMjd2date.c conversions.c math_utils.c -output mjd2date -DMEXCOMPILE
*/

#ifndef MEXMJD2DATE_H
#define MEXMJD2DATE_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXMJD2DATE_H */
