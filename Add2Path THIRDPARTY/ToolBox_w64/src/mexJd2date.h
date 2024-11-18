/******************************************************************************
 *                       Mex Gateway routine for jd2date                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexJd2date.c conversions.c math_utils.c -output jd2date -DMEXCOMPILE
*/

#ifndef MEXJD2DATE_H
#define MEXJD2DATE_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXJD2DATE_H */
