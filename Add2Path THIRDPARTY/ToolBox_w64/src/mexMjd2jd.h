/******************************************************************************
 *                        Mex Gateway routine for mjd2jd                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexMjd2jd.c conversions.c math_utils.c -output mjd2jd -DMEXCOMPILE
*/

#ifndef MEXMJD2JD_H
#define MEXMJD2JD_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXMJD2JD_H */
