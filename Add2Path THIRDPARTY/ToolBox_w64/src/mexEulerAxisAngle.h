/******************************************************************************
 *                  Mex Gateway routine for euler_axis_angle                  *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexEuler_axis_angle.c conversions.c math_utils.c -output euler_axis_angle -DMEXCOMPILE
*/

#ifndef MEXEULERAXISANGLE_H
#define MEXEULERAXISANGLE_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXEULERAXISANGLE_H */
