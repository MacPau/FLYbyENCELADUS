/******************************************************************************
 *                      Mex Gateway routine for kepEq_f                       *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexKep_eq1.c orbital_transfers.c math_utils.c -output kep_eq1 -DMEXCOMPILE
*/

#ifndef MEXKEPEQ_F_H
#define MEXKEPEQ_F_H

/* Include */
#include "mex.h"
#include "keplerianMotion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXKEPEQ_F_H */
