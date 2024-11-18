/******************************************************************************
 *                        Mex Gateway routine for kepEq_t                     *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexTl3.c orbital_transfers.c math_utils.c -output tl3 -DMEXCOMPILE
*/

#ifndef MEXKEPEQ_T_H
#define MEXKEPEQ_T_H

/* Include */
#include "mex.h"
#include "keplerianMotion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXKEPEQ_T_H */
