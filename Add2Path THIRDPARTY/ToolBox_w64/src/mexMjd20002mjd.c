/******************************************************************************
 *                     Mex Gateway routine for mjd20002mjd                    *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexMjd20002mjd.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double mjd2000;     /* Date in MJD2000 (input) */
    double mjd;         /* Date in MJD (output) */
    
    
    /* Check the number of inputs and outputs */
    if(nrhs != 1)
        mexErrMsgTxt("MJD20002MJD requires 1 and only 1 input.");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. MJD20002MJD returns only 1 output.");
    
    
    /* Get values of input */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetM(prhs[0])!=1
                             || mxGetN(prhs[0])!=1
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The only input of MJD20002MJD must be a real scalar!");
    mjd2000 = *mxGetPr(prhs[0]);
    
    
    /* Compute the output */
    mjd = mjd20002mjd(mjd2000);
    
    
    /* Return the output */
    plhs[0] = mxCreateDoubleScalar(mjd);
    

    /* Exit the function */
    return;
}
