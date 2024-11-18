/******************************************************************************
 *                        Mex Gateway routine for jd2mjd                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexJd2mjd.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double jd;          /* Date in JD (input) */
    double mjd;         /* Date in MJD (output) */
    
    
    /* Check the number of inputs and outputs */
    if(nrhs != 1)
        mexErrMsgTxt("JD2MJD requires 1 and only 1 input.");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. JD2MJD returns only 1 output.");
    
    
    /* Get values of input */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetM(prhs[0])!=1
                             || mxGetN(prhs[0])!=1
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The only input of JD2MJD must be a real scalar!");
    jd = *mxGetPr(prhs[0]);
    
    
    /* Compute the output */
    mjd = jd2mjd(jd);
    
    
    /* Return the output */
    plhs[0] = mxCreateDoubleScalar(mjd);
    

    /* Exit the function */
    return;
}
