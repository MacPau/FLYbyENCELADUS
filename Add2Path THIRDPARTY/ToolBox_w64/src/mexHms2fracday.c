/******************************************************************************
 *                     Mex Gateway routine for hms2fracday                    *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexHms2fracday.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double hrs;         /* Number of hours (input) */
    double mn;          /* Number of minutes (input) */
    double sec;         /* Number of seconds (input) */
    double fracday;     /* Fraction of day (output) */
    
    
    /* Check the number of inputs and outputs */
    if(nrhs != 3)
        mexErrMsgTxt("HMS2FRACDAY requires 3 and only 3 inputs.");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. HMS2FRACDAY returns only 1 output.");
    
    
    /* Get values of inputs */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetM(prhs[0])!=1
                             || mxGetN(prhs[0])!=1
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The first input of HMS2FRACDAY must be a real scalar!");
    hrs = *mxGetPr(prhs[0]);
        
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
                             || mxGetM(prhs[1])!=1
                             || mxGetN(prhs[1])!=1
                             || mxGetNumberOfDimensions(prhs[1]) != 2)
        mexErrMsgTxt("The second input of HMS2FRACDAY must be a real scalar!");
    mn = *mxGetPr(prhs[1]);
        
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])
                             || mxGetM(prhs[2])!=1
                             || mxGetN(prhs[2])!=1
                             || mxGetNumberOfDimensions(prhs[2]) != 2)
        mexErrMsgTxt("The third input of HMS2FRACDAY must be a real scalar!");
    sec = *mxGetPr(prhs[2]);
    
    
    /* Compute the output */
    fracday = hms2fracday(hrs,mn,sec);
    
    
    /* Return the output */
    plhs[0] = mxCreateDoubleScalar(fracday);

    

    /* Exit the function */
    return;
}
