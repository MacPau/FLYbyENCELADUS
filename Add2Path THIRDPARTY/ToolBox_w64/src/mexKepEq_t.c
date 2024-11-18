/******************************************************************************
 *                        Mex Gateway routine for tl3                         *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/
 
#include "mexKepEq_t.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declarations */
    double f;       /* Double containing the the true anomaly [rad] (input) */
    double a;       /* Double containing the semi-major axis [L] (input) */
    double ecc;     /* Double containing the eccentricity (input) */
    double mu;      /* Double containing the planetary constant (input)
                       mu = mass * G) [L^3/T^2]. */
    double f0;      /* Double containing a given true anomaly [rad] (input) */
    double t0;      /* Double containing the time corresponding to f0 [T]
                       (optional input) */
    double out;     /* Double containing the time corresponding to f [T] between
                       [-Inf, Inf] (output) */
    
    /* Check the number of inputs and outputs */
    if((nrhs < 5) || nrhs > 6)
        mexErrMsgTxt("KepEq_t requires 5 or 6 inputs.");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. KepEq_t returns only 1 output.");
    
    
    /* Get values of inputs */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetM(prhs[0])!=1
                             || mxGetN(prhs[0])!=1
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The first input of KepEq_t must be a real scalar!");
    f = mxGetScalar(prhs[0]);

    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
                             || mxGetM(prhs[1])!=1
                             || mxGetN(prhs[1])!=1
                             || mxGetNumberOfDimensions(prhs[1]) != 2)
        mexErrMsgTxt("The second input of KepEq_t must be a real scalar!");
    a = mxGetScalar(prhs[1]);
    
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])
                             || mxGetM(prhs[2])!=1
                             || mxGetN(prhs[2])!=1
                             || mxGetNumberOfDimensions(prhs[2]) != 2)
        mexErrMsgTxt("The third input of KepEq_t must be a real scalar!");
    ecc = mxGetScalar(prhs[2]);
 
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])
                             || mxGetM(prhs[3])!=1
                             || mxGetN(prhs[3])!=1
                             || mxGetNumberOfDimensions(prhs[3]) != 2)
        mexErrMsgTxt("The fourth input of KepEq_t must be a real scalar!");
    mu = mxGetScalar(prhs[3]);
    
    if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])
                             || mxGetM(prhs[4])!=1
                             || mxGetN(prhs[4])!=1
                             || mxGetNumberOfDimensions(prhs[4]) != 2)
        mexErrMsgTxt("The fifth input of KepEq_t must be a real scalar!");
    f0 = mxGetScalar(prhs[4]);
    
    /* Optional input */
    if(nrhs < 6)
        t0 = 0;     /* Default value for t0 */
    else
    {
        if (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5])
                                 || mxGetM(prhs[5])!=1
                                 || mxGetN(prhs[5])!=1
                                 || mxGetNumberOfDimensions(prhs[5]) != 2)
            mexErrMsgTxt("The sixth input of KepEq_t must be a real scalar!");
        t0 = mxGetScalar(prhs[5]);
    }
        
    /* Compute the output and check for error */
    if(kepEq_t(f, a, ecc, mu, f0, t0, &out)==1)
        mexErrMsgTxt("KepEq_t: The eccentricity is greater than 1!");
    
    /*Allocate memory and assign out value to the output pointer */
    plhs[0] = mxCreateDoubleScalar(out);
    
    
    /* Exit the function */
    return;
}
