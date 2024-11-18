/******************************************************************************
 *                      Mex Gateway routine for kep_eq1                       *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/
 
#include "mexKepEq_f.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declarations */
    double t;       /* Double containing the the time when the true anomaly is
                       required [T] (input) */
    double a;       /* Double containing the semi-major axis [L] (input) */
    double ecc;     /* Double containing the eccentricity (input) */
    double mu;      /* Double containing the planetary constant (input)
                       mu = mass * G) [L^3/T^2]. */
    double f0;      /* Double containing a given true anomaly [rad] (input) */
    double t0;      /* Double containing the time corresponding to f0 [T]
                       (optional input) */
    int option;     /* Integer containing the printing options
                       (optional input, default = 0) */
    int imax;       /* Integer containing the maximum number of iteration.
                       (optional input, default = 5) */
    double tol;     /* Double containing the convergence tolerance for the
                       Newton loop on abs(E-E0) [rad]
                       (optional input, default = 1e-15) */
    double out;     /* Double containing the the true anomaly at time t [rad]
                       in [-Inf, Inf] (output) */
    int error;      /* Integer containing the error code returned by kep_eq1
                       (optional output) */
    
    
    /* Check the number of inputs and outputs */
    if((nrhs < 6) || nrhs > 9)
        mexErrMsgTxt("KepEq_f requires 6 to 9 inputs.");
    if(nlhs > 2)
        mexErrMsgTxt("Too many outputs. KepEq_f returns 1 or 2 output(s).");
    
    
    /* Get values of inputs */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetM(prhs[0])!=1
                             || mxGetN(prhs[0])!=1
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The first input of KepEq_f must be a real scalar!");
    t = mxGetScalar(prhs[0]);

    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
                             || mxGetM(prhs[1])!=1
                             || mxGetN(prhs[1])!=1
                             || mxGetNumberOfDimensions(prhs[1]) != 2)
        mexErrMsgTxt("The second input of KepEq_f must be a real scalar!");
    a = mxGetScalar(prhs[1]);
    
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])
                             || mxGetM(prhs[2])!=1
                             || mxGetN(prhs[2])!=1
                             || mxGetNumberOfDimensions(prhs[2]) != 2)
        mexErrMsgTxt("The third input of KepEq_f must be a real scalar!");
    ecc = mxGetScalar(prhs[2]);
 
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])
                             || mxGetM(prhs[3])!=1
                             || mxGetN(prhs[3])!=1
                             || mxGetNumberOfDimensions(prhs[3]) != 2)
        mexErrMsgTxt("The fourth input of KepEq_f must be a real scalar!");
    mu = mxGetScalar(prhs[3]);
    
    if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])
                             || mxGetM(prhs[4])!=1
                             || mxGetN(prhs[4])!=1
                             || mxGetNumberOfDimensions(prhs[4]) != 2)
        mexErrMsgTxt("The fifth input of KepEq_f must be a real scalar!");
    f0 = mxGetScalar(prhs[4]);
    
    if (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5])
                             || mxGetM(prhs[5])!=1
                             || mxGetN(prhs[5])!=1
                             || mxGetNumberOfDimensions(prhs[5]) != 2)
        mexErrMsgTxt("The sixth input of KepEq_f must be a real scalar!");
    t0 = mxGetScalar(prhs[5]);
    
    /* Optional inputs */
    if(nrhs < 7)
        option = 0;
    else
    {
        if (!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6])
                                 || mxGetN(prhs[6])*mxGetM(prhs[6]) != 1
                                 || mxGetNumberOfDimensions(prhs[6]) != 2)
            mexErrMsgTxt("The seventh input of KepEq_f must be an integer scalar.");
        option = roundd2i(mxGetScalar(prhs[6]));
    } /* if(nrhs < 7) */
    
    if(nrhs < 8)
        imax = 5;
    else
    {
        if (!mxIsDouble(prhs[7]) || mxIsComplex(prhs[7])
                                 || mxGetN(prhs[7])*mxGetM(prhs[7]) != 1
                                 || mxGetNumberOfDimensions(prhs[7]) != 2)
            mexErrMsgTxt("The eighth input of KepEq_f must be an integer scalar.");
        imax = roundd2i(mxGetScalar(prhs[7]));
    } /* if(nrhs < 8) */
    
    if(nrhs < 9)
        tol = 1e-15;
    else
    {
        if (!mxIsDouble(prhs[8]) || mxIsComplex(prhs[8])
                                 || mxGetM(prhs[8])!=1
                                 || mxGetN(prhs[8])!=1
                                 || mxGetNumberOfDimensions(prhs[8]) != 2)
            mexErrMsgTxt("The ninth input of KepEq_f must be a real scalar!");
        tol = mxGetScalar(prhs[8]);
    } /* if(nrhs < 9) */
        
    /* Compute the outputs */
    error = kepEq_f(t, a, ecc, mu, f0, t0, imax, tol, &out);
    
    
    /* Allocate memory for the output(s), and assign value(s) */
   if(error!=2)
        plhs[0] = mxCreateDoubleScalar(out);
    else
        plhs[0] = mxCreateNumericArray(0, 0, mxINT16_CLASS, mxREAL);
    
    if(nlhs > 1)    /* error flag is required */
        plhs[1] = mxCreateDoubleScalar(error);
    
    /* Print if option == 1 */
    if(option==1)
    {
        if(error==1)
            mexPrintf("The Newton loop did not converge within the required tolerance");
        else if(error==2)
            mexPrintf("The orbit is not elliptic");
    }
        
    /* Exit the function */
    return;
}
