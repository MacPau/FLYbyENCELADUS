/******************************************************************************
 *                      Mex Gateway routine for kep2car                       *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexKep2car.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double *kep;        /* Arrays to hold the keplerian elements (inputs) */
    double mu;          /* planetary gravity constant (input) */
    double p;           /* Semi-latus rectum (input) */
    double *out;        /* array to hold the state vector in cartesian coordinates (output) */
    int orbitType;      /* orbit type returned by KEP2CART */
    mwSize nr, nc;      /* number of rows and columns of the first input */
    
    
    /* Check the number of inputs and outputs */
    if(nrhs > 3)
        mexErrMsgTxt("KEP2CAR requires up to 3 inputs!");
    if(nlhs > 2)
        mexErrMsgTxt("Too many outputs. KEP2CAR can return 2 outputs maximum!");
    if (nrhs < 2)
		mexErrMsgTxt("KEP2CAR requires at least 2 inputs (kep, mu)!");
    
    /* Get values of inputs */
    nr = mxGetM(prhs[0]);
    nc = mxGetN(prhs[0]);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || (!(nr==6 && nc==1) && !(nr==1 && nc==6))
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The first input of KEP2CAR must be a 6-elements vector of reals!");
    kep = mxGetPr(prhs[0]);
        
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
                             || mxGetM(prhs[1])*mxGetN(prhs[1]) != 1
                             || mxGetNumberOfDimensions(prhs[1]) != 2)
        mexErrMsgTxt("The second input of KEP2CAR (mu) must be a real scalar!");
    mu = mxGetScalar(prhs[1]);
    
    if (nrhs > 2){
            if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])
                             || mxGetM(prhs[2])*mxGetN(prhs[2]) != 1
                             || mxGetNumberOfDimensions(prhs[2]) != 2)
        mexErrMsgTxt("The third input of KEP2CAR (p) must be a real scalar!");
        p = mxGetScalar(prhs[2]);
    } else {
        p = 0; /* Not provided, 0 by default */
    }
    
    /*Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(nr,nc,mxREAL);

    /* Pointer to first output */
    out = mxGetPr(plhs[0]);
    
    /* Compute the outputs */
    orbitType = kep2car(kep, mu, p, out);

    if (nrhs < 3 && orbitType == 2) {
        /* This was a parabolic case, but p was not explicitly provided */
        mexErrMsgTxt("Parabolic case: the semi-latus rectum needs to be provided");
    }
    
    if(nlhs>1)      /* OrbitType requested */
        plhs[1] = mxCreateDoubleScalar((double) orbitType);
        
    /* Exit the function */
    return;
}
