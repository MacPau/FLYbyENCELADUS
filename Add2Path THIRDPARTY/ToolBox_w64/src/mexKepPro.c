/******************************************************************************
 *                  Mex Gateway routine for kepPro3                           *
 *                                                                            *
 *                                                 Matteo Ceriotti, 2008      *
 ******************************************************************************/

#include "mexKepPro.h"

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Declarations */
	double *in, dt, mu, tol; /* Inputs */
	int maxiter, npar; /* Inputs */
   	double *out; /* Outputs */
	int error; /* Error code */

	/* Code */
	
	/* Check for proper number of arguments. */
    /* NOTE: You do not need an else statement when using 
        mexErrMsgTxt within an if statement. It will never 
        get to the else statement if mexErrMsgTxt is executed. 
        (mexErrMsgTxt breaks you out of the MEX-file.) 
    */
    
	if (nlhs > 2)
		mexErrMsgTxt("kepPro3_mex: Too many output arguments. 2 is the maximum number.");
	if (nrhs < 3)
		mexErrMsgTxt("kepPro3_mex: at least 3 inputs required.");
	if (nrhs > 5)
		mexErrMsgTxt("kepPro3_mex: Too many input arguments. 5 is the maximum number.");
	
	/* Check size of inputs and read them */
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||	mxGetN(prhs[0])*mxGetM(prhs[0]) != 6 || mxGetNumberOfDimensions(prhs[0]) != 2)
		mexErrMsgTxt("Input 1 (in) must be a real 1x6 vector.");
	in = (double *)mxGetPr(prhs[0]);
	
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||	mxGetN(prhs[1])*mxGetM(prhs[1]) != 1 || mxGetNumberOfDimensions(prhs[1]) != 2)
		mexErrMsgTxt("Input 2 (dt) must be a real scalar.");
	dt = (double)mxGetScalar(prhs[1]);
	
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetN(prhs[2])*mxGetM(prhs[2]) != 1 || mxGetNumberOfDimensions(prhs[2]) != 2)
		mexErrMsgTxt("Input 3 (mu) must be a real scalar.");
	mu = (double)mxGetScalar(prhs[2]);

	if(nrhs > 3){ /* tol provided */
		if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||	mxGetN(prhs[3])*mxGetM(prhs[3]) != 1 || mxGetNumberOfDimensions(prhs[3]) != 2)
			mexErrMsgTxt("Input 4 (tol) must be an integer scalar.");
		tol = (double)mxGetScalar(prhs[3]);
	}
    
	if (nrhs > 4){ /* maxiter provided */
		if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetN(prhs[4])*mxGetM(prhs[4]) != 1 || mxGetNumberOfDimensions(prhs[4]) != 2)
			mexErrMsgTxt("Input 5 (maxiter) must be an integer scalar.");
		maxiter = (mxGetScalar(prhs[4]) > 0) ? (int)(mxGetScalar(prhs[4]) + 0.5) : (int)(mxGetScalar(prhs[4]) - 0.5); /* Rounding to int */
	}
	
	/* Allocates memory for C output */
	plhs[0] = mxCreateDoubleMatrix(1, 6, mxREAL); /* out[6] */
	
	out = mxGetPr(plhs[0]);
	
    /* Call the C subroutine. */
    npar =  nrhs - 3; /* Number of optional arguments provided */
	if(error = kepPro(out, in, dt, mu, npar, tol, maxiter)){
		; /* Here something can be done in case of error */
	}
	
	if (nlhs > 1){ /* error requested */
        plhs[1] = mxCreateDoubleScalar((double)error); /* error */
    }
    
	return;
}
