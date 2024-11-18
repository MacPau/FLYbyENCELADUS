/******************************************************************************
 *                  Mex Gateway routine for lambertMR                         *
 *                                                                            *
 *                                                 Matteo Ceriotti, 2008      *
 ******************************************************************************/

#include "mexLambertMR.h"

/* Prototype of LambertMR
int lambertMR(const double ri[3], const double rf[3], double tof, double mu, int orbittype, int nrev, int ncase, int optionslmr,
	double *a_out, double *p_out, double *e_out, int *error_out, double vi_out[3], double vf_out[3], double *tpar_out, double *theta_out);
*/

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double ri[3], rf[3];
    double tof, mu;
	int orbittype, nrev, ncase, optionslmr;
	double a_out, p_out, e_out, vi_out[3], vf_out[3], tpar_out, theta_out;
	int error_out;
	
	/* Check for proper number of arguments. */
    /* NOTE: You do not need an else statement when using 
        mexErrMsgTxt within an if statement. It will never 
        get to the else statement if mexErrMsgTxt is executed. 
        (mexErrMsgTxt breaks you out of the MEX-file.) 
    */
    
	if (nlhs > 8)
		mexErrMsgTxt("Too many output arguments. 8 is the maximum number.");
	if (nrhs < 4)
		mexErrMsgTxt("LambertMR_mex: at least 4 inputs required.");

	if (nrhs < 8)
		optionslmr = 0;
    else{
    	if (!mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) ||
        mxGetN(prhs[7])*mxGetM(prhs[7]) != 1 || mxGetNumberOfDimensions(prhs[7]) != 2)
            mexErrMsgTxt("Input 8 must be an integer scalar.");
        optionslmr = (mxGetScalar(prhs[7]) > 0) ? (int)(mxGetScalar(prhs[7]) + 0.5) : (int)(mxGetScalar(prhs[7]) - 0.5);
	}
	
	if(nrhs < 7)
    	ncase = 0;
    else{
    	if (!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) ||
        mxGetN(prhs[6])*mxGetM(prhs[6]) != 1 || mxGetNumberOfDimensions(prhs[6]) != 2)
            mexErrMsgTxt("Input 7 must be an integer scalar.");
        ncase = (mxGetScalar(prhs[6]) > 0) ? (int)(mxGetScalar(prhs[6]) + 0.5) : (int)(mxGetScalar(prhs[6]) - 0.5);
	}
	
	if (nrhs < 6)
		nrev = 0;
	else{
    	if (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) ||
        mxGetN(prhs[5])*mxGetM(prhs[5]) != 1 || mxGetNumberOfDimensions(prhs[5]) != 2)
            mexErrMsgTxt("Input 6 must be an integer scalar.");
        nrev = (mxGetScalar(prhs[5]) > 0) ? (int)(mxGetScalar(prhs[5]) + 0.5) : (int)(mxGetScalar(prhs[5]) - 0.5);
	}
	if (nrhs < 5)
		orbittype = 0;
	else{
    	if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) ||
        mxGetN(prhs[4])*mxGetM(prhs[4]) != 1 || mxGetNumberOfDimensions(prhs[4]) != 2)
            mexErrMsgTxt("Input 5 must be an integer scalar.");
        orbittype = (mxGetScalar(prhs[4]) > 0) ? (int)(mxGetScalar(prhs[4]) + 0.5) : (int)(mxGetScalar(prhs[4]) - 0.5);
	}
    
    /* Check to make sure the input arguments are correct and read. */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
        mxGetN(prhs[0])*mxGetM(prhs[0]) != 3 || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("Input 1 must be a real 1x3 vector.");
	cpdouble(mxGetPr(prhs[0]), ri, 3);
	
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
		mxGetN(prhs[1])*mxGetM(prhs[1]) != 3 || mxGetNumberOfDimensions(prhs[1]) != 2)
		mexErrMsgTxt("Input 2 must be a real 1x3 vector.");
	cpdouble(mxGetPr(prhs[1]), rf, 3);
	
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
        mxGetN(prhs[2])*mxGetM(prhs[2]) != 1 || mxGetNumberOfDimensions(prhs[2]) != 2)
		mexErrMsgTxt("Input 3 must be a real scalar.");
	tof = (double)mxGetScalar(prhs[2]);
	
    /* int n */
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
        mxGetN(prhs[3])*mxGetM(prhs[3]) != 1 || mxGetNumberOfDimensions(prhs[3]) != 2)
    	mexErrMsgTxt("Input 4 must be a real scalar.");
    mu = (double)mxGetScalar(prhs[3]);
	
#ifdef DEBUG
    mexPrintf("ri=%.14f, %.14f, %.14f\nrf=%.14f, %.14f, %.14f\ntof=%.14f\nmu=%f\npars=%d %d %d\n",ri[0],ri[1],ri[2],rf[0],rf[1],rf[2],tof,mu,orbittype,nrev,ncase);
#endif

    /* Call the C subroutine. */
    if(lambertMR(ri, rf, tof, mu, orbittype, nrev, ncase, optionslmr,
		&a_out, &p_out, &e_out, &error_out, vi_out, vf_out, &tpar_out, &theta_out))
    {
        ;/*printf("Error in exposin2!\n");*/
    }

#ifdef DEBUG
    mexPrintf("a=%.14f\np=%.14f\ne=%.14f\nerror=%i\nvi=%.14f %.14f %.14f\nvf=%.14f %.14f %.14f\ntpar=%.14f\ntheta=%.14f\n",a_out, p_out, e_out, error_out, vi_out[0], vi_out[1], vi_out[2], vf_out[0], vf_out[1], vf_out[2], tpar_out, theta_out);
#endif

	/* Set output to zeros if error, to be consistent with MATLAB version */
	if(error_out){
		a_out = 0; p_out = 0; e_out = 0;
		vi_out[0] = 0; vi_out[1] = 0; vi_out[2] = 0;
		vf_out[0] = 0; vf_out[1] = 0; vf_out[2] = 0;
		tpar_out = 0; theta_out = 0;
	}
	
	plhs[0] = mxCreateDoubleScalar(a_out); /* a */
	
    /* Assigne the extra outputs if necessary */
    if(nlhs>1){
        plhs[1] = mxCreateDoubleScalar(p_out); /* p */
        if(nlhs>2){
			plhs[2] = mxCreateDoubleScalar(e_out); /* e */
			if(nlhs>3){
                plhs[3] = mxCreateDoubleScalar(error_out); /* error */
                if(nlhs>4){
					plhs[4] = mxCreateDoubleMatrix(1, 3, mxREAL); /* *vi */
					cpdouble(vi_out, mxGetPr(plhs[4]), 3);
					if(nlhs>5){
                        plhs[5] = mxCreateDoubleMatrix(1, 3, mxREAL); /* *vf */
						cpdouble(vf_out, mxGetPr(plhs[5]), 3);
						if(nlhs>6){
							plhs[6] = mxCreateDoubleScalar(tpar_out); /* tpar */
							if(nlhs>7){
								plhs[7] = mxCreateDoubleScalar(theta_out); /* theta */
							}
						}
					}
				}
			}
		}
	}

	return;
}

void cpdouble(const double *in, double *out, int n){
	int i;
	for(i = 0; i < n; i++){
		out[i] = in[i];
	}
	return;
}
