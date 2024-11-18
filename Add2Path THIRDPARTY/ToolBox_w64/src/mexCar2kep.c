/******************************************************************************
 *                      Mex Gateway routine for car2kep                       *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexCar2kep.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double *in;         /* array to hold the state vector in cartesian coordinates (inputs) */
    double mu;          /* planetary gravity constant (input) */
    double *kep;        /* Arrays to hold the keplerian elements (outputs) */
    double p;           /* Parameter */
    double eccAnom;     /* Eccentric (or hyperbolic or parabolic anomaly) */
    double meanAnom;    /* Mean anomaly */
    double dt;          /* Time from the pericentre passage */
    int orbitType;      /* Orbit type */
    mwSize nr, nc;      /* number of rows and columns of the first input */
    
    
    /* Check the number of inputs and outputs */
    if(nrhs != 2)
        mexErrMsgTxt("CAR2KEP requires 2 and only 2 inputs!");
    if(nlhs > 6)
        mexErrMsgTxt("Too many outputs. CAR2KEP can return up to 6 outputs!");
    
    
    /* Get values of inputs */
    nr = mxGetM(prhs[0]);
    nc = mxGetN(prhs[0]);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || (!(nr==6 && nc==1) && !(nr==1 && nc==6))
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The first input of CAR2KEP must be a 6-elements vector of reals!");
    in = mxGetPr(prhs[0]);
        
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
                             || mxGetM(prhs[1])*mxGetN(prhs[1]) != 1
                             || mxGetNumberOfDimensions(prhs[1]) != 2)
        mexErrMsgTxt("The second input of CAR2KEP must be a real scalar!");
    mu = mxGetScalar(prhs[1]);
    
    
    /*Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(nr,nc,mxREAL);

    /* Pointer to first output */
    kep = mxGetPr(plhs[0]);
    
    /* Compute the outputs */
    orbitType = car2kep(in, mu, kep, &eccAnom, &meanAnom, &dt, &p);
    if(orbitType==2)
    {
        mexWarnMsgIdAndTxt("spaceToolboxC:car2kep:parabolicOrbit","CAR2KEP: Parabola. Semi-major axis is Inf.");
        kep[0] = mxGetInf();
    }
    
    /* Assign the extra outputs if necessary */
    if(nlhs>1)      /* Parameter requested */
    {
        plhs[1] = mxCreateDoubleScalar(p);
        
        if(nlhs>2)      /* Eccentric, hyperbolic or parabolic anomaly requested */
        {
            plhs[2] = mxCreateDoubleScalar(eccAnom);
            
            if(nlhs>3)      /* Mean anomaly requested */
            {
                plhs[3] = mxCreateDoubleScalar(meanAnom);
                
                if(nlhs>4)      /* Time from the pericentre passage requested */
                {
                    plhs[4] = mxCreateDoubleScalar(dt);
                    
                    if(nlhs>5)      /* OrbitType requested */
                        plhs[5] = mxCreateDoubleScalar((double) orbitType);
                }
            }
        }
    }

    /* Exit the function */
    return;
}
