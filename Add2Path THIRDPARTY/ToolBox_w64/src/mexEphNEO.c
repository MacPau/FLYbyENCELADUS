/******************************************************************************
 *                    Mex Gateway routine for ephNEO                          *
 *                                                                            *
 *                                                      Matteo Ceriotti, 2009 *
 ******************************************************************************/

#include "mexEphNEO.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double time;   /* Date at which the ephemerides shall be computed (input) */
    int id;        /* Identifier of the asteroid (input) */
    double *kep;   /* Arrays to hold the keplerian elements */
    double *pmass, *pmeanAnom;  /* Pointer to mass and and mean anomaly */
    double mass, meanAnom;      /* mass and and mean anomaly */
    double date;                /* Date of the ephemeris */
    char name[81];              /* Name of the body */
    char action[81];            /* String to identify action: load or unload */
    char filename[257];         /* Name of ephemeris file */
    int error;
    
    /* Check the inputs and outputs */
    if(nrhs > 2)
        mexErrMsgTxt("ephNEO mex requires up to 2 inputs.");
    
    if(nrhs < 1)
        mexErrMsgTxt("ephNEO mex requires at least 1 input.");
    
    if(mxIsChar(prhs[0])){
        /* Either load or unload ephemeris data */
        
        if(nlhs > 0)
            mexErrMsgTxt("ephNEO mex: Too many outputs. When loading or unloading, ephNEO does not return any output.");
        
        mxGetString(prhs[0], action, mxGetNumberOfElements(prhs[0])+1);

        if(!strcmp(action, "load")){
            /* Load */
            if(nrhs > 1){
                /* Filename provided */
                if(!mxIsChar(prhs[1])){
                    mexErrMsgTxt("ephNEO mex: When loading, second input must be filename.");
                }
                
                /* Loads ephemeris data */            
                if(mxGetString(prhs[1], filename, mxGetNumberOfElements(prhs[1])+1)){
                    mexErrMsgTxt("ephNEO mex: Error getting filename string.");
                }
                
                if (loadBodyDatabase(filename)){
                    /* Cannot load database */
                    mexErrMsgTxt("ephNEO mex: Error loading database.");
                }
            } else {
                /* Filename not provided, load default */
                if (loadBodyDatabase(EPHNEO_DEFAULTFILENAME)){
                    /* Cannot load database */
                    mexErrMsgTxt("ephNEO mex: Error loading default database.");
                }
            }
            
        } else if(!strcmp(action, "unload")){
            /* Unload */
            if(nrhs > 1)
                mexErrMsgTxt("ephNEO mex: when unloading, no additional inputs are possible.");
            
            /* Unload ephemeris data */
            if(unloadBodyDatabase()){
                mexWarnMsgIdAndTxt("spaceToolboxC:ephNEO:noDatabaseLoaded","No database loaded to unload.");
            }
        } else {
            mexErrMsgTxt("ephNEO mex: First string must be either \"load\" or \"unload\".");
        }
    } else {
        /* Normal ephemeris request */
    
        if(mxGetM(prhs[0])!=1 || mxGetN(prhs[0])!=1 || mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1)
            mexErrMsgTxt("The 2 inputs of ephNEO must be scalars");
        if(nlhs > 3)
            mexErrMsgTxt("ephNEO mex: Too many outputs. ephNEO can return up to 3 outputs.");

        /* Get values of inputs */
        if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                                 || mxGetN(prhs[0])*mxGetM(prhs[0]) != 1
                                 || mxGetNumberOfDimensions(prhs[0]) != 2)
            mexErrMsgTxt("ephNEO mex: The first input (time) must be a scalar.");
        time = mxGetScalar(prhs[0]);

        if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
                                 || mxGetN(prhs[1])*mxGetM(prhs[1]) != 1
                                 || mxGetNumberOfDimensions(prhs[1]) != 2)
            mexErrMsgTxt("ephNEO mex: The second input (id) must be an integer scalar.");
        id = roundd2i(mxGetScalar(prhs[1]));

        /*Allocate memory and assign output pointer */
        plhs[0] = mxCreateDoubleMatrix(1,6,mxREAL);

        /* Pointer to first output */
        kep = mxGetPr(plhs[0]);

        /* Compute the outputs */
        error = ephNEO(time, id, kep, &mass, &meanAnom, name, &date);

        if (error == 1)
            mexErrMsgTxt("ephNEO mex: The second input (id) is out of bound.");
        else if (error == 2)
            mexErrMsgTxt("ephNEO mex: Error loading database.");
        else if (error == 3)
            mexErrMsgTxt("ephNEO mex: Error computing true anomaly.");
        else if (error > 3)
            mexErrMsgTxt("ephNEO mex: Unknown error.");

        /* Assign the extra outputs if necessary */
        if(nlhs>1)
        {
            /*Allocate memory and assign output pointer */
            plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
            /* Pointer to the second output */
            pmass  = mxGetPr(plhs[1]);
            *pmass = mass;

            if(nlhs>2)
            {
                /*Allocate memory and assign output pointer */
                plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
                /* Pointer to the third output */
                pmeanAnom  = mxGetPr(plhs[2]);
                *pmeanAnom = meanAnom;
            }
        }
    }
    
    /* Exit the function */
    return;
}
