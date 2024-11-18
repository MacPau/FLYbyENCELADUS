/******************************************************************************
 *                       Mex Gateway routine for jd2date                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexJd2date.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double jd;          /* Date expressed in JD (input) */
    double *date;       /* 6-element array of double containing the date
                           converted in [Y M D h m s] (output) */
    int error;          /* error code */

    
    /* Check the number of inputs and outputs */
    if(nrhs != 1)
        mexErrMsgTxt("JD2DATE requires 1 and only 1 inputs.");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. JD2DATE returns only 1 output.");
    
    
    /* Get values of inputs */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetM(prhs[0])!=1
                             || mxGetN(prhs[0])!=1
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The only input of JD2DATE must be a real scalar!");
    jd = *mxGetPr(prhs[0]);
     
    
    /*Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(1,6,mxREAL);

    /* Pointer to first output */
    date = mxGetPr(plhs[0]);

    /* Compute the output */
    error = jd2date(jd,date);
    
    
    /* Check for error */
    if (error>0)
        mexErrMsgTxt("JD2DATE: The function is valid for dates after since 12:00 noon 24 November -4713, Gregorian calendar, i.e. for positive JD dates.");
    
    
    /* Exit the function */
    return;
}
