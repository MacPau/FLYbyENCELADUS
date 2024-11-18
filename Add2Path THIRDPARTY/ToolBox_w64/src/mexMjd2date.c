/******************************************************************************
 *                      Mex Gateway routine for mjd2date                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexMjd2date.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double mjd;         /* Date expressed in MJD (input) */
    double *date;       /* 6-element array of double containing the date
                           converted in [Y M D h m s] (output) */
    int error;          /* error code */

    
    /* Check the number of inputs and outputs */
    if(nrhs != 1)
        mexErrMsgTxt("MJD2DATE requires 1 and only 1 inputs.");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. MJD2DATE returns only 1 output.");
    
    
    /* Get values of inputs */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetM(prhs[0])!=1
                             || mxGetN(prhs[0])!=1
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The only input of MJD2DATE must be a real scalar!");
    mjd = *mxGetPr(prhs[0]);
     
    
    /*Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(1,6,mxREAL);

    /* Pointer to first output */
    date = mxGetPr(plhs[0]);

    /* Compute the output */
    error = mjd2date(mjd,date);
    
    
    /* Check for error */
    if (error>0)
        mexErrMsgTxt("MJD2DATE: The function is valid for dates after since 12:00 noon 24 November -4713, Gregorian calendar, i.e. for MJD>-2400000.5 dates.");
    
    
    /* Exit the function */
    return;
}
