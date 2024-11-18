/******************************************************************************
*                               astroConstants.c                              *
*                  Definition of the astronautical constants                  *
*																			  *
*								 Space Toolbox								  *
******************************************************************************/

#include "astroConstants.h"

int astroConstants(const int n, const int in[], double out[]){


    /* Declare variables */
    int i;       /* counter */
    int error=0; /* error output */
    
    /* For loop on the number of inputs */
    for(i=0;i<n;i++){
        /* Switch loop */
        switch (in[i])
        {
            case 1:
            {
              out[i] = GRAVITY_CONSTANT;
              break;
            }
            case 2:
            {
              out[i] = AU; /* From DITAN */
              break;
            }
            case 3:
            {
              out[i] = R_SUN; /* From DITAN */
              break;
            }
            case 4:
            {
              out[i] = MU_SUN; /* From DITAN */
              break;
            }
            case 5:
            {
              out[i] = V_LIGHT; /* Definition in the SI */
              break;
            }
            case 6:
            {
              out[i] = FREE_FALL; /* Definition in Wertz */
              break;
            }
            case 7:
            {
              out[i] = D_EARTH_MOON; /* Definition in Wertz */
              break;
            }
            case 8:
            {
              out[i] = OBLIQUITY_ECLIPTIC_2000; /* Definition in Wertz */
              break;
            }
            case 11:
            {
                out[i] = MU_ME; /* From DITAN */
                break;
            }
            case 12:
            {
                out[i] = MU_V; /* From DITAN */
                break;
            }
            case 13:
            {
                out[i] = MU_E; /* From DITAN */
                break;
            }
            case 14:
            {
                out[i] = MU_M; /* From DITAN */
                break;
            }
            case 15:
            {
                out[i] = MU_J; /* From DITAN */
                break;
            }
            case 16:
            {
                out[i] = MU_S; /* From DITAN */
                break;
            }
            case 17:
            {
                out[i] = MU_U; /* From DITAN */
                break;
            }
            case 18:
            {
                out[i] = MU_N; /* From DITAN */
                break;
            }
            case 19:
            {
                out[i] = MU_P; /* From DITAN */
                break;
            }
            case 20:
            {
                out[i] = MU_MOON; /* From DITAN */
                break;
            }
            case 21:
            {
                out[i] = R_ME; /* From DITAN */
                break;
            }
            case 22:
            {
                out[i] = R_V; /* From DITAN */
                break;
            }
            case 23:
            {
                out[i] = R_E; /* From DITAN */
                break;
            }
            case 24:
            {
                out[i] = R_M; /* From DITAN */
                break;
            }
            case 25:
            {
                out[i] = R_J; /* From DITAN */
                break;
            }
            case 26:
            {
                out[i] = R_S; /* From DITAN */
                break;
            }
            case 27:
            {
                out[i] = R_U; /* From DITAN */
                break;
            }
            case 28:
            {
                out[i] = R_N; /* From DITAN */
                break;
            }
            case 29:
            {
                out[i] = R_P; /* From DITAN */
                break;
            }
            case 30:
            {
                out[i] = R_MOON; /* From DITAN */
                break;
            }
            case 31:
            {
                out[i] = FLUX_DENSITY; /* From SMAD/Wertz */
                break;
            }
            /* Add an identifier and constant here. Prototype:
             * case $identifier$:
             * {
             *     out[i] = $constant_value$;
             *     break;
             * }                                                             */
            default:
            {
                printf("Constant identifier %d is not defined!\n",in[i]);
                out[i] = 0.0;
                error = 1;
                break;
            }
        }
    }
    return error;
}
