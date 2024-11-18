/******************************************************************************
*                               ephNEO.c                                      *
*              Functions for ephemerides of NEOs read from a file             *
*                                                                             *
*                                                       Matteo Ceriotti, 2009 *
*******************************************************************************/

#include "ephNEO.h"

/* Global variables */
static bodyDatabase bodyDB = {0};
static int dbLoaded = 0;

int loadBodyDatabase(const char *filename){
    int i;
    FILE *fp;
    char dummychar;
    int dummyint;
    
    if (dbLoaded){       /* Database already loaded */
#ifdef DEBUG
        printf("Database already loaded. Unload first.\n", filename);
#endif
        return 2;
    }
    
#ifdef DEBUG
    printf("Loading database %s.\n", filename);
#endif

    /* Load from file */
    if((fp = fopen(filename, "r")) == NULL) {
#ifdef DEBUG
        printf("Cannot open source file.\n");
#endif
        return 3;
    }
    
    fscanf(fp, "%d\n", &bodyDB.nBodies);
    
#ifdef DEBUG
    printf("Found %d bodies in database %s.\n", bodyDB.nBodies, filename);
#endif

    if (!(bodyDB.body = (bodyEphemeris *) malloc(bodyDB.nBodies*sizeof(bodyEphemeris)))){   /* Memory full */
        printf("Memory full.\n");
        fclose(fp);
        bodyDB.nBodies = 0;
        return 1;
    }
    
    dbLoaded = 1;

#ifdef MEXCOMPILE
    /* Makes memory allocated through mxMalloc persistent */
    mexMakeMemoryPersistent(bodyDB.body);

    /* Registers the function to clear the memory when MEX is cleared in MATLAB */
    mexAtExit(unloadBodyDatabase_void);
#endif

    for(i=0; i<bodyDB.nBodies; i++){
#ifdef DEBUG
        printf("Loading body %d.\n", i);
#endif
        fscanf(fp, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s",
            &dummyint,
            &bodyDB.body[i].kep[0],
            &bodyDB.body[i].kep[1],
            &bodyDB.body[i].kep[2],
            &bodyDB.body[i].kep[3],
            &bodyDB.body[i].kep[4],
            &bodyDB.body[i].kep[5],
            &bodyDB.body[i].epoch,
            &bodyDB.body[i].magnitude,
            &bodyDB.body[i].mass,
            &bodyDB.body[i].date,
            bodyDB.body[i].name);
#ifdef DEBUG
        printf("%d %f %f %f %f %f %f %f %f %f %f %s\n",
            dummyint,
            bodyDB.body[i].kep[0],
            bodyDB.body[i].kep[1],
            bodyDB.body[i].kep[2],
            bodyDB.body[i].kep[3],
            bodyDB.body[i].kep[4],
            bodyDB.body[i].kep[5],
            bodyDB.body[i].epoch,
            bodyDB.body[i].magnitude,
            bodyDB.body[i].mass,
            bodyDB.body[i].date,
            bodyDB.body[i].name);
#endif
        /* Keeps reading until EOL */
        do{
            fscanf(fp, "%c", &dummychar);
        }
        while(dummychar!='\n');

        /* printf("%f %f %f %f %f %f %f\n",  bodyDB.body[i].kep[0],bodyDB.body[i].kep[1],bodyDB.body[i].kep[2],bodyDB.body[i].kep[3],bodyDB.body[i].kep[4],bodyDB.body[i].kep[5],bodyDB.body[i].epoch);
        */
    }
    
    fclose(fp);
    
    return 0;
}

int unloadBodyDatabase(){
    if(!dbLoaded){
#ifdef DEBUG
        printf("No database loaded\n");
#endif
        return 1;
    }
        
#ifdef DEBUG
    printf("Unloading database\n");
#endif

    unloadBodyDatabase_void();
    return 0;
}

void unloadBodyDatabase_void(){
    if(bodyDB.body) free(bodyDB.body);
    bodyDB.body = NULL;
    bodyDB.nBodies = 0;
    dbLoaded = 0;
}

int ephNEO(const double mjd2000, const int id, double kep[], double *mass, double *meanAnom, char *name, double *date)
{
    /* Declaration */
	int id_neo;            /* Local value of the identification index of the NEO */
	double a;
	double e;
	double di;
	double omm;
	double om;
	double m;
	double n;
	double t0;
	double t;
	double p;
	double np;
   	double phi;
    double wom; /* In radians */
    double epoch_mjd2000;
    double d;
    int error = 0;      /* returned integer: 2 anomaly error, 1 if invalid id, 0 otherwise */
    
    /* Update input id */
    id_neo = id - EPHNEO_FIRSTNEOID;
    
    /* Check if database loaded */
    if (!dbLoaded){
        /* Database not loaded: load default filename */

#ifdef MEXCOMPILE
        mexWarnMsgIdAndTxt("spaceToolboxC:ephNEO:loadDefaultDatabase","No database loaded: loading default database");
#endif
#ifdef DEBUG
        printf("No database loaded. Loading default database.\n");
#endif

        if (loadBodyDatabase(EPHNEO_DEFAULTFILENAME)){
            /* Cannot load database */
#ifdef DEBUG
            printf("Error loading database.\n");
#endif
            return 2;
        }
    }
    
    /* Check if id_neo is valid */
	if((id_neo > bodyDB.nBodies-1) || id_neo < 0){
#ifdef DEBUG
        printf("id out of bounds.\n");
#endif
        return 1;
	}
    
    /* Reads the NEO database */
    a       = bodyDB.body[id_neo].kep[0] * AU;  /* from [AU] to [km]. Cf. astro_constants.h */
    e       = bodyDB.body[id_neo].kep[1];
    di      = bodyDB.body[id_neo].kep[2];
    omm     = bodyDB.body[id_neo].kep[3];
    om      = bodyDB.body[id_neo].kep[4];
    m       = bodyDB.body[id_neo].kep[5];
    n       = sqrt(MU_SUN / DCUBE(a)); /* MU_SUN defined in astro_constants.h */

    t0      = m * DEG2RAD / n;

    epoch_mjd2000 = bodyDB.body[id_neo].epoch - 51544.5;   /* Convert epoch from MJD to MJD2000 */
    t        = (mjd2000 - epoch_mjd2000) * 86400.0 + t0;

    p   = PI2 * sqrt(DCUBE(a)/MU_SUN);
    np  = floor(t/p);
    t   = t - p * np;
    
    (*meanAnom)   = n*t;
    if((*meanAnom)>PI)
        (*meanAnom) -= PI2;
    

    if(meanAn2eccAn(*meanAnom, e, &phi))
        return 3;

    if(eccAn2trueAn(phi, e, &wom))
        return 3;

    
    kep[0] = a;
    kep[1] = e;
    kep[2] = di * DEG2RAD;
    kep[3] = omm * DEG2RAD;
    kep[4] = om * DEG2RAD;
    kep[5] = wom;

    /* Mass */
    if(bodyDB.body[id_neo].mass != 0){
        (*mass) = bodyDB.body[id_neo].mass;
#ifdef DEBUG
        printf("Mass=%f\n", *mass);
#endif
    } else {
#ifdef DEBUG
        printf("bodyDB.body[id_neo].magnitude=%f\n", bodyDB.body[id_neo].magnitude);
#endif
        /* Polynomial fit to the H-d diagram which gives the diameter of the NEO
         * as a function of its magnitude.
         * Polynomial of fifth order, Note: NOT VERY ACCURATE! */
        d = - 2.522e-2 * pow(bodyDB.body[id_neo].magnitude,5)
            + 3.2961   * pow(bodyDB.body[id_neo].magnitude,4)
            - 1.7249e2 * DCUBE(bodyDB.body[id_neo].magnitude)
            + 4.5231e3 * DSQR(bodyDB.body[id_neo].magnitude)
            - 5.9509e4 * bodyDB.body[id_neo].magnitude
            + 3.1479e5;
#ifdef DEBUG
        printf("d=%f\n", d);
#endif        /* Estimated mass of the NEO computed considering a density of
         * 2 kg/dm^3 */
        (*mass) = 4.0 * PI / 3.0 * DCUBE(0.5 * d) * 2.0 * 1000.0;
#ifdef DEBUG
        printf("Mass est with d=%f\n", *mass);
#endif
    }
    
    strcpy(name, bodyDB.body[id_neo].name);
    
    /* Date */
    (*date) = bodyDB.body[id_neo].date;
    
    return error;
}
