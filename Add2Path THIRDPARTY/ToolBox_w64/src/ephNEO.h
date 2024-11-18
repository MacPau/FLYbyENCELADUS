/******************************************************************************
*                               ephNEO.h                                      *
*              Functions for ephemerides of NEOs read from a file             *
*                                                                             *
*                                                       Matteo Ceriotti, 2009 *
*******************************************************************************/

#ifndef EPHNEO_H
#define EPHNEO_H

/* Include libraries
 * ----------------- */
#ifdef MEXCOMPILE
#include "mex.h"
#endif

#include <string.h>

#include "mathUtils.h"         /* In-house C mathematical utilities */
#include "astroConstants.h"    /* In-house astronautical constants */
#include "anomConv.h"    /* In-house conversion between anomalies */

/* Mex compiler switch
 * ------------------- */
#ifdef MEXCOMPILE
#define printf mexPrintf
#define malloc mxMalloc
#define calloc mxCalloc
#define free mxFree
#endif

#ifndef MEXCOMPILE
#include <stdio.h>
#endif

/* #define DEBUG */

/* Definitions */
#define EPHNEO_DEFAULTFILENAME "ephNEO_data.txt" /* Default filename for body database */
#define EPHNEO_FIRSTNEOID 12 /* ID of the first NEO */

/* Constants 
 * --------- */

/* Data types */
typedef struct {
    double kep[6];
    double epoch;
    double magnitude;
    double mass;
    double date;
    char name[81];
} bodyEphemeris;

typedef struct {
    int nBodies;
    bodyEphemeris *body;
} bodyDatabase;

/* Prototype
 * --------- */

int ephNEO(const double mjd2000, const int id, double kep[], double *mass, double *meanAnom, char *name, double *date);
/* ephNEO - Ephemerides of NEOs read from a text file database.
* 
* PROTOTYPE:
*   int ephNEO(const double mjd2000, const int id,
*       double kep[], double *mass, double *meanAnom,
*       char *name, double *date);
* 
* DESCRIPTION:
*   This function returns the orbital parameters, the mass, the mean
*   anomaly and the name of NEOs. Each NEO is identified by an id.
*   NEO data is read from a database of ephemerides in text file, and stored
*   in the structure bodyDatabase.
*   The database is read the first time the function is invoked and then
*   stored into persistent memory.
*   If no database is loaded, it invokes the loadBodyDatabase function with
*   the default filename EPHNEO_DEFAULTFILENAME.
*   The ID of the first NEO in the database is conventionally defined by
*   EPHNEO_FIRSTNEOID.
*
*   Global variables:
*       static bodyDatabase bodyDB;
*       static int dbLoaded;
*   must be defined.
*
* INPUT:
*   (double) mjd2000        Epoch in MJD2000
*                           (Modified Julian Day since 01/01/2000, 12:00
*                           MJD2000 = MJD-51544.5)
*   (int) id                NEO identifier. It is a number starting from
*                           EPHNEO_FIRSTNEOID, because lower ids are reserved
*
* OUTPUT:
*   (double) kep[6]         Pointer to a pre-allocated vector of 6 doubles
*                           containing the keplerian elements at epoch.
*                               kep[0] = semimajor axis in km
*                               kep[1] = eccentricity
*                               kep[2] = inclination in rad
*                               kep[3] = right ascension of the
*                                   ascending node in rad
*                               kep[4] = argument of perigee in rad
*                               kep[5] = true anomaly in rad
*   (double *) mass         Pointer to a double containing the mass of
*                           the NEO [kg]. It can be read from the
*                           database, or, if not available, estimated by
*                           an approximate equation based on the magnitude.
*   (double *) meanAnom     Pointer to a double containing the mean
*                           anomaly at time [rad].
*   (int) ephNEO            Error code:
*                               0: no error.
*                               1: identifier not valid.
*                               2: error loading database.
*                               3: error converting to true anomaly.
*     
* NON-STANDARD LIBRARIES:
*   mathUtils.h, astroConstants
*
* ORIGINAL VERSION:
*   Paolo de Pascale, 01/11/2004, MATLAB, NeoEphemeris.m
*   
* AUTHOR:
*   Matteo Ceriotti, 09/12/2009
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

int unloadBodyDatabase();
/* unloadBodyDatabase - Check and unload database of bodies.
*
* PROTOTYPE:
*   int unloadBodyDatabase()
* 
* DESCRIPTION:
*   Used by ephNEO.
*
*   Global variables:
*       static bodyDatabase bodyDB;
*       static int dbLoaded;
*   must be defined.
*
* INPUT:
*   (none)
*     
* OUTPUT:
*   (int) unloadBodyDatabase    Error code: 0 if database unloaded,
*                                   1 if no database was loaded (does nothing)
*
* AUTHOR:
*   Matteo Ceriotti, 09/12/2009
*
* -------------------------------------------------------------------------
*/

void unloadBodyDatabase_void();
/* unloadBodyDatabase - Unload database of bodies with no check.
*
* PROTOTYPE:
*   void unloadBodyDatabase_void()
* 
* DESCRIPTION:
*   Used by unloadBodyDatabase.
*
*   Global variables:
*       static bodyDatabase bodyDB;
*       static int dbLoaded;
*   must be defined.
*
* INPUT:
*   (none)
*     
* OUTPUT:
*   (none)
*
* AUTHOR:
*   Matteo Ceriotti, 09/12/2009
*
* -------------------------------------------------------------------------
*/

int loadBodyDatabase(const char *filename);
/* loadBodyDatabase - Load database of bodies from file
*
* PROTOTYPE:
*   int loadBodyDatabase(const char *filename)
* 
* DESCRIPTION:
*   Opens the specified file and populates the structure bodyDatabase with
*   ephemeris data.
*   Used by ephNEO.
*   The data file should be:
*       numberOfNEOs
*       n a e i Omega omega meanAnom epoch mag mass refDate name comment
*       n a e i Omega omega meanAnom epoch mag mass refDate name comment
*       n a e i Omega omega meanAnom epoch mag mass refDate name comment
*       ...
*   Where:
*       numberOfNEOs    Number of objects in database.
*       n               Id of NEO. The first ID should be matching
*                       EPHNEO_FIRSTNEOID. Not actually used.
*       a               Semimajor axis, AU.
*       e               Eccentricity.
*       i               Inclination, deg.
*       Omega           Right ascension of ascending node, deg.
*       omega           Anomaly of pericentre, deg.
*       meanAnom        Mean anomaly at epoch, deg.
*       epoch           Epoch of ephemeris, MJD.
*       mag             Magnitude. Used to compute the mass, if the one in
*                       the database is 0.
*       mass            Mass, kg. If 0, then mag is used to estimate the mass.
*       refDate         Date at which the ephemeris was added, or any other
*                       reference date.
*       name            Name (max 80 characters, no spaces allowed)
*       comment         Any comment, even with spaces.
*   Note:
*       All floating point numbers, except numberOfNEOs and n, must contain
*       the decimal separator ".": so zero must be "0.0", etc. Exponential 
*       notation is possible.
*   
*   Global variables:
*       static bodyDatabase bodyDB;
*       static int dbLoaded;
*   must be defined.
*
* INPUT:
*   (char *) filename           Database filename.
*     
* OUTPUT:
*   (int) loadBodyDatabase      Error code:
*                                   0: No error.
*                                   1: Cannot allocate dynamic memory (memory full?).
*                                   2: Database already loaded.
*                                   3: Cannot open database file.
*
* AUTHOR:
*   Matteo Ceriotti, 09/12/2009
*
* -------------------------------------------------------------------------
*/

 #endif /* EPHNEO_H */
