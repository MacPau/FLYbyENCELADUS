/******************************************************************************
*                                ephemerides.c                                *
*                             Ephemerides functions                           *
*                                                                             *
*                                Space Toolbox                                *
******************************************************************************/ 


#include "ephemerides.h"

int uplanet(const double mjd2000, const int ibody, double kep[])
{
    /* Declaration */
    double au = 149597870.66; /* Astronomical unit (not the one from astro_constants) */
    double t, tt, ttt;  /* Date, date^2 and date^3 */
    double xm;
    double phi;         /* eccentric anomaly */
    double g, g_primo;	/* variables to compute the eccentric anomaly */
    
    int i;              /* Counter */
    int error;          /* returned integer: 1 if ibody=11 (Moon),
                           2 if idbody>11, 0 otherwise */
    
    /*
     * t = JULIAN CENTURIES SINCE 31/12/1899 at 12:00
     *                                                                   */
    t   = (mjd2000 + 36525)/36525.00;
    tt  = t*t;
    ttt = t*tt;

    /*
     *  CLASSICAL PLANETARY ELEMENTS ESTIMATION IN MEAN ECLIPTIC OF DATE
     *                                                                   */
    switch (ibody)
    {

        /*
         *  MERCURY
         *                                                               */
        case 1:
        {
            kep[0] = 0.38709860;
            kep[1] = 0.205614210 + 0.000020460*t - 0.000000030*tt;
            kep[2] = 7.002880555555555560 + 1.86083333333333333e-3*t - 1.83333333333333333e-5*tt;
            kep[3] = 4.71459444444444444e+1 + 1.185208333333333330*t + 1.73888888888888889e-4*tt;
            kep[4] = 2.87537527777777778e+1 + 3.70280555555555556e-1*t +1.20833333333333333e-4*tt;
            xm     = 1.49472515288888889e+5 + 6.38888888888888889e-6*t;
            kep[5] = 1.02279380555555556e2 + xm*t;
            break;
        }
        /*
         *  VENUS
         *                                                               */
        case 2:
        {
            kep[0] = 0.72333160;
            kep[1] = 0.006820690 - 0.000047740*t + 0.0000000910*tt;
            kep[2] = 3.393630555555555560 + 1.00583333333333333e-3*t - 9.72222222222222222e-7*tt;
            kep[3] = 7.57796472222222222e+1 + 8.9985e-1*t + 4.1e-4*tt;
            kep[4] = 5.43841861111111111e+1 + 5.08186111111111111e-1*t -1.38638888888888889e-3*tt;
            xm     = 5.8517803875e+4 + 1.28605555555555556e-3*t;
            kep[5] = 2.12603219444444444e2 + xm*t;
            break;
        }
        /*
         *  EARTH
         *                                                               */
        case 3:
        {
            kep[0] = 1.000000230;
            kep[1] = 0.016751040 - 0.000041800*t - 0.0000001260*tt;
            kep[2] = 0.00;
            kep[3] = 0.00;
            kep[4] = 1.01220833333333333e+2 + 1.7191750*t + 4.52777777777777778e-4*tt + 3.33333333333333333e-6*ttt;
            xm     = 3.599904975e+4 - 1.50277777777777778e-4*t - 3.33333333333333333e-6*tt;
            kep[5] = 3.58475844444444444e2 + xm*t;
            break;
        }
        /*
         *  MARS
         *                                                               */
        case 4:
        {
            kep[0] = 1.5236883990;
            kep[1] = 0.093312900 + 0.0000920640*t - 0.0000000770*tt;
            kep[2] = 1.850333333333333330 - 6.75e-4*t + 1.26111111111111111e-5*tt;
            kep[3] = 4.87864416666666667e+1 + 7.70991666666666667e-1*t - 1.38888888888888889e-6*tt - 5.33333333333333333e-6*ttt;
            kep[4] = 2.85431761111111111e+2 + 1.069766666666666670*t +  1.3125e-4*tt + 4.13888888888888889e-6*ttt;
            xm     = 1.91398585e+4 + 1.80805555555555556e-4*t + 1.19444444444444444e-6*tt;
            kep[5] = 3.19529425e2 + xm*t;
            break;
        }
        /*
         *  JUPITER
         *                                                               */
        case 5:
        {
            kep[0] = 5.2025610;
            kep[1] = 0.048334750 + 0.000164180*t  - 0.00000046760*tt -0.00000000170*ttt;
            kep[2] = 1.308736111111111110 - 5.69611111111111111e-3*t +  3.88888888888888889e-6*tt;
            kep[3] = 9.94433861111111111e+1 + 1.010530*t + 3.52222222222222222e-4*tt - 8.51111111111111111e-6*ttt;
            kep[4] = 2.73277541666666667e+2 + 5.99431666666666667e-1*t + 7.0405e-4*tt + 5.07777777777777778e-6*ttt;
            xm     = 3.03469202388888889e+3 - 7.21588888888888889e-4*t + 1.78444444444444444e-6*tt;
            kep[5] = 2.25328327777777778e2 + xm*t;
            break;
        }
        /*
         *  SATURN
         *                                                               */
        case 6:
        {
            kep[0] = 9.5547470;
            kep[1] = 0.055892320 - 0.00034550*t - 0.0000007280*tt + 0.000000000740*ttt;
            kep[2] = 2.492519444444444440 - 3.91888888888888889e-3*t - 1.54888888888888889e-5*tt + 4.44444444444444444e-8*ttt;
            kep[3] = 1.12790388888888889e+2 + 8.73195138888888889e-1*t -1.52180555555555556e-4*tt - 5.30555555555555556e-6*ttt;
            kep[4] = 3.38307772222222222e+2 + 1.085220694444444440*t + 9.78541666666666667e-4*tt + 9.91666666666666667e-6*ttt;
            xm     = 1.22155146777777778e+3 - 5.01819444444444444e-4*t - 5.19444444444444444e-6*tt;
            kep[5] = 1.75466216666666667e2 + xm*t;
            break;
        }
        /*
         *  URANUS
         *                                                               */
        case 7:
        {
            kep[0] = 19.218140;
            kep[1] = 0.04634440 - 0.000026580*t + 0.0000000770*tt;
            kep[2] = 7.72463888888888889e-1 + 6.25277777777777778e-4*t + 3.95e-5*tt;
            kep[3] = 7.34770972222222222e+1 + 4.98667777777777778e-1*t + 1.31166666666666667e-3*tt;
            kep[4] = 9.80715527777777778e+1 + 9.85765e-1*t - 1.07447222222222222e-3*tt - 6.05555555555555556e-7*ttt;
            xm     = 4.28379113055555556e+2 + 7.88444444444444444e-5*t + 1.11111111111111111e-9*tt;
            kep[5] = 7.26488194444444444e1 + xm*t;
            break;
        }
        /*
         *  NEPTUNE
         *                                                               */
        case 8:
        {
            kep[0] = 30.109570;
            kep[1] = 0.008997040 + 0.0000063300*t - 0.0000000020*tt;
            kep[2] = 1.779241666666666670 - 9.54361111111111111e-3*t - 9.11111111111111111e-6*tt;
            kep[3] = 1.30681358333333333e+2 + 1.0989350*t + 2.49866666666666667e-4*tt - 4.71777777777777778e-6*ttt;
            kep[4] = 2.76045966666666667e+2 + 3.25639444444444444e-1*t + 1.4095e-4*tt + 4.11333333333333333e-6*ttt;
            xm     = 2.18461339722222222e+2 - 7.03333333333333333e-5*t;
            kep[5] = 3.77306694444444444e1 + xm*t;
            break;
        }
        /*
         *  PLUTO
         *                                                               */
        case 9:
        {
            kep[0] = 39.481686778174627;
            kep[1] = 2.4467e-001;
            kep[2] = 17.150918639446061;
            kep[3] = 110.27718682882954;
            kep[4] = 113.77222937912757;
            xm     = 4.5982945101558835e-008;
            kep[5] = 1.5021e+001 + xm*mjd2000*86400;
            break;
        }
        /*
         *  SUN
         *                                                               */
        case 10:
        {
            for (i=0;i<6;i++)
                kep[i] = 0.0;
            
            break;
        }
        /*
         *  MOON (AROUND EARTH)
         *                                                               */
    /*     case 11:
     *     {
     *         kep[0] = 0.0025695549067660;
     *         kep[1] = 0.0549004890;
     *         kep[2] = 5.145396388888888890;
     *         kep[3] = 2.59183275e+2  - 1.93414200833333333e+3*t + 2.07777777777777778e-3*tt + 2.22222222222222222e-6*ttt;
     *         kep[3] = mod(kep[3] + 108e3, 360);
     *         kep[4] = 7.51462805555555556e+1    + 6.00317604166666667e+3*t - 1.24027777777777778e-2*tt - 1.47222222222222222e-5*ttt;
     *         kep[4] = mod(kep[4], 360);
     *         xm     = 4.77198849108333333e+5    + 9.19166666666666667e-3*t  + 1.43888888888888889e-5*tt;
     *         kep[5] = 2.96104608333333333e2     + xm                    *t;
     *                                                                   */

        default:
        {
            if (ibody==11)
            {
                error = 1;
                goto exit;
            }
            else
            {
                error = 2;
                goto exit;
            }
        }
    }
    
    /*
     *  CONVERSION OF AU INTO KM, DEG INTO RAD AND DEFINITION OF  XMU
     *                                                                   */
 
    kep[0] = kep[0] * au;       /* a [km] */
    kep[2] = kep[2] * DEG2RAD;  /* Transform from deg to rad */
    kep[3] = kep[3] * DEG2RAD;  /* Transform from deg to rad */
    kep[4] = kep[4] * DEG2RAD;  /* Transform from deg to rad */
    kep[5] = kep[5] * DEG2RAD;  /* Transform from deg to rad */
    kep[5] = fmod(kep[5],PI2);  /* kep[5] within [0, 2*pi[ */
    /* XMU  = (xm*DEG2RAD/(864000*365250))^2*kep[0]^3; */
    phi  = kep[5];              /* phi is the eccentric anomaly */
                                /* uses kep[5]=M as a first guess */
     
    for (i=0;i<5;i++)
    {
        g       = kep[5]-(phi-kep[1]*sin(phi));
        g_primo = (-1.0 + kep[1]*cos(phi));
        phi     = phi-g/g_primo;   /* Computes the eccentric anomaly kep */
    }

    kep[5] = 2.0*atan(sqrt((1.0+kep[1])/(1.0-kep[1]))*tan(phi/2.0));

    error = 0;
    
exit:    
    return error;

}


int ephMoon_pos(const double mjd2000, double xp[])
{
    /* Declaration */
    double t_tdb;
    double angles[10];
    double L_ecl, phi_ecl, p, eps, r;
    double s[10],c[10];
    
    int i;      /* counter */
    
    /* Convert the date */
    t_tdb = mjd2000 / 36525.0;
    
    /* Compute the position vector at epoch date*/
    /*---------------------------------------------*/
    
    /* Definition of the angles */
    angles[0]  = 134.9 + 477198.85*t_tdb;   /* L_ecl 1 and p 1 */
    angles[1]  = 259.2 - 413335.38*t_tdb; 	/* L_ecl 2 and p 2 */
    angles[2]  = 235.7 + 890534.23*t_tdb; 	/* L_ecl 3 and p 3 */
    angles[3]  = 269.9 +  954397.7*t_tdb; 	/* L_ecl 4 and p 4 */
    angles[4]  = 357.5 +  35999.05*t_tdb; 	/* L_ecl 5 */
    angles[5]  = 186.6 + 966404.05*t_tdb; 	/* L_ecl 6 */
    angles[6]  =  93.3 + 483202.03*t_tdb; 	/* phi_ecl 1 */
    angles[7]  = 228.2 + 960400.87*t_tdb; 	/* phi_ecl 2 */
    angles[8]  = 318.3 +   6003.18*t_tdb; 	/* phi_ecl 3 */
    angles[9]  = 217.6 -  407332.2*t_tdb; 	/* phi_ecl 4 */
    
    /* Transform angles in radians */
    for(i=0;i<10;i++)
    {
        angles[i] = angles[i] * DEG2RAD;
        c[i] = cos(angles[i]);
        s[i] = sin(angles[i]);
    }
    
    
    /* Compute L_ecl */
    L_ecl =     218.32
           + 481267.883*t_tdb
           +      6.29*sin(angles[0])
           -      1.27*sin(angles[1])
           +      0.66*sin(angles[2])
           +      0.21*sin(angles[3])
           -      0.19*sin(angles[4])
           -      0.11*sin(angles[5]);
    
    /* Compute phi_ecl */
    phi_ecl =  5.13*sin(angles[6])
             + 0.28*sin(angles[7])
             - 0.28*sin(angles[8])
             - 0.17*sin(angles[9]);
    
    /* Compute p */
    p =  0.9508
       + 0.0518*cos(angles[0])
       + 0.0095*cos(angles[1])
       + 0.0078*cos(angles[2])
       + 0.0028*cos(angles[3]);
    
    /* Compute eps */
    eps =  23.439291
         -  0.0130042*t_tdb
         -  1.64e-7*DSQR(t_tdb)
         +  5.04e-7*DCUBE(t_tdb);
    
    /* Transform L_ecl, phi_ecl, p and eps in radians */
    L_ecl   = L_ecl * DEG2RAD;
    phi_ecl = phi_ecl * DEG2RAD;
    p       = p * DEG2RAD;
    eps     = eps * DEG2RAD;
    
    /* Compute r */
    r = 1/sin(p) * R_E;
    
    /* Compute the position vector xp */
    xp[0] = r * cos(L_ecl)*cos(phi_ecl);
    xp[1] = r * (cos(eps)*cos(phi_ecl)*sin(L_ecl) - sin(eps)*sin(phi_ecl));
    xp[2] = r * (sin(eps)*cos(phi_ecl)*sin(L_ecl) + cos(eps)*sin(phi_ecl));
    
    return 0;
}

int ephMoon(const double mjd2000, double xp[], double vp[])
{
    /* Declaration */
    double xp2[3];  /* Pointer to array of cartesian position of the Moon */
    double mjd2000_2;   /* Date to be 1 second after date to compute the
                           velocity */
    
    mjd2000_2 = mjd2000 + 1.0/86400.0;     /* date MJD2000 + 1 second */
    
    /* Compute the cartesian position of the Moon at date */
    ephMoon_pos(mjd2000, xp);
    
    /* Compute the cartesian position of the Moon at date + 1 second */
    ephMoon_pos(mjd2000_2, xp2);
    
    /* Compute the velocity of the Moon [km/s] */
    dmatdiff(xp2,xp,3,vp);
    
    return 0;
}

int ephSS_car(const int ibody, const double mjd2000, double r[], double v[])
{
    /* Declaration */
    double kep[6];          /* Keplerian elements */
    double car[6];          /* Cartesian coordinates */
    double mass, meanAnom, date;  /* extra parameters needed to call ephNEO */
    int error = 0;          /* returned integer: 0 if no error
                               1 if the error occurs in UPLANET
                               2 if the error occurs in NEOEPHEMERIS */
	char name[81];
    
    /* Call the appropriate ephemeris */
    if(ibody<11)            /* uplanet needed */
    {
        if (uplanet(mjd2000, ibody, kep)!=0)
        {
            error = 1;
            goto exit;
        }
    }
    else if(ibody==11)      /* moon_eph needed */
        ephMoon(mjd2000, r, v);
    else                /* NeoEphemeris needed */
    {
        if (ephNEO(mjd2000, ibody, kep, &mass, &meanAnom, name, &date)!=0)
        {
            error = 2;
            goto exit;
        }
    }
    
    
    /* If ephMoon has been called, the keplerian elements mush be
     * transformed into cartesian coordinates                            */
    
    if(ibody!=11 && error==0)           /* Planet or asteroid, Sun-centered */
    {
        /* Transform from Keplerian to cartesian */
        kep2car(kep, MU_SUN, 0., car);
        dcopy(car,r,3);
        dcopy(car+3,v,3);

    }

exit:
    return error;
}


int ephSS_kep(const int ibody, const double mjd2000, double kep[])
{
    /* Declaration */
    double r[3], v[3];  /* Cartesian position and velocity vectors */
    double car[6];      /* vector of both cartesian position and velocity */
    double eccAnom, meanAnom, dt, p, mass, date;    /* extra parameters needed to call
													CART2KEP or ephNEO */
	char name[81];		/* Used for ephNEO */
    int error = 0;          /* returned integer: 0 if no error
                               1 if the error occurs in uplanet
                               2 if the error occurs in ephNEO */

    
    if(ibody<11)            /* uplanet needed */
    {
        if (uplanet(mjd2000, ibody, kep)!=0)
        {
            error = 1;
            goto exit;
        }
    }
    else if(ibody==11)      /* moon_eph needed */
    {
        /* Returns the cartesian position and velocity */
        ephMoon(mjd2000, r, v);
        dcopy(r,car,3);
        dcopy(v,car+3,3);
        
	    /* Transform from cartesian to Keplerian */
        if(error == 0)    
            car2kep(car, MU_E, kep, &eccAnom, &meanAnom, &dt, &p);
    }
    else /* NeoEphemeris needed */
    {
        if (ephNEO(mjd2000, ibody, kep, &mass, &meanAnom, name, &date)!=0)
        {
            error = 2;
            goto exit;
        }
    }

exit:
    return error;
}
