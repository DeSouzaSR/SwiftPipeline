#!/usr/bin/env python
#-*- coding:utf-8 -*-

import math

# Global variable equivalent to swift.inc
NPLMAX = 51   # max number of planets, including the Sun
NTPMAX = 1001 # max number of test particles
    
#Size of the test particle integer status flag
#NSTATP Number of status parameters
NSTATP = 3
#NSTAT Number of status parameters
NSTAT = NSTATP + NPLMAX - 1 # include one for @ planet  

# Size of the test particle integer status flag
# integer NSTATR    
NSTATR = NSTAT # io_init_tp assumes NSTAT==NSTATR

# Convergence criteria for danby
# real DANBYAC , DANBYB
DANBYAC = 1.0e-14
DANBYB = 1.0e-13

# loop limits in the Laguerre attempts
# integer NLAG1, NLAG2
NLAG1 = 50
NLAG2 = 400

# A small number
# real*8 TINY
TINY = 4.0e-15

# trig stuff
# real*8 PI,TWOPI,PIBY2,DEGRAD
PI = 3.14159265358979e0
TWOPI = 2.0e0 * PI
PIBY2 = PI/2.0e0
DEGRAD = 180.0e0 / PI

def orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm):
    """
    *****************************************************************************
    *                          ORBEL_EL2XV.F
    *****************************************************************************
    *     PURPOSE: To compute cartesian positions and velocities given
    *               central mass, ialpha ( = +1 for hyp., 0 for para. and
    *               -1 for ellipse), and orbital elements.
    C       input:
    c            gm       ==> G times central mass (real scalar)
    c            ialpha   ==> conic section type ( see PURPOSE, integer scalar)
    C            a        ==> semi-major axis or pericentric distance if a parabola
    c                          (real scalar)
    c            e        ==> eccentricity (real scalar)
    C            inc      ==> inclination  (real scalar)
    C            capom    ==> longitude of ascending node (real scalar)
    C            omega    ==> argument of perihelion (real scalar)
    C            capm     ==> mean anomoly(real scalar)
    *       
    c       Output:
    c            x,y,z    ==>  position of object (real scalars)
    c            vx,vy,vz ==>  velocity of object (real scalars)
    c
    *     ALGORITHM:  See Fitzpatrick "Principles of Cel. Mech."
    *     REMARKS: All angles are in RADIANS
    *       
    *     AUTHOR:  M. Duncan.
    *     DATE WRITTEN:  May 11, 1992.
    *     REVISIONS: May 26 - now use better Kepler solver for ellipses
    *                 and hyperbolae called EHYBRID.F and FHYBRID.F
    ***********************************************************************
    """

    if e < 0.0:
        print(" ERROR in orbel_el2xv: e<0, setting e=0!")
        e = 0.0

    # Check for inconsistencies between ialpha and e
    em1 = e - 1.0e0
    if ((ialpha == 0) and (abs(em1) > TINY))  or \
        ((ialpha < 0) and (e > 1.0e0)) or \
        ((ialpha > 0) and (e < 1.0e0)):
        print("ERROR in orbel_el2xv: ialpha and e inconsistent")
        print("ialpha = ", ialpha)
        print("e = ", e)

    # Generate rotation matrices (on p. 42 of Fitzpatrick)
    sp, cp = orbel_scget(omega)
    so, co = orbel_scget(capom)
    si, ci = orbel_scget(inc)
    
    d11 = cp*co - sp*so*ci
    d12 = cp*so + sp*co*ci
    d13 = sp*si
    d21 = -sp*co - cp*so*ci
    d22 = -sp*so + cp*co*ci
    d23 = cp*si

    # Get the other quantities depending on orbit type ( i.e. IALPHA)
    if ialpha == -1:
        cape = orbel_ehybrid(e,capm)
        scap, ccap = orbel_scget(cape)
        sqe = math.sqrt(1.0e0 - e*e)
        sqgma = math.sqrt(gm * a)
        xfac1 = a*(ccap - e)
        xfac2 = a*sqe*scap
        ri = 1.0e0 / (a*(1.0e0 - e*ccap))
        vfac1 = -ri * sqgma * scap
        vfac2 = ri * sqgma * sqe * ccap

    elif ialpha == 1:
        capf = orbel_fhybrid(e,capm)
        shcap, chcap = orbel_schget(capf)
        sqe = math.sqrt(e*e - 1.0e0 )
        sqgma = math.sqrt(gm*a)
        xfac1 = a*(e - chcap)
        xfac2 = a*sqe*shcap
        ri = 1.0e0 / (a*(e*chcap - 1.0e0))
        vfac1 = -ri * sqgma * shcap
        vfac2 = ri * sqgma * sqe * chcap

    else:
        zpara = orbel_zget(capm)
        sqgma = math.sqrt(2.0e0*gm*a)
        xfac1 = a*(1.0e0 - zpara*zpara)
        xfac2 = 2.0e0*a*zpara
        ri = 1.0e0/(a*(1.0e0 + zpara*zpara))
        vfac1 = -ri * sqgma * zpara
        vfac2 = ri * sqgma

    x =  d11*xfac1 + d21*xfac2
    y =  d12*xfac1 + d22*xfac2
    z =  d13*xfac1 + d23*xfac2
    vx = d11*vfac1 + d21*vfac2
    vy = d12*vfac1 + d22*vfac2
    vz = d13*vfac1 + d23*vfac2

    return x, y, z, vx, vy, vz

    # End orbel_el2xv

def orbel_scget(angle):
    """
    ***********************************************************************
    c                         ORBEL_SCGET.F
    ***********************************************************************
    *     PURPOSE:  Given an angle, efficiently compute sin and cos.
    *
    *        Input:
    *             angle ==> angle in radians (real scalar)
    *        
    *        Output:
    *             sx    ==>  sin(angle)  (real scalar)
    *             cx    ==>  cos(angle)  (real scalar)
    *
    *     ALGORITHM: Obvious from the code 
    *     REMARKS: The HP 700 series won't return correct answers for sin
    *       and cos if the angle is bigger than 3e7. We first reduce it
    *       to the range [0,2pi) and use the sqrt rather than cos (it's faster)
    *       BE SURE THE ANGLE IS IN RADIANS - NOT DEGREES!
    *     AUTHOR:  M. Duncan.
    *     DATE WRITTEN:  May 6, 1992.
    *     REVISIONS: 
    ***********************************************************************
    """

    PI3BY2 = 1.5e0*PI
    
    nper = int(angle / TWOPI)
    x = angle - nper * TWOPI
    
    if x < 0.0e0:
        x = x + TWOPI
    sx = math.sin(x)
    cx = math.sqrt(1.0e0 - sx*sx)
       
    if x > PIBY2 and x < PI3BY2:
        cx = -cx

    return sx, cx
    # end orbel_scget

def orbel_ehybrid(e,m):
    """
    ***********************************************************************
    c                    ORBEL_EHYBRID.F
    ***********************************************************************
    *     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
    *
    *             Input:
    *                           e ==> eccentricity anomaly. (real scalar)
    *                           m ==> mean anomaly. (real scalar)
    *             Returns:
    *              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
    *
    *     ALGORITHM: For e < 0.18 uses fast routine ESOLMD 
    *                For larger e but less than 0.8, uses EGET
    *                For e > 0.8 uses EHIE
    *     REMARKS: Only EHIE brings M and E into range (0,TWOPI)
    *     AUTHOR: M. Duncan 
    *     DATE WRITTEN: May 25,1992.
    *     REVISIONS: 2/26/93 hfl
    ***********************************************************************
    """
    if e < 0.18e0:
        orbel_ehybrid = orbel_esolmd(e,m)
    else:
        if e <= 0.8e0:
            orbel_ehybrid = orbel_eget(e,m)
        else:
            orbel_ehybrid = orbel_ehie(e,m) 

    return orbel_ehybrid

    # End orbel_ehybrid

def orbel_fhybrid(e, n):
    """
    ***********************************************************************
    c                    ORBEL_FHYBRID.F
    ***********************************************************************
    *     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
    *
    *             Input:
    *                           e ==> eccentricity anomaly. (real scalar)
    *                           n ==> hyperbola mean anomaly. (real scalar)
    *             Returns:
    *               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
    *
    *     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON 
    *                For larger N, uses FGET
    *     REMARKS: 
    *     AUTHOR: M. Duncan 
    *     DATE WRITTEN: May 26,1992.
    *     REVISIONS: 
    *     REVISIONS: 2/26/93 hfl
    ***********************************************************************
    """
    
    abn = n
    
    if n < 0.0e0:
        abn = -abn

    if abn < 0.636e0*e -0.6e0:
        orbel_fhybrid = orbel_flon(e,n)
    else:
        orbel_fhybrid = orbel_fget(e,n)

    return orbel_fhybrid

    # end orbel_fhybrid

def orbel_zget(q):
    """
    ***********************************************************************
    c                    ORBEL_ZGET.F
    ***********************************************************************
    *     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola 
    *          given Q (Fitz. notation.)
    *
    *             Input:
    *                           q ==>  parabola mean anomaly. (real scalar)
    *             Returns:
    *                  orbel_zget ==>  eccentric anomaly. (real scalar)
    *
    *     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
    *     REMARKS: For a parabola we can solve analytically.
    *     AUTHOR: M. Duncan 
    *     DATE WRITTEN: May 11, 1992.
    *     REVISIONS: May 27 - corrected it for negative Q and use power
    *             series for small Q.
    ***********************************************************************
    """
    iflag = 0
    
    if q < 0.0e0:
        iflag = 1
        q = -q
    
    if q < 1.0e-3:
        orbel_zget = q*(1.0e0 - (q*q/3.0e0)*(1.0e0 -q*q))
    else:
        x = 0.5e0*(3.0e0*q + sqrt(9.0e0*(q**2) +4.0e0))
        tmp = x**(1.0e0 / 3.0e0)
          orbel_zget = tmp - 1.0e0 / tmp
    
    if iflag == 1:
        orbel_zget = -orbel_zget
        q = -q
       
    return orbel_zget
    # end orbel_zget

def orbel_esolmd(e,m):
    """
    ***********************************************************************
    c                    ORBEL_ESOLMD.F
    ***********************************************************************
    *     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
    *
    *             Input:
    *                           e ==> eccentricity anomaly. (real scalar)
    *                           m ==> mean anomaly. (real scalar)
    *             Returns:
    *                orbel_esolmd ==>  eccentric anomaly. (real scalar)
    *
    *     ALGORITHM: Some sort of quartic convergence from Wisdom. 
    *     REMARKS: ONLY GOOD FOR SMALL ECCENTRICITY SINCE IT ONLY
    *         ITERATES ONCE. (GOOD FOR PLANET CALCS.)
    *               ALSO DOES NOT PUT M OR E BETWEEN 0. AND 2*PI 
    *     INCLUDES: needs SCGET.F
    *     AUTHOR: M. Duncan 
    *     DATE WRITTEN: May 7, 1992.
    *     REVISIONS: 2/26/93 hfl
    ***********************************************************************
    """
    

