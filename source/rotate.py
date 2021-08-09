# -*- coding: utf-8 -*-

###############################################################################
###############################################################################
##                                                                           ##
##     _    ___  ___   ___ ___ ___                                           ##
##    | |  | __ /   \ / __| _ | __|                                          ##
##    | |__| __| ( ) | (_ |  _|__ \                                          ##
##    |____|___ \___/ \___|_| \___/                                          ##
##                                    v 1.2 (Stable)                         ##
##                                                                           ##
##    Module for rotation matrices supporting ICRF-to-ITRF conversion.       ##
##    Uses the IAU1976 Theory of Precession and IAU1980 Theory of Nutation.  ##
##                                                                           ##
##    References:                                                            ##
##    https://gssc.esa.int/navipedia/index.php/ICRF_to_CEP                   ##
##    https://gssc.esa.int/navipedia/index.php/CEP_to_ITRF                   ##
##    https://gssc.esa.int/navipedia/index.php/Julian_Date                   ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 10-Jun-2021                                              ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##    Documentation: https://leogps.readthedocs.io/en/latest/                ##
##                                                                           ##
###############################################################################
###############################################################################

import os
import math
import datetime
import numpy as np
from os.path import dirname, abspath, join

##############################################################################
##############################################################################
###                                                                        ###
###             FUNCTIONS BELOW FOR DIRECTION COSINE MATRICES              ###
###                                                                        ###
##############################################################################
##############################################################################

def _dcmX(t):
    '''Direction cosine matrix of local X-axis, with input "t" in radians'''
    
    dcm = np.array([[ 1.0,    0.0,       0.0       ],
                    [ 0.0,    math.cos(t), math.sin(t) ],
                    [ 0.0, -1*math.sin(t), math.cos(t) ]])

    return dcm

def _dcmY(t):
    '''Direction cosine matrix of local Y-axis, with input "t" in radians'''
    
    dcm = np.array([[ math.cos(t), 0.0, -1*math.sin(t) ],
                    [ 0.0,       1.0,    0.0       ],
                    [ math.sin(t), 0.0,    math.cos(t) ]])

    return dcm

def _dcmZ(t):
    '''Direction cosine matrix of local Z-axis, with input "t" in radians'''
    
    dcm = np.array([[    math.cos(t), math.sin(t), 0.0 ],
                    [ -1*math.sin(t), math.cos(t), 0.0 ],
                    [    0.0,       0.0,       1.0 ]])
    
    return dcm

##############################################################################
##############################################################################
###                                                                        ###
###          FUNCTIONS THAT RETURNS THE DAY-OF-WEEK AND GPS WEEK           ###
###                                                                        ###
##############################################################################
##############################################################################

# We define a function that returns the day-of-week and the GPS week.
def _gpsweekday(t):
    
    # Logic below calculates the desired GPS day and week number.
    wkday = (t.weekday() + 1) % 7 # Weekday from Python to GPST
    GPST_epoch = datetime.date(1980,1,6) # Date of GPST epoch
    user_epoch = t.date() # Get the date of the input time
    GPST_epoch_Monday = GPST_epoch - datetime.timedelta(GPST_epoch.weekday())
    user_epoch_Monday = user_epoch - datetime.timedelta(user_epoch.weekday())
    wwww = int(((user_epoch_Monday-GPST_epoch_Monday).days/7)-1) # GPS week
    
    # Algorithmic correction to the above logic.
    if wkday == 0:
        wwww += 1
        
    return str(wkday), str(wwww)

##############################################################################
##############################################################################
###                                                                        ###
###         FUNCTION BELOW TO COMPUTE PRECESSION ROTATION MATRIX           ###
###                                                                        ###
##############################################################################
##############################################################################

def precession(t):
    '''Computes the precession matrix using the IAU 1976 Precession Model, 
    where the matrix is applied in the direction from ICRF to the CEP.
    
    Parameters
    ----------
    t : datetime.datetime
        Current time of observation in GPST.
    
    Returns
    -------
    precession : numpy.ndarray
        A 3x3 DCM that rotates the ICRF by the precession offset.
    
    '''
    
    # Get a degree-to-radian multiplier.
    d2r = np.pi / 180.0
    
    # Convert the current epoch in GPST to Terrestrial Time (TDT)
    tt = t + datetime.timedelta(seconds=51)
    
    # Get the J2000 reference epoch
    j2000 = datetime.datetime(2000,1,1,12,0,0)
    
    # From the IAU 1976 precession model, we have the axial rotations:
    # p = (2306".2181 * T) + (1".09468 * T^2) + (0".018203 * T^3)
    # q = (2004".3109 * T) − (0".42665 * T^2) − (0".041833 * T^3)
    # r = (2306".2181 * T) + (0".30188 * T^2) + (0".017998 * T^3)
    # Where T is the time expressed in Julian centuries (36,525 Earth days),
    # between the reference epoch J2000 and the current time of observation.
    
    pqr = np.array([[2306.2181, 1.09468, 0.018203],
                    [2004.3109, 0.42665, 0.041833],
                    [2306.2181, 0.30188, 0.017998]])
    
    # Compute the time elapsed T in Julian centuries
    T = ( ( (tt - j2000).days ) + ( (tt - j2000).seconds / 86400 ) ) / 36525
    
    p = np.dot( pqr[0], np.array([ T, T**2, T**3 ]) ) * d2r / 3600
    q = np.dot( pqr[1], np.array([ T, T**2, T**3 ]) ) * d2r / 3600
    r = np.dot( pqr[2], np.array([ T, T**2, T**3 ]) ) * d2r / 3600
    
    # Compute the final precession matrix
    precession = _dcmZ( -1*p ) @ _dcmY( q ) @ _dcmZ( -1*r )
    
    return precession

##############################################################################
##############################################################################
###                                                                        ###
###          FUNCTION BELOW TO COMPUTE NUTATION ROTATION MATRIX            ###
###                                                                        ###
##############################################################################
##############################################################################

def nutation(t):
    '''Computes the nutation matrix using the IAU 1980 Nutation Model, 
    where the matrix is applied in the direction from ICRF to the CEP.
    
    Parameters
    ----------
    t : datetime.datetime
        Current time of observation in GPST.
    
    Returns
    -------
    nutation : numpy.ndarray
        A 3x3 DCM that rotates the ICRF by the nutation offset.
    
    '''
    
    # Get the coefficients to the IAU 1980 Nutation Model
    # Headers k1, k2, k3, k4, k5, A0 ("), A1 ("), B0 ("), B1 (")
    
    N = [[  0.0,  0.0,  0.0,  0.0,  1.0, -17.1996, -174.2,  9.2025,  8.9  ],
         [  0.0,  0.0,  2.0, -2.0,  2.0,  -1.3187,   -1.6,  0.5736, -3.1  ],
         [  0.0,  0.0,  2.0,  0.0,  2.0,  -0.2274,   -0.2,  0.0977, -0.5  ],
         [  0.0,  0.0,  0.0,  0.0,  2.0,   0.2062,    0.2, -0.0895,  0.5  ],
         [  0.0, -1.0,  0.0,  0.0,  0.0,  -0.1426,    3.4,  0.0054, -0.1  ],
         [  1.0,  0.0,  0.0,  0.0,  0.0,   0.0712,    0.1, -0.0007,  0.0  ],
         [  0.0,  1.0,  2.0, -2.0,  2.0,  -0.0517,    1.2,  0.0224, -0.6  ],
         [  0.0,  0.0,  2.0,  0.0,  1.0,  -0.0386,   -0.4,    0.02,  0.0  ],
         [  1.0,  0.0,  2.0,  0.0,  2.0,  -0.0301,    0.0,  0.0129, -0.1  ],
         [  0.0, -1.0,  2.0, -2.0,  2.0,   0.0217,   -0.5, -0.0095,  0.3  ],
         [ -1.0,  0.0,  0.0,  2.0,  0.0,   0.0158,    0.0, -0.0001,  0.0  ],
         [  0.0,  0.0,  2.0, -2.0,  1.0,   0.0129,    0.1,  -0.007,  0.0  ],
         [ -1.0,  0.0,  2.0,  0.0,  2.0,   0.0123,    0.0, -0.0053,  0.0  ],
         [  1.0,  0.0,  0.0,  0.0,  1.0,   0.0063,    0.1, -0.0033,  0.0  ],
         [  0.0,  0.0,  0.0,  2.0,  0.0,   0.0063,    0.0, -0.0002,  0.0  ],
         [ -1.0,  0.0,  2.0,  2.0,  2.0,  -0.0059,    0.0,  0.0026,  0.0  ],
         [ -1.0,  0.0,  0.0,  0.0,  1.0,  -0.0058,   -0.1,  0.0032,  0.0  ],
         [  1.0,  0.0,  2.0,  0.0,  1.0,  -0.0051,    0.0,  0.0027,  0.0  ],
         [ -2.0,  0.0,  0.0,  2.0,  0.0,  -0.0048,    0.0,  0.0001,  0.0  ],
         [ -2.0,  0.0,  2.0,  0.0,  1.0,   0.0046,    0.0, -0.0024,  0.0  ],
         [  0.0,  0.0,  2.0,  2.0,  2.0,  -0.0038,    0.0,  0.0016,  0.0  ],
         [  2.0,  0.0,  2.0,  0.0,  2.0,  -0.0031,    0.0,  0.0013,  0.0  ],
         [  2.0,  0.0,  0.0,  0.0,  0.0,   0.0029,    0.0, -0.0001,  0.0  ],
         [  1.0,  0.0,  2.0, -2.0,  2.0,   0.0029,    0.0, -0.0012,  0.0  ],
         [  0.0,  0.0,  2.0,  0.0,  0.0,   0.0026,    0.0, -0.0001,  0.0  ],
         [  0.0,  0.0,  2.0, -2.0,  0.0,  -0.0022,    0.0,     0.0,  0.0  ],
         [ -1.0,  0.0,  2.0,  0.0,  1.0,   0.0021,    0.0,  -0.001,  0.0  ],
         [  0.0,  2.0,  0.0,  0.0,  0.0,   0.0017,   -0.1,     0.0,  0.0  ],
         [  0.0,  2.0,  2.0, -2.0,  2.0,  -0.0016,    0.1,  0.0007,  0.0  ],
         [ -1.0,  0.0,  0.0,  2.0,  1.0,   0.0016,    0.0, -0.0008,  0.0  ],
         [  0.0,  1.0,  0.0,  0.0,  1.0,  -0.0015,    0.0,  0.0009,  0.0  ],
         [  1.0,  0.0,  0.0, -2.0,  1.0,  -0.0013,    0.0,  0.0007,  0.0  ],
         [  0.0, -1.0,  0.0,  0.0,  1.0,  -0.0012,    0.0,  0.0006,  0.0  ],
         [  2.0,  0.0, -2.0,  0.0,  0.0,   0.0011,    0.0,     0.0,  0.0  ],
         [ -1.0,  0.0,  2.0,  2.0,  1.0,   -0.001,    0.0,  0.0005,  0.0  ],
         [  1.0,  0.0,  2.0,  2.0,  2.0,  -0.0008,    0.0,  0.0003,  0.0  ],
         [  0.0, -1.0,  2.0,  0.0,  2.0,  -0.0007,    0.0,  0.0003,  0.0  ],
         [  0.0,  0.0,  2.0,  2.0,  1.0,  -0.0007,    0.0,  0.0003,  0.0  ],
         [  1.0,  1.0,  0.0, -2.0,  0.0,  -0.0007,    0.0,     0.0,  0.0  ],
         [  0.0,  1.0,  2.0,  0.0,  2.0,   0.0007,    0.0, -0.0003,  0.0  ],
         [ -2.0,  0.0,  0.0,  2.0,  1.0,  -0.0006,    0.0,  0.0003,  0.0  ],
         [  0.0,  0.0,  0.0,  2.0,  1.0,  -0.0006,    0.0,  0.0003,  0.0  ],
         [  2.0,  0.0,  2.0, -2.0,  2.0,   0.0006,    0.0, -0.0003,  0.0  ],
         [  1.0,  0.0,  0.0,  2.0,  0.0,   0.0006,    0.0,     0.0,  0.0  ],
         [  1.0,  0.0,  2.0, -2.0,  1.0,   0.0006,    0.0, -0.0003,  0.0  ],
         [  0.0,  0.0,  0.0, -2.0,  1.0,  -0.0005,    0.0,  0.0003,  0.0  ],
         [  0.0, -1.0,  2.0, -2.0,  1.0,  -0.0005,    0.0,  0.0003,  0.0  ],
         [  2.0,  0.0,  2.0,  0.0,  1.0,  -0.0005,    0.0,  0.0003,  0.0  ],
         [  1.0, -1.0,  0.0,  0.0,  0.0,   0.0005,    0.0,     0.0,  0.0  ],
         [  1.0,  0.0,  0.0, -1.0,  0.0,  -0.0004,    0.0,     0.0,  0.0  ],
         [  0.0,  0.0,  0.0,  1.0,  0.0,  -0.0004,    0.0,     0.0,  0.0  ],
         [  0.0,  1.0,  0.0, -2.0,  0.0,  -0.0004,    0.0,     0.0,  0.0  ],
         [  1.0,  0.0, -2.0,  0.0,  0.0,   0.0004,    0.0,     0.0,  0.0  ],
         [  2.0,  0.0,  0.0, -2.0,  1.0,   0.0004,    0.0, -0.0002,  0.0  ],
         [  0.0,  1.0,  2.0, -2.0,  1.0,   0.0004,    0.0, -0.0002,  0.0  ],
         [  1.0,  1.0,  0.0,  0.0,  0.0,  -0.0003,    0.0,     0.0,  0.0  ],
         [  1.0, -1.0,  0.0, -1.0,  0.0,  -0.0003,    0.0,     0.0,  0.0  ],
         [ -1.0, -1.0,  2.0,  2.0,  2.0,  -0.0003,    0.0,  0.0001,  0.0  ],
         [  0.0, -1.0,  2.0,  2.0,  2.0,  -0.0003,    0.0,  0.0001,  0.0  ],
         [  1.0, -1.0,  2.0,  0.0,  2.0,  -0.0003,    0.0,  0.0001,  0.0  ],
         [  3.0,  0.0,  2.0,  0.0,  2.0,  -0.0003,    0.0,  0.0001,  0.0  ],
         [ -2.0,  0.0,  2.0,  0.0,  2.0,  -0.0003,    0.0,  0.0001,  0.0  ],
         [  1.0,  0.0,  2.0,  0.0,  0.0,   0.0003,    0.0,     0.0,  0.0  ],
         [ -1.0,  0.0,  2.0,  4.0,  2.0,  -0.0002,    0.0,  0.0001,  0.0  ],
         [  1.0,  0.0,  0.0,  0.0,  2.0,  -0.0002,    0.0,  0.0001,  0.0  ],
         [ -1.0,  0.0,  2.0, -2.0,  1.0,  -0.0002,    0.0,  0.0001,  0.0  ],
         [  0.0, -2.0,  2.0, -2.0,  1.0,  -0.0002,    0.0,  0.0001,  0.0  ],
         [ -2.0,  0.0,  0.0,  0.0,  1.0,  -0.0002,    0.0,  0.0001,  0.0  ],
         [  2.0,  0.0,  0.0,  0.0,  1.0,   0.0002,    0.0, -0.0001,  0.0  ],
         [  3.0,  0.0,  0.0,  0.0,  0.0,   0.0002,    0.0,     0.0,  0.0  ],
         [  1.0,  1.0,  2.0,  0.0,  2.0,   0.0002,    0.0, -0.0001,  0.0  ],
         [  0.0,  0.0,  2.0,  1.0,  2.0,   0.0002,    0.0, -0.0001,  0.0  ],
         [  1.0,  0.0,  0.0,  2.0,  1.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  1.0,  0.0,  2.0,  2.0,  1.0,  -0.0001,    0.0,  0.0001,  0.0  ],
         [  1.0,  1.0,  0.0, -2.0,  1.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  0.0,  1.0,  0.0,  2.0,  0.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  0.0,  1.0,  2.0, -2.0,  0.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  0.0,  1.0, -2.0,  2.0,  0.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  1.0,  0.0, -2.0,  2.0,  0.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  1.0,  0.0, -2.0, -2.0,  0.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  1.0,  0.0,  2.0, -2.0,  0.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  1.0,  0.0,  0.0, -4.0,  0.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  2.0,  0.0,  0.0, -4.0,  0.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  0.0,  0.0,  2.0,  4.0,  2.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  0.0,  0.0,  2.0, -1.0,  2.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [ -2.0,  0.0,  2.0,  4.0,  2.0,  -0.0001,    0.0,  0.0001,  0.0  ],
         [  2.0,  0.0,  2.0,  2.0,  2.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  0.0, -1.0,  2.0,  0.0,  1.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  0.0,  0.0, -2.0,  0.0,  1.0,  -0.0001,    0.0,     0.0,  0.0  ],
         [  0.0,  0.0,  4.0, -2.0,  2.0,   0.0001,    0.0,     0.0,  0.0  ],
         [  0.0,  1.0,  0.0,  0.0,  2.0,   0.0001,    0.0,     0.0,  0.0  ],
         [  1.0,  1.0,  2.0, -2.0,  2.0,   0.0001,    0.0, -0.0001,  0.0  ],
         [  3.0,  0.0,  2.0, -2.0,  2.0,   0.0001,    0.0,     0.0,  0.0  ],
         [ -2.0,  0.0,  2.0,  2.0,  2.0,   0.0001,    0.0, -0.0001,  0.0  ],
         [ -1.0,  0.0,  0.0,  0.0,  2.0,   0.0001,    0.0, -0.0001,  0.0  ],
         [  0.0,  0.0, -2.0,  2.0,  1.0,   0.0001,    0.0,     0.0,  0.0  ],
         [  0.0,  1.0,  2.0,  0.0,  1.0,   0.0001,    0.0,     0.0,  0.0  ],
         [ -1.0,  0.0,  4.0,  0.0,  2.0,   0.0001,    0.0,     0.0,  0.0  ],
         [  2.0,  1.0,  0.0, -2.0,  0.0,   0.0001,    0.0,     0.0,  0.0  ],
         [  2.0,  0.0,  0.0,  2.0,  0.0,   0.0001,    0.0,     0.0,  0.0  ],
         [  2.0,  0.0,  2.0, -2.0,  1.0,   0.0001,    0.0, -0.0001,  0.0  ],
         [  2.0,  0.0, -2.0,  0.0,  1.0,   0.0001,    0.0,     0.0,  0.0  ],
         [  1.0, -1.0,  0.0, -2.0,  0.0,   0.0001,    0.0,     0.0,  0.0  ],
         [ -1.0,  0.0,  0.0,  1.0,  1.0,   0.0001,    0.0,     0.0,  0.0  ],
         [ -1.0, -1.0,  0.0,  2.0,  1.0,   0.0001,    0.0,     0.0,  0.0  ],
         [  0.0,  1.0,  0.0,  1.0,  0.0,   0.0001,    0.0,     0.0,  0.0  ]]
    
    # Get a degree-to-radian multiplier.
    d2r = np.pi / 180.0
    
    # Convert the current epoch in GPST to Terrestrial Time (TDT)
    tt = t + datetime.timedelta(seconds=51)
    
    # Get the J2000 reference epoch
    j2000 = datetime.datetime(2000,1,1,12,0,0)
    
    # Compute the time elapsed T in Julian centuries
    T = ( ( (tt - j2000).days ) + ( (tt - j2000).seconds / 86400 ) ) / 36525
    
    # Compute alpha-1 term corresponding to Moon's mean anomaly
    a1  = (485868.2490360000) * 1.0
    a1 += (715923.2177999020) * T
    a1 += (1325.00 * 1296000) * T
    a1 += (31.87920000000000) * T**2
    a1 += (0.051635000000000) * T**3
    a1 += (-0.00024470000000) * T**4
    a1  = (a1 / 3600) * d2r
    
    # Compute alpha-2 term corresponding to Sun's mean anomaly
    a2  = (1287104.793050000) * 1.0
    a2 += (1292581.048099995) * T
    a2 += (99.0000 * 1296000) * T
    a2 += (-0.55320000000000) * T**2
    a2 += (0.000136000000000) * T**3
    a2 += (-0.00001149000000) * T**4
    a2  = (a2 / 3600) * d2r
    
    # Compute alpha-3 term corresponding to Moon's mean argument of latitude
    a3  = (335779.5262320000) * 1.0
    a3 += (295262.8478000164) * T
    a3 += (1342.00 * 1296000) * T
    a3 += (-12.7512000000000) * T**2
    a3 += (-0.00103700000000) * T**3
    a3 += (0.000004170000000) * T**4
    a3  = (a3 / 3600) * d2r
    
    # Compute alpha-4 term corresponding to Moon's mean elongation from sun
    a4  = (1072260.703690000) * 1.0
    a4 += (1105601.209000110) * T
    a4 += (1236.00 * 1296000) * T
    a4 += (-6.37060000000000) * T**2
    a4 += (0.006593000000000) * T**3
    a4 += (-0.00003169000000) * T**4
    a4  = (a4 / 3600) * d2r
    
    # Compute alpha-5 term corresponding to longitude of ascending lunar node
    a5  = (450160.3980360000) * 1.0
    a5 += (-482890.543100000) * T
    a5 += (-5.0000 * 1296000) * T
    a5 += (7.472200000000000) * T**2
    a5 += (0.007702000000000) * T**3
    a5 += (-0.00005939000000) * T**4
    a5  = (a5 / 3600) * d2r
        
    # Initialise a nested list of coefficients for the IAU 1980 Nutation Model
    N_len = len(N)
    
    # Compute the first rotation angle, which is to be rotated about X-axis.
    # This angle 
    p = ( 23.439291 - 0.013004*T - 1.639E-7*(T**2) + 5.036E-7*(T**3) ) * d2r
    
    # Compute the second rotation angle, which is to be rotated about Z-axis.
    # This angle is the negative longitudinal nutation, which moves the X-axis
    # into its final position pointing towards the true equinox of date.
    q = 0.0
    for n in range(0, N_len):
        k1  = N[n][0]
        k2  = N[n][1]
        k3  = N[n][2]
        k4  = N[n][3]
        k5  = N[n][4]
        AA0 = N[n][5] * d2r / (3600)
        AA1 = N[n][6] * d2r / (3600*10000)
        delaunay_products = np.dot([k1,k2,k3,k4,k5], [a1,a2,a3,a4,a5])
        q += (AA0 + AA1*T) * math.sin(delaunay_products)
    
    # Compute the third rotation angle, which is to be rotated about X-axis.
    # This angle is the difference between the mean and true obliquity angle.
    # This angle, combined with the negative of the mean obliquity, gives the
    # rotation about Earth's new x-axis by the negative of the true obliquity.
    r = 0.0
    for n in range(0, N_len):
        k1  = N[n][0]
        k2  = N[n][1]
        k3  = N[n][2]
        k4  = N[n][3]
        k5  = N[n][4]
        BB0 = N[n][7] * d2r / (3600)
        BB1 = N[n][8] * d2r / (3600*10000)
        delaunay_products = np.dot([k1,k2,k3,k4,k5], [a1,a2,a3,a4,a5])
        r += (BB0 + BB1*T) * math.cos(delaunay_products)
    
    # Compute the final nutation matrix
    nutation = _dcmX(-1*(p+r)) @ _dcmZ(-1*q) @ _dcmX(p)
    
    return nutation

##############################################################################
##############################################################################
###                                                                        ###
###           FUNCTION BELOW TO COMPUTE DIURNAL ROTATION MATRIX            ###
###                                                                        ###
##############################################################################
##############################################################################

def diurnal(t, N):
    '''Computes the diurnal rotation (Earth's rotation about its own axis) to
    correct for the sidereal time motion since the ITRF is a rotating frame
    while the CEP is a non-rotating frame. 
    
    Parameters
    ----------
    t : datetime.datetime
        Current time of observation in GPST.
    N : numpy.ndarray
        A 3x3 DCM that rotates the ICRF by the nutation offset.

    Returns
    -------
    diurnal_rotation : numpy.ndarray
        A 3x3 DCM that rotates that converts the CEP frame into the rotating
        Earth frame due to sidereal motion, using the Greenwich Mean Sidereal
        Time (GMST) at 0-hours GPST (not UTC or UT1).

    '''
    
    # Note that the original formula in the ESA Navipedia uses UTC. However, 
    # LEOGPS will perform the diurnal rotation in GPST instead. Algorithm-
    # wise, there should be no change, since input `t` is already in GPST.
    
    # Get a degree-to-radian multiplier.
    d2r = np.pi / 180.0
    
    # Get the J2000 reference epoch
    j2000 = datetime.datetime(2000,1,1,12,0,0)
    
    # Compute the time elapsed T in Julian centuries
    T = (((t-j2000).days) + (((t-j2000).seconds) / 86400)) / 36525
    
    # Get the current GPST converted to radians
    gpst = (t.hour/24) + (t.minute/1440) + ((t.second)/(86400))
    gpst = gpst * 360.0 * d2r
    
    # Compute the Greenwich Mean Sidereal Time Offset (with respect to GPST)
    gmst_offset  = (6/24) + (41/1440) + (50.54841/86400)
    gmst_offset += (8640184.812866/86400) * T
    gmst_offset += (0.093104/86400) * T**2
    gmst_offset -= (6.2E-6) * T**3
    gmst_offset  = gmst_offset * 360.0 * d2r # Conversion to radians
    
    # Compute the Greenwich Mean Sidereal Time (with respect to GPST)
    gmst = gpst + gmst_offset
    
    # Compute the difference between the hour angle of the true equinox of
    # date and the mean equinox of date, due to the rotation from nutation
    ae = np.arctan2( N[0][1], N[0][0] )
    
    # Compute the final rotational argument
    theta = gmst - ae
    
    return _dcmZ( theta )

##############################################################################
##############################################################################
###                                                                        ###
###      FUNCTION BELOW TO COMPUTE DIURNAL ROTATION MATRIX DERIVATIVE      ###
###                                                                        ###
##############################################################################
##############################################################################

def diurnal_dot(t, diurnal):
    '''Computes the derivative of the diurnal rotation (Earth's rotation
    about its own axis) to correct for the sidereal time motion since the
    ITRF is a rotating frame while the CEP is a non-rotating frame.
    
    Parameters
    ----------
    t : datetime.datetime
        Current time of observation in GPST.
    diurnal : numpy.ndarray
        A 3x3 DCM that rotates that converts the CEP frame into the rotating
        Earth frame due to sidereal motion, using the Greenwich Mean Sidereal
        Time (GMST) at 0-hours GPST (not UTC or UT1).
    
    Returns
    -------
    diurnal_dot : numpy.ndarray
        A 3x3 DCM that rotates that converts the CEP frame into the rotating
        Earth frame due to sidereal motion, using the Greenwich Mean Sidereal
        Time (GMST) at 0-hours GPST (not UTC or UT1).

    '''
    
    # Compute Earth's rotation rate.
    earth_rotation = ( 2*math.pi ) / 86164.0
    
    # Get the elements of the diurnal rotation DCM.
    sin_term = diurnal[0][1]
    cos_term = diurnal[0][0]
    
    # Compute the first order derivative of the diurnal rotation DCM.
    diurnal_dot = np.array([[ -1*sin_term,    cos_term, 0.0 ],
                            [ -1*cos_term, -1*sin_term, 0.0 ],
                            [    0.0,         0.0,      0.0 ]])
    
    return earth_rotation * diurnal_dot
    
##############################################################################
##############################################################################
###                                                                        ###
###        FUNCTION BELOW TO COMPUTE POLAR MOTION ROTATION MATRIX          ###
###     NOTE THAT ONLY THE POLE WANDER COMPUTATION USES THE .ERP FILES     ###
###                                                                        ###
##############################################################################
##############################################################################

def polewander(t):
    '''Computes the pole wander rotation (motion of the Earth's spin axis).
    
    Parameters
    ----------
    t : datetime.datetime
        Current time of observation in GPST.
    
    Returns
    -------
    polewander : numpy.ndarray
        A 3x3 DCM that accounts for pole wander rotation of the ITRF.

    '''
    
    # First, we need to find the current Modified Julian Date (MJD)
    mjd_ref = datetime.datetime(1858,11,17,0,0,0)
    mjd_str = str((t - mjd_ref).days)
    
    # Before we start running the polar rotation matrix computation, we 
    # first need to access the ERP file to obtain the X-P & Y-P values.
    cwd = dirname(dirname(abspath(__file__))) # Current working directory (str)
    iwd = join(cwd, 'input') # Inputs files directory (str)
    os.chdir(iwd)
    
    # Get the GPS week number as it is part of the ERP file name.
    wd, wwww = _gpsweekday(t)
    fname = 'COD' + wwww + '7' + '.ERP'
    
    # Get a degree-to-radian multiplier.
    d2r = np.pi / 180.0
    
    # Search for the MJD in the ERP file, and get the X-P & Y-P values.
    erp_file = open(fname,'r')
    xp_yp_scaleflag = False
    xp_yp_scalevalue = 0.0
    for line in erp_file:
        linespl = line.split()
        if len(linespl) > 0:
            
            # First, we need to get the scale factor for Xp-Yp.
            if xp_yp_scaleflag == False:
                if 'MJD' in linespl[0]:
                    xp_yp_scaleflag = True
            else:
                if 'E' in linespl[0]:
                    xp_yp_scalevalue = '1' + linespl[0].replace('"','')
                    xp_yp_scalevalue = float(xp_yp_scalevalue)
            
            # Next, we need to get the Xp and Yp actual values.
            if mjd_str in linespl[0]:
                xp = xp_yp_scalevalue * float(linespl[1]) * d2r / 3600
                yp = xp_yp_scalevalue * float(linespl[2]) * d2r / 3600
    
    polewander = _dcmY(-1*xp) @ _dcmX(-1*yp)
    
    return polewander

##############################################################################
##############################################################################

# Example usage, let us assume the user has a CSV file in J2000 frame, with
# each row having 7 elements: time, positions x, y, z, and velocities x, y, z.
# Now, this would be an example routine on how to call the following functions
# to perform coordinate frame transformation from ICRF to CEP and CEP to ITRF.

# if __name__ == '__main__' :
    
#     import csv                              # Import CSV library
#     input_file  = 'J2000.csv'               # File name for input
#     output_file = 'ITRF.csv'                # File name for output
#     ti = datetime.datetime(2020,1,15,4,0,0) # Set an initial epoch
#     ts = datetime.timedelta(seconds=60)     # Set a time step value (s)
#     P = precession(ti)                      # Prepare precession matrix
#     N = nutation(ti)                        # Prepare nutation matrix
#     output = open(output_file, 'w')         # Open up output file
    
#     with open( input_file ) as csvf:        # Begin looping through CSV
        
#         csvr = csv.reader(csvf, delimiter=',')
        
#         for row in csvr:
            
#             if len(row) > 0:
                
#                 px = float(row[1])          # X-Axis Position in J2000 frame
#                 py = float(row[2])          # Y-Axis Position in J2000 frame
#                 pz = float(row[3])          # Z-Axis Position in J2000 frame
#                 vx = float(row[4])          # X-Axis Velocity in J2000 frame
#                 vy = float(row[5])          # Y-Axis Velocity in J2000 frame
#                 vz = float(row[6])          # Z-Axis Velocity in J2000 frame
                
#                 S  = diurnal( ti, N )       # Diurnal Rotation DCM 
#                 Sd = diurnal_dot( ti, S )   # Diurnal Rotation Dot DCM
#                 M  = polewander( ti )       # Pole Wander Rotation DCM
                
#                 pos = np.array([px,py,pz])  # Position Vector J2000
#                 vel = np.array([vx,vy,vz])  # Velocity Vector J2000
                
#                 pos_CEP = N @ P @ pos       # Position Vector CEP
#                 vel_CEP = N @ P @ vel       # Velocity Vector CEP
                
#                 pos_ITRF = M @ S @ pos_CEP                      # ITRF
#                 vel_ITRF = M @ ((Sd @ pos_CEP) + (S @ vel_CEP)) # ITRF
                
#                 line  = str(pos_ITRF[0]) + ', ' # Write X-Axis Position ITRF
#                 line += str(pos_ITRF[1]) + ', ' # Write Y-Axis Position ITRF
#                 line += str(pos_ITRF[2]) + ', ' # Write Z-Axis Position ITRF
#                 line += str(vel_ITRF[0]) + ', ' # Write X-Axis Velocity ITRF
#                 line += str(vel_ITRF[1]) + ', ' # Write Y-Axis Velocity ITRF
#                 line += str(vel_ITRF[2]) + '\n' # Write Z-Axis Velocity ITRF
#                 output.write(line)
                
#                 ti = ti + ts # Update the time step
                
#     output.close() # Close the output file
    