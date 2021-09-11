#!/usr/bin/env python3

###############################################################################
###############################################################################
##                                                                           ##
##     _    ___  ___   ___ ___ ___                                           ##
##    | |  | __ /   \ / __| _ | __|                                          ##
##    | |__| __| ( ) | (_ |  _|__ \                                          ##
##    |____|___ \___/ \___|_| \___/                                          ##
##                                    v 1.3 (Stable)                         ##
##                                                                           ##
##    Provides functions that correct for relativistic effects:              ##
##                                                                           ##
##    1. SHAPIRO DELAY:                                                      ##
##    Due to the space time curvature produced by the gravitational field,   ##
##    the Euclidean range travelled by the signal, which is computed by      ##
##    'posvel.py' must be corrected by the extra distance travelled.         ##
##    Typically, Shapiro effects corrupt with about ~2cm ranging error.      ##
##                                                                           ##
##    2. CLOCK ADVANCE:                                                      ##
##    The rate of advance of two identical clocks, one in the LEO satellite  ##
##    and the other on the GPS satellite, will differ due to differences in  ##
##    the gravitational potential (general relativity) and to the relative   ##
##    speed between them (special relativity).                               ##
##                                                                           ##
##    Formulae below are adapted from the ESA GNSS Guidebook.                ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 08-Jun-2021.                                             ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##    Documentation: https://leogps.readthedocs.io/en/latest/                ##
##                                                                           ##
###############################################################################
###############################################################################

# Import global libraries
import numpy as np

# Import locla libraries
from source import consts

gm  = consts.GM  # GRAVITY CONSTANT * EARTH MASS
c   = consts.C   # VELOCITY OF LIGHT (M/S)
err = consts.ERR # # EARTH INERTIAL ROTATION RATE (RAD/S)

def clockadv(gpspos,gpsvel):
    '''Returns the clock advance effect on signal range as a path length
    correction in meters.

    Parameters
    ----------
    gpspos : numpy.ndarray
        1x3 position of GPS satellite (ITRF)
    gpsvel : numpy.ndarray
        1x3 velocity of GPS satellite (ITRF)

    Returns
    -------
    correction : float
        Units in meters

    '''
    
    # This effect is more dominant than the Shapiro Delay effect.
    
    return -2*np.dot(gpspos,gpsvel)/c # Units in meters, thus c is not squared

# SHAPIRO DELAY EFFECT.
def shapiro(leopos,gpspos):
    '''Returns the Shapiro delay effect on signal range as a path length
    correction in meters.

    Parameters
    ----------
    leopos : numpy.ndarray
        1x3 position of LEO satellite (ITRF)
    gpspos : numpy.ndarray
        1x3 position of GPS satellite (ITRF)

    Returns
    -------
    correction : float
        Units in meters

    '''
    # leo = geocentric distance between LEO satellite and Earth
    # gps = geocentric distance between GPS satellite and Earth
    # rel = relative distance between the receiving LEO and GPS satellite
    
    leo = np.linalg.norm(leopos)
    gps = np.linalg.norm(gpspos)
    rel = np.linalg.norm(gpspos-leopos)
    coeff = 2*gm/(c**2)
    ratio = ( gps + leo + rel ) / ( gps + leo - rel )
    shapiro_path_correction  = coeff * np.log( ratio )
    
    return shapiro_path_correction # Units in meters
