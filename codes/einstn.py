#!/usr/bin/env python3
'''
###############################################################################
###############################################################################
##                                                                           ##
##     _    ___  ___   ___ ___ ___                                           ##
##    | |  | __ /   \ / __| _ | __|                                          ##
##    | |__| __  ( ) | (_ |  _|__ \                                          ##
##    |____|___ \___/ \___|_| \___/                                          ##
##                                    v 1.1 (Stable)                         ##
##                                                                           ##
## FILE DESCRIPTION:                                                         ##
##                                                                           ##
## Provides a suite of functions that correct for relativistic effects:      ##
##                                                                           ##
## 1. SHAPIRO DELAY:                                                         ##
## --> Due to the space time curvature produced by the gravitational field,  ##
## --> the Euclidean range travelled by the signal, which is computed by     ##
## --> 'posvel.py' must be corrected by the extra distance travelled.        ##
## --> Typically, Shapiro effects corrupt with about ~2cm ranging error.     ##
##                                                                           ##
## 2. CLOCK ADVANCE:                                                         ##
## --> The rate of advance of two identical clocks, one in the LEO satellite ##
## --> and the other on the GPS satellite, will differ due to differences in ##
## --> the gravitational potential (general relativity) and to the relative  ##
## --> speed between them (special relativity).                              ##
##                                                                           ##
## INPUTS:                                                                   ##
##                                                                           ##
## Position and velocities of LEO and GPS satellites as numpy arrays.        ##
##                                                                           ##
## OUTPUT:                                                                   ##
##                                                                           ##
## Relativistic orrections to the pseudorange in meters.                     ##
##                                                                           ##
## REMARKS:                                                                  ##
##                                                                           ##
## THANK GOD FOR DR ALBERT EINSTEIN!                                         ##
##                                                                           ##
## AUTHOR MODIFIED: 26-06-2019, by Samuel Y.W. Low                           ##
##                                                                           ##
###############################################################################
###############################################################################
'''

# Formulas below adapted from ESA's GNSS Guidebook

import math
import numpy as np

# IMPORT LOCAL LIBRARIES
from codes import consts

gm  = consts.GM  # GRAVITY CONSTANT * EARTH MASS
c   = consts.C   # VELOCITY OF LIGHT (M/S)
err = consts.ERR # # EARTH INERTIAL ROTATION RATE (RAD/S)

''' CLOCK ADVANCE EFFECT '''

def clockadv(gpspos,gpsvel):
    
    # This effect is more dominant than the Shapiro Delay effect.
    
    return -2*np.dot(gpspos,gpsvel)/c # Units in meters, thus c is not squared

''' SHAPIRO DELAY EFFECT '''

def shapiro(leopos,gpspos):
    
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
