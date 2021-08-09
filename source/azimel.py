#!/usr/bin/env python3

###############################################################################
###############################################################################
##                                                                           ##
##     _    ___  ___   ___ ___ ___                                           ##
##    | |  | __ /   \ / __| _ | __|                                          ##
##    | |__| __| ( ) | (_ |  _|__ \                                          ##
##    |____|___ \___/ \___|_| \___/                                          ##
##                                    v 1.2 (Stable)                         ##
##                                                                           ##
##                                                                           ##
##    Conversion of ECEF XYZ coordinates to WGS84 LLA coordinates, and       ##
##    then the calculation of the azimuth and elevation with respect         ##
##    from the LEO to the GPS satellite.                                     ##
##                                                                           ##
##    Inputs: 3x1 ITRF position vector of a LEO and a GPS satellite.         ##
##    Outputs: Azimuth and elevation (rad) from LEO to GPS satellite.        ##
##                                                                           ##
##    The original efficient algorithm for conversion from ECEF to LLA       ##
##    was taken from the following sources, by original author Olson, D.K.   ##
##    (1996): "Converting Earth-Centered Earth-Fixed Coordinates to          ##
##    Geodetic Coordinates". IEEE Transactions on Aerospace and Electronic   ##
##    Systems, 32 (1996) 473-476.                                            ##
##                                                                           ##
##    Original Algorithm By: Olson, D.K. (1996)                              ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 08-Jun-2021.                                             ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##    Documentation: https://leogps.readthedocs.io/en/latest/                ##
##                                                                           ##
###############################################################################
###############################################################################

import math
import numpy as np

def azimel(leopos,gpspos):
    '''Returns azimuth and elevation angles (rad) from the LEO satellite to
    the GPS satellite when providing an ECEF position coordinate of them both.
    
    Parameters
    ----------
    leopos : list
        Position coordinate of LEO in ECEF [X,Y,Z]
    gpspos : list
        Position coordinate of GPS in ECEF [X,Y,Z]

    Returns
    -------
    az : float
        Azimuth angle (rad, -pi to +pi)
    el : float
        Azimuth angle (rad, -pi to +pi)

    '''
    
    # WGS-84 ellipsoid parameters
    a = 6378137
    f = 1 / 298.257223563
     
    # Derived parameters
    e2 = f * (2 - f)
    a1 = a * e2
    a2 = a1 * a1
    a3 = a1 * e2 / 2
    a4 = 2.5 * a2
    a5 = a1 + a3
    a6 = 1 - e2
    
    # Convert ECEF (meters) to LLA (radians and meters)
    
    def ecef2lla(leopos):
    
        w = math.sqrt(leopos[0] * leopos[0] + leopos[1] * leopos[1])
        z = leopos[2]
        zp = abs(z)
        w2 = w * w
        r2 = z * z + w2
        r  = math.sqrt(r2)
        s2 = z * z / r2
        c2 = w2 / r2
        u = a2 / r
        v = a3 - a4 / r
        
        if c2 > 0.3:
            s = (zp / r) * (1 + c2 * (a1 + u + s2 * v) / r)
            if s > 1:
                s = 1
            if s < -1:
                s = -1
            lat = math.asin(s)
            ss = s * s
            c = math.sqrt(1 - ss)
            
        else:
            c = (w / r) * (1 - s2 * (a5 - u - c2 * v) / r)
            if c > 1:
                c = 1
            if c < -1:
                c = -1
            lat = math.acos(c)
            ss = 1 - c * c
            s = math.sqrt(ss)
            
        g = 1 - e2 * ss
        rg = a / math.sqrt(g)
        rf = a6 * rg
        u = w - rg * c
        v = zp - rf * s
        f = c * u + s * v
        m = c * v - s * u
        p = m / (rf / g + f)
        lat = lat + p
        
        if z < 0:
            lat = -lat
            
        return np.array([lat, math.atan2(leopos[1], leopos[0]), f + m * p / 2])

    # Get latitude / longitude coordinates
    
    lla = ecef2lla(leopos)
    
    lat,lon = lla[0],lla[1]
    
    # Construct the ECEF-to-ENU DCM
    
    rotmat1 = np.array([-1*math.sin(lon),
                        math.cos(lon),
                        0.0])
    
    rotmat2 = np.array([-1 * math.sin(lat) * math.cos(lon),
                        -1 * math.sin(lat) * math.sin(lon),
                        math.cos(lat)])
    
    rotmat3 = np.array([math.cos(lat) * math.cos(lon),
                        math.cos(lat) * math.sin(lon),
                        math.sin(lat)])
    
    rotmat = np.array([rotmat1,rotmat2,rotmat3])
    
    # Rotate the relative position vector
    xRel = np.array(gpspos) - np.array(leopos)
    
    xENU = rotmat @ xRel
    
    xENUnorm = np.linalg.norm(xENU)
    
    # Get the angles
    
    az = np.arctan2( xENU[0] , xENU[1] ) # Azimuth angle
    
    el = np.arcsin( xENU[2] / xENUnorm ) # Elevation angle
    
    return az, el