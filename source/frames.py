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
##    Module for conversion of coordinate frames (ICRF, ITRF, and LVLH)      ##
##    - ICRF - International Celestial Reference Frame (ECI)                 ##
##    - ITRF - International Terrestrial Reference Frame (ECEF)              ##
##    - LVLH - Local Vertical Local Horizontal Frame (Hill Frame, VCI)       ##
##                                                                           ##
##    Uses the IAU1976 Theory of Precession and IAU1980 Theory of Nutation.  ##
##                                                                           ##
##    References:                                                            ##
##    https://gssc.esa.int/navipedia/index.php/ICRF_to_CEP                   ##
##    https://gssc.esa.int/navipedia/index.php/CEP_to_ITRF                   ##
##    https://gssc.esa.int/navipedia/index.php/Julian_Date                   ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 09-Aug-2021                                              ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##    Documentation: https://leogps.readthedocs.io/en/latest/                ##
##                                                                           ##
###############################################################################
###############################################################################

# Import global libraries
import datetime
import numpy as np

# Import local libraries
#from source import rotate

from source import rotate

##############################################################################
##############################################################################

def icrf2cep(t, r, v = np.zeros(3)):
    '''Transformation of the international celestial reference frame (ICRF)
    to the conventional ephemeris pole frame (the True-Of-Epoch frame), by 
    correcting precession and nutation. This transformation is performed using 
    a composite of two orthogonal rotation matrices P and N. This function 
    will return two vectors (position and velocity).
    
    Parameters
    ----------
    t : datetime.datetime
        Current time of observation in GPST.
    r : numpy.ndarray
        Position vector (1x3) in ICRF frame.
    v : numpy.ndarray, optional
        Velocity vector (1x3) in ICRF frame.
    
    Returns
    -------
    r_cep : numpy.ndarray
        Position vector in CEP frame.
    v_cep : numpy.ndarray
        Velocity vector in CEP frame.
    
    '''
    
    P = rotate.precession(t) # Precession rotation DCM
    N = rotate.nutation(t) # Nutation rotation DCM
    
    if sum( abs( v ) ) == 0.0:
        r_cep = N @ P @ r
        return r_cep, np.zeros(3)
    else:
        r_cep = N @ P @ r
        v_cep = N @ P @ v
        return r_cep, v_cep

##############################################################################
##############################################################################

def cep2itrf(t, r, v = np.zeros(3)):
    '''Transformation of the conventional ephemeris pole frame (CEP) to the
    international terrestrial reference frame (ITRF) by accounting for the 
    diurnal rotation of the Earth, and accounting for the motion of the poles
    that matches the CEP to the ITRF. This function will return two vectors
    (position and velocity).

    Parameters
    ----------
    t : datetime.datetime
        Current time of observation in GPST.
    r : numpy.ndarray
        Position vector (1x3) in CEP frame.
    v : numpy.ndarray, optional
        Velocity vector (1x3) in CEP frame.

    Returns
    -------
    r_itrf : numpy.ndarray
        Position vector in ITRF frame.
    v_itrf : numpy.ndarray
        Velocity vector in ITRF frame.
    
    '''
    
    N = rotate.nutation(t)
    S = rotate.diurnal( t, N ) # Diurnal Rotation DCM 
    M = rotate.polewander( t ) # Pole Wander Rotation DCM
    
    if sum( abs( v ) ) == 0.0:
        r_itrf = M @ S @ r
        return r_itrf, np.zeros(3)
    else:
        Sd = rotate.diurnal_dot( t, S )
        r_itrf = M @ S @ r
        v_itrf = M @ ((Sd @ r) + (S @ v))
        return r_itrf, v_itrf
    
##############################################################################
##############################################################################

def itrf2cep(t, r, v = np.zeros(3)):
    '''Transformation of the international terrestrial reference frame (ITRF) 
    to the conventional ephemeris pole frame (CEP) by discounting for the 
    diurnal rotation of the Earth, and discounting the motion of the poles,
    from the ITRF to CEP. This function will return two vectors (position 
    and velocity).

    Parameters
    ----------
    t : datetime.datetime
        Current time of observation in GPST.
    r : numpy.ndarray
        Position vector (1x3) in ITRF frame.
    v : numpy.ndarray, optional
        Velocity vector (1x3) in ITRF frame.

    Returns
    -------
    r_itrf : numpy.ndarray
        Position vector in CEP frame.
    v_itrf : numpy.ndarray
        Velocity vector in CEP frame.
    
    '''
    
    N  = rotate.nutation(t)
    S  = rotate.diurnal( t, N )
    M  = rotate.polewander( t )
    Si = S.transpose()
    Mi = M.transpose()
    
    if sum( abs( v ) ) == 0.0:
        r_cep = Si @ Mi @ r
        return r_cep, np.zeros(3)
    else:
        Sd = rotate.diurnal_dot( t, S )
        r_cep = Si @ Mi @ r
        v_cep = Si @ (( Mi @ v ) - ( Sd @ r_cep ))
        return r_cep, v_cep

##############################################################################
##############################################################################

def cep2icrf(t, r, v = np.zeros(3)):
    '''Transformation of the conventional ephemeris pole frame (the True-Of-
    Epoch frame) to the international celestial reference frame (ICRF), by
    discounting precession and nutation. This transformation is performed 
    via a the inverse of the precession and nutation matrices P and N. This 
    function will return two vectors (position and velocity).
    
    Parameters
    ----------
    t : datetime.datetime
        Current time of observation in GPST.
    r : numpy.ndarray
        Position vector (1x3) in CEP frame.
    v : numpy.ndarray, optional
        Velocity vector (1x3) in CEP frame.
    
    Returns
    -------
    r_icrf : numpy.ndarray
        Position vector in ICRF frame.
    v_icrf : numpy.ndarray
        Velocity vector in ICRF frame.
    
    '''
    
    Pi = rotate.precession(t).transpose()
    Ni = rotate.nutation(t).transpose()
    
    if sum( abs( v ) ) == 0.0:
        r_icrf = Pi @ Ni @ r
        return r_icrf, np.zeros(3)
    else:
        r_icrf = Pi @ Ni @ r
        v_icrf = Pi @ Ni @ v
        return r_icrf, v_icrf
    
##############################################################################
##############################################################################

def itrf2icrf(t, r, v = np.zeros(3)):
    '''Transformation of the international terrestrial reference frame (ITRF)
    to the international celestial reference frame (ICRF), by calling the two
    functions in sequence: `itrf2cep()` and `cep2icrf()`.
    
    Parameters
    ----------
    t : datetime.datetime
        Current time of observation in GPST.
    r : numpy.ndarray
        Position vector (1x3) in ITRF frame.
    v : numpy.ndarray, optional
        Velocity vector (1x3) in ITRF frame.
    
    Returns
    -------
    r_icrf : numpy.ndarray
        Position vector in ICRF frame.
    v_icrf : numpy.ndarray
        Velocity vector in ICRF frame.
    
    '''
    
    r_cep,  v_cep  = itrf2cep( t, r, v )
    r_icrf, v_icrf = cep2icrf( t, r_cep, v_cep )
    return r_icrf, v_icrf

##############################################################################
##############################################################################

def icrf2itrf(t, r, v = np.zeros(3)):
    '''Transformation of the international celestial reference frame (ICRF) to
    the international terrestrial reference frame (ITRF), by calling the two
    functions in sequence: `icrf2cep` and `cep2itrf()`.
    
    Parameters
    ----------
    t : datetime.datetime
        Current time of observation in GPST.
    r : numpy.ndarray
        Position vector (1x3) in ICRF frame.
    v : numpy.ndarray, optional
        Velocity vector (1x3) in ICRF frame.
    
    Returns
    -------
    r_icrf : numpy.ndarray
        Position vector in ITRF frame.
    v_icrf : numpy.ndarray
        Velocity vector in ITRF frame.
    
    '''
    
    r_cep,  v_cep  = icrf2cep( t, r, v )
    r_itrf, v_itrf = cep2itrf( t, r_cep, v_cep )
    return r_itrf, v_itrf

##############################################################################
##############################################################################

def icrf2hill(baseline, rc, vc):
    '''Takes in a relative position vector, or baseline vector, as well as
    the chief position and velocity vectors. All inputs in ICRF. Transforms
    the relative position vector, or baseline vector, to the satellite local
    vertical local horizontal Euler-Hill Frame of the chief spacecraft.
    
    Parameters
    ----------
    baseline : numpy.ndarray
        Relative position vector (1x3) in ICRF frame.
    rc : numpy.ndarray
        Position vector (1x3) of Chief in ICRF frame.
    vc : numpy.ndarray
        Velocity vector (1x3) of Chief in ICRF frame.
    
    Returns
    -------
    hill_baseline : numpy.ndarray
        Relative position vector (1x3) of Deputy in Euler-Hill frame.
    
    '''
    
    # Construct the Euler-Hill frame basis vectors
    if sum( abs( vc ) ) != 0.0:
        
        h = np.cross(rc, vc)            # Angular momentum
        r_hat = rc / np.linalg.norm(rc) # Local X-axis
        h_hat = h / np.linalg.norm(h)   # Local Z-axis
        y_hat = np.cross(h_hat,r_hat)   # Local Y-axis
        
        # Compute the Hill DCM and transform the chief and deputy states.
        hill_dcm = np.array([ r_hat, h_hat, y_hat ])
        return hill_dcm @ baseline
        
    # Else, a Hill DCM cannot be created if velocity vector is invalid.
    else:
        return np.zeros(3)
    

##############################################################################
##############################################################################


# if __name__ == '__main__' :
    
#     import csv                              # Import CSV library
#     input_file  = 'OUTPUT.csv'               # File name for input
#     output_file = 'OUT2.csv'              # File name for output
#     ti = datetime.datetime(2020,1,15,4,0,0) # Set an initial epoch
#     ts = datetime.timedelta(seconds=60)     # Set a time step value (s)

#     output = open(output_file, 'w')         # Open up output file
    
#     with open(input_file) as csvf:          # Begin looping through CSV
        
#         csvr = csv.reader(csvf, delimiter=',')
        
#         for row in csvr:
            
#             if len(row) > 0:
                
#                 px = float(row[0])          # X-Axis Position in J2000 frame
#                 py = float(row[1])          # Y-Axis Position in J2000 frame
#                 pz = float(row[2])          # Z-Axis Position in J2000 frame
#                 vx = float(row[3])          # X-Axis Velocity in J2000 frame
#                 vy = float(row[4])          # Y-Axis Velocity in J2000 frame
#                 vz = float(row[5])          # Z-Axis Velocity in J2000 frame
                
#                 pos = np.array([px,py,pz])  # Position Vector J2000
#                 vel = np.array([vx,vy,vz])  # Velocity Vector J2000
                
#                 pos_CEP,  vel_CEP  = itrf2cep(ti, pos,     vel    )
#                 pos_ITRF, vel_ITRF = cep2icrf(ti, pos_CEP, vel_CEP)
                
#                 line  = str(pos_ITRF[0]) + ', ' # Write X-Axis Position ITRF
#                 line += str(pos_ITRF[1]) + ', ' # Write Y-Axis Position ITRF
#                 line += str(pos_ITRF[2]) + ', ' # Write Z-Axis Position ITRF
#                 line += str(vel_ITRF[0]) + ', ' # Write X-Axis Velocity ITRF
#                 line += str(vel_ITRF[1]) + ', ' # Write Y-Axis Velocity ITRF
#                 line += str(vel_ITRF[2]) + '\n' # Write Z-Axis Velocity ITRF
#                 output.write(line)
                
#                 ti = ti + ts # Update the time step
                
#     output.close() # Close the output file
    