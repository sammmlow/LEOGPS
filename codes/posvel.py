#!/usr/bin/env python3
'''
###############################################################################
###############################################################################
##                                                                           ##
##     _    ___  ___   ___ ___ ___                                           ##
##    | |  | __ /   \ / __| _ | __|                                          ##
##    | |__| __  ( ) | (_ |  _|__ \                                          ##
##    |____|___ \___/ \___|_| \___/                                          ##
##                                    v 1.0 (Stable)                         ##
##                                                                           ##
## FILE DESCRIPTION:                                                         ##
##                                                                           ##
## Code-based position estimation via weighted least squares estimation.     ##
## Doppler-based estimation of velocities only if Doppler is available.      ##
## Ionospheric delay for L1: Offset using the GRAPHIC linear combination.    ##
## Ionospheric delay for L2: Dual-frequency iono-free linear combination.    ##
## Clock sync done by extracting interpolating clock biases and drifts.      ##
## Relativistic effects estimated: Shapiro effect + clock delay/advance.     ##
## This script 'posvel.py' is called up for every epoch in 'leorun.py'       ##
##                                                                           ##
## INPUTS:                                                                   ##
##                                                                           ##
## epoch = datetime.datetime(2018,10,16,07,30,30), datetime object example   ##
## goodsats = [1,2,3... ...30,31,32], a list of good GPS satellites as ints  ##
## gps = {1:{'px':123, 'py':123, 'pz':123,...}, ... 32:(...)}                ##
## rxi = {1:('C1':xxx,'P2':xxx,'L1':xxx,'L2':xxx) ... 32:(...)}              ##
## inps = dictionary of inputs from user-specified parameters in config.txt  ##
## iters = number of iterations to run multivariate linear regression        ##
##                                                                           ##
## OUTPUT:                                                                   ##
##                                                                           ##
## Three arrays of position XYZ, velocity XYZ, and an array for G/P/TDOP     ##
##                                                                           ##
## REMARKS: Use RINEX v2.xx format, GPS only, GLONASS/GALILEO unsupported    ##
##                                                                           ##
## AUTHOR MODIFIED: 29-02-2020, by Samuel Y.W. Low                           ##
##                                                                           ##
###############################################################################
###############################################################################
'''

import datetime
import numpy as np

# IMPORT LOCAL LIBRARIES
from codes import consts
from codes import einstn

''' This is a position-velocity estimator for the SINGLE-EPOCH case. '''

def posvel(epoch, goodsats, gps, rxi, inps, nm, iters = 6):
    
    # gps = {1:{'px':123, 'py':123, 'pz':123,...}, ... 32:(...)}
    # rxi = {1:('C1':xxx,'P2':xxx,'L1':xxx,'L2':xxx) ... 32:(...)}
    
    ''' Before starting, we need at least 4 GPS satellites to solve XYZ! '''
    
    posf = np.array([0.0, 0.0, 0.0, 0.0]) # Default zero vectrix
    velf = np.array([0.0, 0.0, 0.0, 0.0]) # Default zero vectrix
    dopf = np.array([0.0, 0.0, 0.0]) # Default zero vectrix
    clkb = np.array([0.0]) # Default zero vectrix
    
    if len(rxi) < 4:
        print('Error! For the epoch of ' + str(nm) + ' ' + str(epoch))
        print('Only ' + str(len(rxi)) + ' GPS satellites observed!')
        print('Unable to obtain a solution fix, will output float zeros! \n')
        return posf, velf, dopf, clkb
    
    ''' First, initialise variables and constants '''
    
    C    = consts.C             # Speed of light in m/s
    err  = consts.ERR           # Earth rotation rate
    F1   = consts.FREQ1         # L1 frequency 1575420000 Hz
    F2   = consts.FREQ2         # L1 frequency 1227600000 Hz
    WL1  = C/F1                 # L1 wavelength
    WL2  = C/F2                 # L2 wavelength
    
    ''' Now, extract inputs from *config.txt* '''
    
    rel        = inps['relativity'] # Enable relativistic correction?
    erp        = inps['earthrotation'] # Offset Earth rotation?
    freqnum    = inps['freq'] # Check if single or dual frequency
    offsetX    = float(inps['antoffsetX']) # Antenna offset X
    offsetY    = float(inps['antoffsetY']) # Antenna offset Y
    offsetZ    = float(inps['antoffsetZ']) # Antenna offset Z
    
    ''' Initialise the core parameters used for estimation later on. '''
    
    Xp,Yp,Zp,Tp = 1.0, 1.0, 1.0, 1.0 # Initialise XYZT, for position.
    Xv,Yv,Zv,Tv = 0.0, 0.0, 0.0, 0.0 # Initialise XYZT, for velocity.
    
    ''' Now, we extract what information we need into a dictionary. '''
    
    gps_dict = {} # Main dictionary we reference from.
    
    # Now parse through all the SVs found in the RINEX file and GPS data.
    for p in rxi:
        
        if freqnum == 1: # Get RINEX observables if single frequency.
            
            # Checking for pseudorange information...
            if 'P1' in rxi[p]:
                obs_P1 = rxi[p]['P1'] # P1 ...
            elif 'C1' in rxi[p]:
                obs_P1 = rxi[p]['C1'] # ... or C1
            else:
                print('C1 or P1 pseudorange is missing this epoch.')
                print(str(epoch) +' '+ str(nm) + ' for SV' + str(p))
                print(' ')
                return posf, velf, dopf, clkb
            if 'D1' in rxi[p]:
                if rxi[p]['D1'] != 'NaN':
                    obs_D1 = (-1)*rxi[p]['D1'] # D1
                else:
                    obs_D1 = 0.0 # Invalid D1
            else:
                print('D1 doppler values are missing this epoch.')
                print(str(epoch) +' '+ str(nm) + ' for SV' + str(p))
                print(' ')
                return posf, velf, dopf, clkb
            if 'L1' in rxi[p]:
                obs_L1 = rxi[p]['L1'] * WL1
            else:
                print('L1 carrier phases are missing this epoch.')
                print(str(epoch) +' '+ str(nm) + ' for SV' + str(p))
                print(' ')
                return posf, velf, dopf, clkb
            
            # GRAPHIC single frequency combination (T. P. Yunck, 1993)
            pr_rnx = 0.5 * (obs_P1 + obs_L1)
            pv_rnx = obs_D1 * WL1
            
        # Are we using dual frequency data?
        elif freqnum == 2:
            
            # Checking for pseudorange information...
            if 'P1' in rxi[p] and 'P2' in rxi[p]:
                obs_P1 = rxi[p]['P1'] # P1
                obs_P2 = rxi[p]['P2'] # P2
            elif 'C1' in rxi[p] and 'P2' in rxi[p]:
                obs_P1 = rxi[p]['C1'] # C1
                obs_P2 = rxi[p]['P2'] # P2
            else:
                print('Dual frequency P1/P2 values are missing.')
                print(str(epoch) +' '+ str(nm) + ' for SV' + str(p))
                print(' ')
                return posf, velf, dopf, clkb
            if 'D1' in rxi[p] and 'D2' in rxi[p]:
                if rxi[p]['D1'] != 'NaN' and rxi[p]['D2'] != 'NaN':
                    obs_D1 = (-1)*rxi[p]['D1'] # D1
                    obs_D2 = (-1)*rxi[p]['D2'] # D2
                else:
                    obs_D1 = 0.0 # Invalid D1
                    obs_D2 = 0.0 # Invalid D2
            else:
                print('Dual frequency D1/D2 values are missing.')
                print(str(epoch) +' '+ str(nm) + ' for SV' + str(p))
                print(' ')
                return posf, velf, dopf, clkb
            
            # Iono-free dual frequency linear combination.
            pr_rnx = (2.546*obs_P1 - 1.546*obs_P2)
            pv_rnx = (2.546*obs_D1*WL1 - 1.546*obs_D2*WL2)
            
        else:
            print('The input text file has an invalid frequency!')
            return posf, velf, dopf, clkb
        
        # Now, we get all ephemeris data, time biases, and RINEX observables.
        gps_dict[p] = {}
        gps_dict[p]['pos'] = np.array([gps[p]['px'],gps[p]['py'],gps[p]['pz']])
        gps_dict[p]['pos'] = gps_dict[p]['pos'] * 1000 # Convert to meters.
        gps_dict[p]['vel'] = np.array([gps[p]['vx'],gps[p]['vy'],gps[p]['vz']])
        gps_dict[p]['vel'] = gps_dict[p]['vel'] * 1000 # Convert to m/s
        gps_dict[p]['clb'] = C * gps[p]['clkb'] # P-th GPS clock bias (m)
        gps_dict[p]['cld'] = C * gps[p]['clkd'] # P-th GPS clock drift (m)
        gps_dict[p]['rng'] = 0.0 # Initialise p-th GPS pseudorange
        gps_dict[p]['dop'] = 0.0 # Initialise p-th GPS doppler
        
        # Check for valid observables.
        if pr_rnx != 0.0:
            gps_dict[p]['rng'] = pr_rnx
        if pv_rnx != 0.0:
            gps_dict[p]['dop'] = pv_rnx + ( C * gps[p]['clkd'] )
    
    ''' Begin two runs of iterative least squares via gradient descent. '''
        
    # Initialise the LEO position vector
    leopos = np.array([Xp,Yp,Zp]) # Initial position vector
    leovel = np.array([Xv,Yv,Zv]) # Initial velocity vector
    
    for i in range(iters):
        
        A_mat = [] # Geometry or design matrix
        bp_mat = [] # Position-time vectrix
        p_mat = [] # Array to store the SVN IDs of GPS satellites
        
        # Generate the arrays holding the least-squares variables for velocity
        # estimation, only in the final iteration (velocity estimation does
        # not require iterative least squares, just the normal least squares)
        if i == (iters-1):
            bv_mat = [] # Pseudorange rate vectrix of LEO
            gv_mat = [] # Velocity vectrix of GPS satellite
            Av_mat = [] # Geometry matrix for velocity estimation
        
        # Begin parsing ephemeris, clock data, and observables.
        for p in gps_dict:
            
            # Store the SVN ID
            p_mat.append(p)
            
            # Initialise the GPS parameters
            gpspos = gps_dict[p]['pos'] # Position of p-th GPS satellite
            gpsvel = gps_dict[p]['vel'] # Velocity of p-th GPS satellite
            gpsrng = gps_dict[p]['rng'] # Pseudorange (code-based) (m)
            gpsdop = gps_dict[p]['dop'] # Doppler (pseudorange rate)
            gpsclb = gps_dict[p]['clb'] # Clock bias of p-th GPS satellite
            
            # Apply measurement model to obtain true range estimate
            # Note that an ionosphere free linear combination was already
            # applied (GRAPHIC was used even in the L1 only case)
            gpsrng_obsv  = gpsrng     # Observed pseudorange
            gpsrng_obsv += gpsclb     # GPS satellite clock bias
            gpsrng_obsv -= Tp         # LEO satellite clock bias
            
            # Compute the signal time of flight from GPS satellite to LEO
            sigTOF = gpsrng_obsv / C
            
            # The GPS positions are given in ECEF but at the transmission 
            # time of the GPS satellites, and at the reception time of the
            # LEO satellites. Two compensating calculations must be applied
            # here, the first being 
            
            # even if the GPS satellite was inertially not moving, it would
            # have experienced a change in position vector by virtue of Earth's
            # rotation when its "stationarity" was converted into the
            # rotating fixed frame. In the second thing, the GPS satellite is
            # also moving inertially and so this also has to be compensated for.
            
            if erp == 'True':
                
                # Rotation rate * TOF
                r = err * sigTOF 
                
                # Get the Earth rotation direction cosine matrix, with
                # small angle approximations for simplicity.
                e_rate = np.array([[  1.0,   r, 0.0 ],
                                   [ -1*r, 1.0, 0.0 ],
                                   [  0.0, 0.0, 1.0 ]])
                
                # Rotate the coordinate frame for GPS satellite position
                gpspos = np.matmul(e_rate, gpspos)
                
            # Compensate for GPS satellite motion during signal TOF
            gpspos = gpspos - (gpsvel * sigTOF)
            
            # Calculate the relativistic effects if 'rel' flag is true.
            if rel == 'True':
                
                # Compute the relative Clock Advance (~15m for statics)
                clockadv  = einstn.clockadv(gpspos,gpsvel)
                clockadv -= einstn.clockadv(leopos,leovel)
                
                # Compute the Shapiro delay (~0.015m usually)
                shapiro = einstn.shapiro(leopos,gpspos) 
                
            else:
                shapiro, clockadv = 0.0, 0.0
            
            # Update the range measurement model with relativistic effects
            gpsrng_obsv += shapiro    # Relativistic shapiro path delay
            gpsrng_obsv += clockadv   # Relativistic clock advance
            
            # Relative position vector between GPS & LEO satellite
            delpos = leopos - gpspos # Relative position vector
            delpos_norm = np.linalg.norm(delpos)
            
            # Compute the "true range" from current LEO position estimate
            gpsrng_comp = delpos_norm
            
            # Compute the residual (observed minus the computed range)
            gpsrng_resd = gpsrng_obsv - gpsrng_comp
            
            # Calculate partial derivatives of the i-th row in Jacobian
            partderiv = delpos / delpos_norm
            partderiv = np.append(partderiv,1)
            
            # Construct the geometry matrix A and pseudorange residual 'b'
            A_mat.append(partderiv) # Add this SV into Jacobian
            bp_mat.append(gpsrng_resd) # Add this SV into residual vector
            
            # Update the pseudorange rate vector and matrices.
            # We only need to do it for iter = 0, as velocity estimation
            # can be done with a directly least-squares, whereas position
            # estimation is done with an iterative least squares.
            
            # In the final iteration:
            if i == (iters-1):
                
                # Setup least squares equations to solve for Doppler.
                if gpsdop != 0.0:
                    
                    bv_mat.append(gpsdop) # Add the pseudorange rate
                    gv_mat.append(np.append(gpsvel,0.0)) # GPS velocity
                    Av_mat.append(partderiv) # Add this SV into Jacobian
        
        # Convert the matrices and vectors into numpy arrays
        A_mat = np.array(A_mat)
        bp_mat = np.array(bp_mat)
        
        # Now we get the solution for this delta...
        PosSln = np.linalg.lstsq(A_mat,bp_mat,rcond=None)
        
        # Adjust the new position values
        Xp = Xp + PosSln[0][0] # Update coordinate X
        Yp = Yp + PosSln[0][1] # Update coordinate Y
        Zp = Zp + PosSln[0][2] # Update coordinate Z
        Tp = Tp + PosSln[0][3] # Update receiver clock bias
        
        # Update the LEO position vector
        leopos = np.array([Xp,Yp,Zp]) # Initial position vector
    
    ''' Now that we have the positions and the geometry matrix...
    ... we are able to determine the velocity of the LEO via Doppler '''
    
    bv_mat = np.array(bv_mat)
    gv_mat = np.array(gv_mat)
    Av_mat = np.array(Av_mat)
    bv2_mat = np.array([]) 
    
    # Dot product of GPS satellite velocities and the Jacobian matrix
    for q in range(0,len(bv_mat)):
        bv2_elem = np.dot(gv_mat[q],Av_mat[q])
        bv2_mat = np.append(bv2_mat,bv2_elem)
    
    # Update the combined pseudorange rate residuals vector
    bv_mat = bv_mat + bv2_mat
    
    # Solve for the velocities now
    VelSln = np.linalg.lstsq(Av_mat,bv_mat,rcond=None)
    
    # Adjust the new velocity values
    Xv = VelSln[0][0] # Update coordinate X
    Yv = VelSln[0][1] # Update coordinate Y
    Zv = VelSln[0][2] # Update coordinate Z
    Tv = VelSln[0][3] # Update time drift T
    
    # Update the LEO velocity vector
    leovel = np.array([Xv,Yv,Zv]) # Initial position vector
    
    ''' Then, we shift the final solution by the antenna offsets. '''
    
    Xp = Xp - offsetX # Update coordinate X
    Yp = Yp - offsetY # Update coordinate Y
    Zp = Zp - offsetZ # Update coordinate Z
    
    ''' We would also like to calculate the PDOP, TDOP, and GDOP values '''
    
    H_mat = np.linalg.inv(np.matmul(A_mat.transpose(),A_mat))
    
    H11 = H_mat[0][0] # X component of H matrix
    H22 = H_mat[1][1] # Y component of H matrix
    H33 = H_mat[2][2] # Z component of H matrix
    H44 = H_mat[3][3] # T component of H matrix
    
    GDOP = (H11 + H22 + H33 + H44)**0.5
    PDOP = (H11 + H22 + H33)**0.5
    TDOP = H44**0.5
    
    # Check for invalid DOP values (too high values)
    if GDOP > 100 or PDOP > 100 or TDOP > 100:
        print('Extremely high DOP values found in ' + str(epoch))
        print('DOP values will be reset to zero (considered invalid)')
        print('GDOP = '+str(GDOP)+' PDOP = '+str(PDOP)+' TDOP = '+str(TDOP))
        print('\n')
        GDOP, PDOP, TDOP = 0.0, 0.0, 0.0
        Xp, Yp, Zp, Tp = 0.0, 0.0, 0.0, 0.0
        Xv, Yv, Zv, Tv = 0.0, 0.0, 0.0, 0.0
    
    # Check for invalid DOP values (NaN values computed somehow)
    if np.isnan(GDOP) or np.isnan(PDOP) or np.isnan(TDOP):
        print('GDOP values are computed to be NaN in ' + str(epoch))
        print('DOP values will be reset to zero (considered invalid)')
        print('GDOP = '+str(GDOP)+' PDOP = '+str(PDOP)+' TDOP = '+str(TDOP))
        print('\n')
        GDOP, PDOP, TDOP = 0.0, 0.0, 0.0
        Xp, Yp, Zp, Tp = 0.0, 0.0, 0.0, 0.0
        Xv, Yv, Zv, Tv = 0.0, 0.0, 0.0, 0.0
        
    # Check if the LEO is "buried in Earth" or "above the GPS constellation"...
    if np.linalg.norm(np.array([Xp,Yp,Zp])) > 6371140 + 20000000:
        print('LEO altitude is above the GPS constellation in ' + str(epoch))
        print('LEO computed positions to be reset to zero (invalid). \n')
        Xp, Yp, Zp, Tp = 0.0, 0.0, 0.0, 0.0
        Xv, Yv, Zv, Tv = 0.0, 0.0, 0.0, 0.0
    
    elif np.linalg.norm(np.array([Xp,Yp,Zp])) < 6371140:
        print('LEO altitude is below the Earth surface in ' + str(epoch))
        print('LEO computed positions to be reset to zero (invalid). \n')
        Xp, Yp, Zp, Tp = 0.0, 0.0, 0.0, 0.0
        Xv, Yv, Zv, Tv = 0.0, 0.0, 0.0, 0.0
    
    ''' Finally, we output the (hopefully) converged solution... '''
    
    posf = np.array([Xp,Yp,Zp,Tp])
    velf = np.array([Xv,Yv,Zv,Tv])
    dopf = np.array([GDOP,PDOP,TDOP])
    clkb = np.array([Tp])
    
    return posf, velf, dopf, clkb
