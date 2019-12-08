#!/usr/bin/env python3
'''
###############################################################################
###############################################################################
##                                                                           ##
##     _    ___  ___   ___ ___ ___                                           ##
##    | |  | __ /   \ / __| _ | __|                                          ##
##    | |__| __  ( ) | (_ |  _|__ \                                          ##
##    |____|___ \___/ \___|_| \___/                                          ##
##                                    v 0.1 (Alpha)                          ##
##                                                                           ##
## FILE DESCRIPTION:                                                         ##
##                                                                           ##
## Code-based position estimation via weighted least squares estimation.     ##
## Doppler-based estimation of velocities only if Doppler is available.      ##
## Ionospheric delay for L1: Klobuchar Model (NOT IMPLEMENTED YET!).         ##
## Ionospheric delay for L2: Dual-frequency iono-free linear combination.    ##
## Clock sync done by extracting interpolating clock biases and drifts.      ##
## Relativistic effects estimated: Shapiro effect + clock delay/advance.     ##
## This script 'posest.py' is called up for every epoch in 'mainpp.py'       ##
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
## AUTHOR MODIFIED: 30-11-2019, by Samuel Y.W. Low                           ##
##                                                                           ##
###############################################################################
###############################################################################
'''

import numpy as np

# IMPORT LOCAL LIBRARIES
from codes import consts
from codes import einstn

''' This is a position-velocity estimator for the SINGLE-EPOCH case. '''

def posvel(epoch, goodsats, gps, rxi, inps, iters = 6):
    
    # gps = {1:{'px':123, 'py':123, 'pz':123,...}, ... 32:(...)}
    # rxi = {1:('C1':xxx,'P2':xxx,'L1':xxx,'L2':xxx) ... 32:(...)}
    
    ''' Before starting, we need at least 4 GPS satellites to solve XYZ! '''
    if len(rxi) < 4:
        print('Error! For the epoch of ' + str(epoch))
        print('Only ' + str(len(rxi)) + ' GPS satellites observed!')
        print('Unable to obtain a solution fix, will output float zeros!')
        posf = np.array([0.0, 0.0, 0.0, 0.0])
        velf = np.array([0.0, 0.0, 0.0, 0.0])
        dopf = np.array([0.0, 0.0, 0.0, 0.0])
        return posf, velf, dopf
    
    ''' First, initialise variables and constants '''
    
    C = consts.C # Speed of light in m/s
    F1 = consts.FREQ1 # L1 frequency 1575420000 Hz
    F2 = consts.FREQ2 # L1 frequency 1227600000 Hz
    WL1 = C/F1 # L1 wavelength
    WL2 = C/F2 # L2 wavelength
    err = consts.ERR # Earth rotation rate
    
    ''' Now, extract inputs from *config.txt* '''
    
    rel     = inps['relativity'] # Enable relativistic correction?
    erp     = inps['earthrotation'] # Offset Earth rotation?
    freqnum = inps['freq'] # Check if single or dual frequency
    offsetX = float(inps['antoffsetX']) # Antenna offset X
    offsetY = float(inps['antoffsetY']) # Antenna offset Y
    offsetZ = float(inps['antoffsetZ']) # Antenna offset Z
    
    ''' Initialise the core parameters used for estimation later on. '''
    
    Xp,Yp,Zp,Tp = 1.0, 1.0, 1.0, 1.0 # Initialise XYZT, for position.
            
    ''' Begin iterations of gradient descent for least squares estimate. '''
    
    for i in range(iters):
        
        A_mat = [] # Geometry or design matrix
        bp_mat = [] # Position-time vectrix
        
        if i == 0:
            bv_mat = [] # Pseudorange rate vectrix of LEO
            gv_mat = [] # Velocity vectrix of GPS satellite
            Av_mat = [] # Geometry matrix for velocity estimation
        
        # Begin parsing ephemeris, clock data, and observables for all PRNs
        for p in rxi:
            
            ''' Extract ephemeris and clock data from gps '''
            
            # Get the positions of this k-th GPS satellite (m)
            gpspos = np.array([gps[p]['px'],
                               gps[p]['py'],
                               gps[p]['pz']])
            
            # Get the velocity of this k-th GPS satellite (m/s)
            gpsvel = np.array([gps[p]['vx'],
                               gps[p]['vy'],
                               gps[p]['vz']])
            
            # Get the clock biases and drifts
            clkb = C * gps[p]['clkb']
            clkd = C * gps[p]['clkd']
            
            ''' Now, extract the RINEX observables from rxi '''
            
            # Are we stuck with old fashioned single frequency data?
            if freqnum == 1:
                
                # Checking for pseudorange information...
                if 'P1' in rxi[p]:
                    obs_P1 = rxi[p]['P1'] # P1 ...
                elif 'C1' in rxi[p]:
                    obs_P1 = rxi[p]['C1'] # ... or C1
                else:
                    print('C1 or P1 pseudorange is missing this epoch.')
                    print(str(epoch) + ' for PRN ID ' + str(p))
                    return False         
                                
                # Finally, checking for doppler frequency information...
                if 'D1' in rxi[p]:
                    if rxi[p]['D1'] != 'NaN':
                        obs_D1 = (-1)*rxi[p]['D1'] # D1
                    else:
                        obs_D1 = 0 # Invalid D1
                else:
                    print('D1 doppler values are missing this epoch.')
                    print(str(epoch) + ' for PRN ID ' + str(p))
                    return False
                
                # We would need to do an ionospheric correction here for F1
                # (TO BE UPDATED!)
                
                # Let's get the time-bias corrected pseudorange and rates.
                pr_rnx = obs_P1 + clkb
                pv_rnx = (obs_D1 * WL1) + clkd
                
            # Are we using dual frequency data?
            if freqnum == 2:
                
                # Checking for pseudorange information...
                if 'P1' in rxi[p] and 'P2' in rxi[p]:
                    obs_P1 = rxi[p]['P1'] # P1
                    obs_P2 = rxi[p]['P2'] # P2
                elif 'C1' in rxi[p] and 'P2' in rxi[p]:
                    obs_P1 = rxi[p]['C1'] # C1
                    obs_P2 = rxi[p]['P2'] # P2
                else:
                    print('Dual frequency P1/P2 values are missing.')
                    print(str(epoch) + ' for PRN ID ' + str(p))
                    return False
                
                # Finally, checking for doppler frequency information...
                if 'D1' in rxi[p] and 'D2' in rxi[p]:
                    if rxi[p]['D1'] != 'NaN' and rxi[p]['D2'] != 'NaN':
                        obs_D1 = (-1)*rxi[p]['D1'] # D1
                        obs_D2 = (-1)*rxi[p]['D2'] # D2
                    else:
                        obs_D1 = 0 # Invalid D1
                        obs_D2 = 0 # Invalid D2
                    
                else:
                    print('Dual frequency D1/D2 values are missing.')
                    print(str(epoch) + ' for PRN ID ' + str(p))
                    return False
                
                # Let's get the time-bias corrected pseudorange and rates.
                pr_rnx = (2.546*obs_P1 - 1.546*obs_P2) + clkb
                pv_rnx = (2.546*obs_D1*WL1 - 1.546*obs_D2*WL2) + clkd
                
            else:
                print('The input text file has an invalid frequency!')
                return False
                
            # Now, adjust the GPS positions due to Earth's rotation
            if erp == 'True':
                
                # Signal transit time from GPS satellite
                sigT = pr_rnx / C 
                
                # Get the Earth rotation direction cosine matrix
                e_rate = np.array([[ 1.0,            err*sigT,      0.0 ],
                                   [ (-1)*err*sigT,  1.0,           0.0 ],
                                   [ 0.0,            0.0,           1.0 ]])
                
                # Rotate the coordinate frame for GPS satellite position
                gpspos = np.matmul(e_rate, gpspos)

            # Relative position vector between GPS positions & LEO satellite
            leopos = np.array([Xp,Yp,Zp]) # Initial position vector
            delpos = leopos - gpspos # Relative position vector
            
            # Calculate the relativistic effects if necessary
            if rel == 'True':
                shapiro = einstn.shapiro(leopos,gpspos) # Shapiro delay
                clockadv = einstn.clockadv(gpspos,gpsvel) # Clock advance
            else:
                shapiro = 0.0
                clockadv = 0.0

            # Calculate the computed pseudorange (that is meant to converge)
            pr_comp = np.linalg.norm(delpos) - shapiro - clockadv + Tp
            pr_resd = pr_rnx - pr_comp
            
            # Calculate partial derivatives of the i-th row in the Jacobian
            partderiv = delpos / pr_comp
            partderiv = np.append(partderiv,1)
            
            # Construct the geometry matrix A and pseudorange residual 'b'
            A_mat.append(partderiv) # Add this PRN ID entry into Jacobian
            bp_mat.append(pr_resd) # Add this PRN ID entry into residual vector
            
            # Update the pseudorange rate vector and matrices.
            # We only need to do it for iter = 0, because velocity estimation
            # can be done with a directly least-squares, whereas position
            # estimation is done with an iterative least squares.
            
            if i == 0:
                if obs_D1 != 0:
                    bv_mat.append(pv_rnx)
                    gv_mat.append(np.append(gpsvel,0.0))
                    Av_mat.append(partderiv)
            
        # Convert the matrices and vectors into numpy arrays
        A_mat = np.array(A_mat)
        bp_mat = np.array(bp_mat)

        # Now we get the solution for this delta...
        PosSln = np.linalg.lstsq(A_mat,bp_mat,rcond=None)
        
        # Adjust the new position values
        Xp = Xp + PosSln[0][0] - offsetX # Update coordinate X
        Yp = Yp + PosSln[0][1] - offsetY # Update coordinate Y
        Zp = Zp + PosSln[0][2] - offsetZ # Update coordinate Z
        Tp = Tp + PosSln[0][3] # Update time offset T
    
    ''' Now that we have the positions and the geometry matrix...
    ... we are able to determine the velocity of the LEO via Doppler '''
    
    bv_mat = np.array(bv_mat)
    gv_mat = np.array(gv_mat)
    Av_mat = np.array(Av_mat)
    bv2_mat = np.array([]) 
    
    # Matrix dot product of GPS satellite velocities and the Jacobian matrix
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
    Tv = VelSln[0][3] # Update time offset T
    
    ''' We would also like to calculate the PDOP, TDOP, and GDOP values '''
    
    H_mat = np.linalg.inv(np.matmul(A_mat.transpose(),A_mat))
    
    H11 = H_mat[0][0] # X component of H matrix
    H22 = H_mat[1][1] # Y component of H matrix
    H33 = H_mat[2][2] # Z component of H matrix
    H44 = H_mat[3][3] # T component of H matrix
    
    GDOP = (H11 + H22 + H33 + H44)**0.5
    PDOP = (H11 + H22 + H33)**0.5
    TDOP = H44**0.5
    
    ''' Finally, we output the (hopefully) converged solution... '''
    
    posf = np.array([Xp,Yp,Zp,Tp])
    velf = np.array([Xv,Yv,Zv,Tv])
    dopf = np.array([GDOP,PDOP,TDOP])
    
    return posf, velf, dopf