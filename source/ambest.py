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
## FILE DESCRIPTION:                                                         ##
##                                                                           ##
## This module contains functions that support epoch-wise processing for     ##
## integer ambiguity resolution step. The chief function in the module is    ##
## the 'ambest' function, which outputs the precise relative baseline vector ##
## between the two spacecraft. This processing is done only for one epoch,   ##
## and thus it is a snapshot-style of baseline estimation rather than a      ##
## sequential or batch of carrier phase observations. All other supporting   ##
## functions are called within 'ambest'. At the user-level, it is advised to ##
## modify contents only within 'ambest' unless the user wishes to modify the ##
## core integer ambiguity resolution algorithm.                              ##
##                                                                           ##
## INPUTS:                                                                   ##
##                                                                           ##
## epoch = datetime.datetime object                                          ##
## gps = {1:{'px':123, 'py':123, 'pz':123,...}, ... 32:(...)}                ##
## rx1 = {SV:{'C1':123,'L1':123,'L4':321,'flag':'none'},...}                 ##
## rx2 = {SV:{'C1':123,'L1':123,'L4':321,'flag':'none'},...}                 ##
## pos1 = [PosX,PosY,PosZ,Bias] of LEO satellite 1                           ##
## pos2 = [PosX,PosY,PosZ,Bias] of LEO satellite 2                           ##
##                                                                           ##
##                                                                           ##
## OUTPUT:                                                                   ##
##                                                                           ##
## Precisely determined relative baseline between the two LEOs in formation. ##
##                                                                           ##
## REMARKS:                                                                  ##
##                                                                           ##
## This is really the key ingredient in LEOGPS for precise determination of  ##
## baselines. However, as of the current version v1.2 Alpha, it is unable to ##
## find a suitable method that fixes the floats to integers as a solution,   ##
## and thus uses the pseudorange as an estimate for the float solution only. ##
## This part of the code and algorithm is still under development.           ##
##                                                                           ##
## AUTHOR MODIFIED: 30-May-2021, by Samuel Y.W. Low                          ##
##                                                                           ##
###############################################################################
###############################################################################


import numpy as np

# IMPORT LOCAL LIBRARIES
from source import consts
from source import azimel
from source import ambfix

# Define the main ambiguity estimation function below

def ambest( epoch, gps, rx1, rx2, pos1, pos2, inps, sigma=0.002, fix=False ):
    
    # Before anything else, check if the position vectors are valid.
    if np.sum(pos1) == 0.0 or np.sum(pos2) == 0.0:
        return np.array([0.0, 0.0, 0.0])
    
    # Get all input variables first.
    freqnum = inps['freq'] # Single or dual-frequency?
    C = consts.C # Speed of light
    F1 = consts.FREQ1 # L1 frequency 1575420000 Hz
    F2 = consts.FREQ2 # L2 frequency 1227600000 Hz
    
    if freqnum == 1:
        wavelength = C/F1 # Normal L1 wavelength
    elif freqnum == 2:
        wavelength = C/(F1-F2) # Wide laned wavelength
    else:
        print('The input frequency mode is invalid! It is either 1 or 2.')
        return None
    
    # The float and integer estimation of the ambiguities happens here.
    
    # Numpy-rize the position vectors
    pos1 = np.array(pos1)
    pos2 = np.array(pos2)       
    
    # Get list of observable GPS satellites, with local function getgps.
    glist = getgps(gps, rx1, rx2)
    
    # If there are not enough common satellites between rx1 and rx2...
    if len(glist) < 4:
        print('Error in baseline estimation at ' + str(epoch))
        print('Not enough common satellites for double differencing! \n')
        return np.array([0.0,0.0,0.0])
    
    # Get the double difference reference satellite SV number.
    grefr = getref(gps, rx1, rx2, pos1, pos2)
    
    # Get the index of the double difference reference satellite.
    grefi = glist.index(grefr)
    
    # Get the observations (code phase, carrier phase, from both LEO 1 & 2).
    obs = getobs( freqnum, gps, rx1, rx2 )
    
    # Get the DD transformation matrix.
    diff = difference( len(glist), grefi )
    
    # Form the zero-difference code and carrier observations of both LEOs
    codeZD = np.concatenate((np.array(obs['LEO1 Code']),
                             np.array(obs['LEO2 Code'])))
    carrZD = np.concatenate((np.array(obs['LEO1 Carr']),
                             np.array(obs['LEO2 Carr'])))
        
    # Form the double-difference code and carrier observations of both LEOs
    codeDD = np.matmul(diff,codeZD) # Converts ZD values into DD values
    carrDD = np.matmul(diff,carrZD) # Converts ZD values into DD values
    
    # Now, we can estimate the double difference float ambiguity estimates!
    ahat = np.array([]) # Initialise an empty array first.
    for x in range(0,len(codeDD)): # Parse through the elements.
        estimate = carrDD[x] - ( codeDD[x] / wavelength ) # Estimate ambiguity
        ahat = np.append(ahat,estimate) # Append the array
    
    # Check if the user has opted to fix the ambiguities using LAMBA.
    if fix == False:
        zfix = ahat # Choose the float ambiguities, no fixing by default.
    
    # Else if true, perform LAMBDA integer ambiguity resolution.
    elif fix == True:
        # Now, our final step is to estimate the covariance matrix of DD data.
        # Currently, the covariance matrix is just a scaled identity matrix.
        # You will have to change this manually as of v1.2, in this code.
        covZD = sigma*np.identity(2*len(glist))
        covDD = np.matmul(np.matmul(diff,covZD),diff.transpose())
        Qahat = covDD
        afix, sqnorm = ambfix.LAMBDA(ahat, Qahat)
        zfix = afix[0] # Choose the best ambiguities.
    
    # Else, if non-binary input, print a warning message, output zero vector.
    else:
        print('Error! Ambiguity fixing option is invalid (neither T/F!) \n')
        return np.zeros(3)
    
    ###########################################################################
    ###########################################################################
    
    # We now proceed to solve the least squares problem.
    A = getgeom( gps, glist, grefr, grefi, pos1, pos2, inps ) # Geometry matrix
    B = wavelength * (carrDD - zfix) # Carrier phase minus integer ambiguities
    
    # Solve the equation now!
    X = np.linalg.lstsq(A,B,rcond=None)
    
    # To prevent bugs, check if there is a 'NaN' value recorded.
    if np.isnan(np.sum(X[0])):
        print('WARNING! Baseline estimate NaN entry recorded in '+str(epoch))
        print('Check the RINEX and input files for any errors! \n')
        return np.array([0.0, 0.0, 0.0])
        
    # If it was processed in batch-mode, we could record multiple epochs.
    # We could then find the mean of the float ambiguities.
    # We could also then calculate (empirically) the covariances.
    
    return X[0] # Return the corrected relative positions.

# Define a function that outputs the observed GPS SVs by both LEOs

def getgps( gps, rx1, rx2 ):

    Glist = [] # Ordered array of observed GPS satellites.
    
    for k in gps:
                    
        # Check if k is observed in both RINEX files?
        if k in rx1 and k in rx2:
            
            # We also have to screen for cycle slips.
            flag1 = rx1[k]['flag']
            flag2 = rx2[k]['flag']
            
            if flag1 != 'slip' and flag2 != 'slip':
            
                # Add this SV to the list of GPS satellites.
                Glist.append(k)
    
    Glist.sort() # Sort it, just in case...
    
    return Glist # Ordered array of observed GPS satellites.

# Define a function that finds a DD reference based on elevation

def getref( gps, rx1, rx2, pos1, pos2 ):
    
    # This function selects the reference GPS satellite based on elevation.
    # Input variables are all epoch-specific observations.
    
    leopos = (0.5 * (pos1 + pos2))[:3] # Formation center
    G_elev = {} # Dictionary variable to hold elevation angle information.
    G_list = getgps(gps, rx1, rx2) # GPS satellites
        
    for k in G_list:

        gpspos = np.array([gps[k]['px'], gps[k]['py'], gps[k]['pz']]) * 1000
        G_elev[k] = azimel.azimel(leopos,gpspos)[1] # In radians

    # Reference GPS satellite SV ID is chosen based on highest elevation.

    Grefr = max(G_elev,key=G_elev.get)
    
    return Grefr

# Define a function that finds the DD geometry matrix
    
def getgeom( gps, glist, grefr, grefi, pos1, pos2, inps ):
    
    DD_geom_mat = np.array([]) # Matrix for DD geometry vector.
    leopos = (0.5 * (pos1 + pos2))[:3] # Formation center
    erp = inps['earthrotation'] # Offset Earth rotation?
    
    for k in glist:
        
        # Get the double-differenced geometry vector.
        if k != grefr:
            
            # Get the GPS satellite positions
            gpsk = np.array([gps[k]['px'],
                             gps[k]['py'],
                             gps[k]['pz']]) * 1000
            gpsr = np.array([gps[grefr]['px'],
                             gps[grefr]['py'],
                             gps[grefr]['pz']]) * 1000
            
            # Get the GPS satellite velocities
            gpskv = np.array([gps[k]['vx'],
                              gps[k]['vy'],
                              gps[k]['vz']]) * 1000
            gpsrv = np.array([gps[grefr]['vx'],
                              gps[grefr]['vy'],
                              gps[grefr]['vz']]) * 1000
            
            # Now, we need to compensate for the GPS satellite positions due
            # to Earth's rotation and the movement of the GPS satellites
            # during the signal time-of-flight.
            if erp == 'True':
                
                # Get the signal time of flights for the k-th GPS SV
                sigTOF_k = np.linalg.norm( gpsk - leopos ) / consts.C
                
                # Get the signal time of flights for the k-th GPS SV
                sigTOF_r = np.linalg.norm( gpsr - leopos ) / consts.C
                
                # Rotation rate * TOF
                Rk = consts.ERR * sigTOF_k
                Rr = consts.ERR * sigTOF_r
                
                # Get the Earth rotation direction cosine matrices,
                # using small angle approximations for simplicity.
                dcm_k = np.array([[  1.0,   Rk, 0.0 ],
                                  [ -1*Rk, 1.0, 0.0 ],
                                  [  0.0,  0.0, 1.0 ]])
                dcm_r = np.array([[  1.0,   Rr, 0.0 ],
                                  [ -1*Rr, 1.0, 0.0 ],
                                  [  0.0,  0.0, 1.0 ]])
                
                # Rotate the coordinate frame for GPS satellite positions
                gpsk = np.matmul(dcm_k, gpsk)
                gpsr = np.matmul(dcm_r, gpsr)
                
                # Compensate for GPS satellite motion during signal TOF
                gpsk = gpsk - (gpskv * sigTOF_k)
                gpsr = gpsr - (gpsrv * sigTOF_r)
            
            # First, obtain the undifferenced relative vector of nth GPS.
            ZD_geom1 = gpsk - leopos
            ZD_geom1 = ZD_geom1 / np.linalg.norm(ZD_geom1)
            
            # Next, get undifferenced relative vector of reference kth GPS.
            ZD_geom2 = gpsr - leopos
            ZD_geom2 = ZD_geom2 / np.linalg.norm(ZD_geom2)
            
            # Calculate the geometry vector used in DD least squares.
            DD_geom_vec = ZD_geom2 - ZD_geom1
            
            # Add this geometry vector to the geometry array.
            DD_geom_mat = np.append( DD_geom_mat, DD_geom_vec )
    
    # Reshape the geometry matrix so that it has 3 columns for XYZ.
    DD_geom_mat = DD_geom_mat.reshape((len(glist)-1,3))
    
    return DD_geom_mat

# Define a function that extracts code and carrier observations

def getobs( freqnum, gps, rx1, rx2 ):
    
    # Output is a dictionary with 4 keys:
    # obs['LEO1 Code'] = Code pseudorange observations of LEO1
    # obs['LEO1 Carr'] = Carrier phase observations of LEO1
    # obs['LEO2 Code'] = Code pseudorange observations of LEO2
    # obs['LEO2 Carr'] = Carrier phase observations of LEO2
    
    obs = {} # Initialise a dictionary
    F1 = consts.FREQ1 # L1 frequency 1575420000 Hz
    F2 = consts.FREQ2 # L2 frequency 1227600000 Hz
    LEO1_Code = [] # Code pseudorange observations of LEO1
    LEO1_Carr = [] # Carrier phase observations of LEO1
    LEO2_Code = [] # Code pseudorange observations of LEO2
    LEO2_Carr = [] # Carrier phase observations of LEO2
    
    # First, we find all SVs common to both RINEX observations.
    glist = getgps(gps, rx1, rx2)
    
    # Consider the case of L1 only...
    
    if freqnum == 1:
        
        for k in glist:
            
            # Get the code phase values
            
            if 'P1' in rx1[k] and 'P1' in rx2[k]:
                LEO1_Code.append(rx1[k]['P1'])
                LEO2_Code.append(rx2[k]['P1'])
            elif 'C1' in rx1[k] and 'C1' in rx2[k]:
                LEO1_Code.append(rx1[k]['C1'])
                LEO2_Code.append(rx2[k]['C1'])
            else:
                print('Error! Neither C1 nor P1 values in RINEX!')
                return None
            
            # Then get the carrier phase values
            
            LEO1_Carr.append(rx1[k]['L1'])
            LEO2_Carr.append(rx2[k]['L1'])
    
    # If we have L1/L2, we can employ wide-lane float ambiguity resolution.
    
    elif freqnum == 2:
        
        for k in glist:
            
            # Get the code phase values
            
            if 'P1' in rx1[k] and 'P1' in rx2[k]:
                R1 = (F1*rx1[k]['P1']+F2*rx1[k]['P2'])/(F1+F2)
                LEO1_Code.append(R1) # Narrow lane pseudorange
                R2 = (F1*rx2[k]['P1']+F2*rx2[k]['P2'])/(F1+F2)
                LEO2_Code.append(R2) # Narrow lane pseudorange
            elif 'C1' in rx1[k] and 'C1' in rx2[k]:
                R1 = (F1*rx1[k]['C1']+F2*rx1[k]['P2'])/(F1+F2)
                LEO1_Code.append(R1) # Narrow lane pseudorange
                R2 = (F1*rx2[k]['C1']+F2*rx2[k]['P2'])/(F1+F2)
                LEO2_Code.append(R2) # Narrow lane pseudorange
            else:
                print('Error! Neither C1 nor P1 values in RINEX!')
                return None
            
            # Then get the carrier phase values
            
            L1 = (F1*rx1[k]['L1']-F2*rx1[k]['L2'])/(F1-F2)
            LEO1_Carr.append(L1) # Wide lane carrier phase
            L2 = (F1*rx2[k]['L1']-F2*rx2[k]['L2'])/(F1-F2)
            LEO2_Carr.append(L2) # Wide lane carrier phase
        
    else:
        
        print('The input frequency mode is invalid! It is either 1 or 2.')
        return None
    
    obs['LEO1 Code'] = LEO1_Code
    obs['LEO1 Carr'] = LEO1_Carr
    obs['LEO2 Code'] = LEO2_Code
    obs['LEO2 Carr'] = LEO2_Carr
    
    return obs

# Define a matrix transforming 02x single to double difference values

def difference( s, r ):
    
    # Generates the transformation matrix D at a particular time t.
    # Transforms a vector of zero difference phase observables [2s x 1]
    # Into a vector of double difference phase observables [(s-1) x 1]
    # s = number of observed GPS satellites in ONE epoch
    # r = index of chosen reference satellite for double differencing
    # It is assumed that there are only 2 LEO satellites (i.e. 2 receivers)
    # It is assumed that across q epochs, s remains unchanged.
    
    # We begin by constructing the single difference transformation matrix
    # For an input of a [2s x 1] array of undifferenced observations
    
    SD_mat_LEO1 = np.identity(s) # Take the phase values of LEO1...
    SD_mat_LEO2 = -1*np.identity(s) # ... and subtract LEO2 from them.
    SD_mat = np.concatenate( ( SD_mat_LEO1 , SD_mat_LEO2 ), axis = 1 )

    DD_ins = -1*np.ones(s-1) # Inserted at index k, to subtract reference.
    DD_mat = np.insert( np.identity(s-1), r, DD_ins, axis=1) # DD_ins

    Transf = np.matmul(DD_mat, SD_mat) # Transform matrix from ZD to DD
    
    return Transf
