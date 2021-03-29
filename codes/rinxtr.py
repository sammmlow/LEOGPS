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
## Extraction of RINEX observations C1/P1, P2, L1, L2, D1, D2                ##
## If Doppler observables D1/D2 are not found in the RINEX observation file, ##
## they will be estimated by a first-order derivative of L1/L2 phase values. ##
##                                                                           ##
## INPUTS:                                                                   ##
##                                                                           ##
## RINEX observation file (v2.xx) of LEO satellite.                          ##
## Only GPS observables (no multi-GNSS) are supported until new updates      ##
##                                                                           ##
## OUTPUT:                                                                   ##
##                                                                           ##
## RINEX observations as a dictionary of epochs, each epoch with a sub       ##
## dictionary of GPS satellites based on SVIDs, and each SVID with a         ##
## sub dictionary of the various observations (C1/P1,P2,L1,L2,D1,D2).        ##
## Output = {epoch1:{5:{'C1':123,'L1':123, ... 'L4':321,'flag':'none'}...}...##
##           epoch2:{3:{'C1':123,'L1':123, ... 'L4':321,'flag':'slip'}...}...##
##           ... ... ... ... ... ...                                         ##
##           epochX:{3:{'C1':123,'L1':123, ... 'L4':321,'flag':'none'}...}}  ##
##                                                                           ##
## REMARKS:                                                                  ##
##                                                                           ##
## You may enable a code-carrier smoothing Hatch filter in the config file.  ##
## You may also change the length of the cycle slip filter, and the filter,  ##
## tolerance in terms of the number of standard deviations in the config.txt ##
## Ensure that RINEX observation files follow 4-letter ID naming convention  ##
## Followed by the DOY + 0, with the file extension .YYO                     ##
##                                                                           ##
## AUTHOR MODIFIED: 14-03-2020, by Samuel Y.W. Low                           ##
##                                                                           ##
###############################################################################
###############################################################################
'''

import datetime

# IMPORT LOCAL LIBRARIES
from codes import phasep
from codes import dopest

''' Function that begins parsing the RINEX file with the file name as input '''

def rinxtr(namepath, inps, goodsats, tstart, tstop, rnxstep):
    
    freqnum   = inps['freq'] # Single frequency or dual frequency processing?
    hatchfilt = inps['hatchfilter'] # Enable hatch filtering?
    hatchlen = -1*inps['hatchlength'] # Length of hatch filter

    # Create a meta data dictionary that holds header information
    rnxmeta = {} # {'NumObs': N, 'TypeObs': ['L1','L2','C1','P2']}
    
    # Create the main data dictionary that holds ALL the observation info
    rnxdata = {} # {epoch1:{SV1:[P1,P2,L1,L2],SV2:[...] , epoch2:{...}}}
    
    # Initialise a trigger to record observations
    trigger = False
    
    # Initialise other parameters to be used in this scenario.
    obsi = []
    
    # Initialise a parameter used for triggering the recording of observations.
    obsready = False
    
    # Initialise an integer for SV ID numbers.
    SVNum = 0
    
    # A trigger to mark that LEOGPS is reading the epoch header.
    reading_header = False
    
    # First we open the file corresponding to input 'namepath'
    with open(namepath) as file:  
        
        # Create a readlines object to parse through
        rnxlines = file.readlines()

        # Parse through the first line to check RINEX version
        if '2.' not in rnxlines[0]:
            print('This is the wrong RINEX version!')
            print('This program accepts only RINEX version 2.xx!')
            return False

        # Parse through the lines to find observation symbols
        for obsline in rnxlines:
            
            if 'TYPES OF OBSERV' in obsline:
                
                obs = obsline[:60].split()
                obsi = obsi + obs # Concatenate observations
            
            elif 'END OF HEADER' in obsline:
                
                print('Observables found in RINEX file:')
                print(obsi)
                print('\n')
                
                # Save only L1-related observations for single frequency.
                if freqnum == 1:
                    obsz = [x for x in obsi[1:] if '1' in x]
                    rnxmeta['TypeObs'] = obsz
                    rnxmeta['NumObs'] = len(obsz)
                    
                # Save L1 and L2 related observations for dual frequency.
                if freqnum == 2:
                    obsz = [x for x in obsi[1:] if '1' in x or '2' in x]
                    rnxmeta['TypeObs'] = obsz
                    rnxmeta['NumObs'] = len(obsz)
                
                rnxmeta['NumObs'] = int(obsi[0]) # Number of observations
                rnxmeta['TypeObs'] = obsi[1:] # List of observation types
                
                # Perform a sanity check; NumObs must be the length of TypeObs
                if rnxmeta['NumObs'] != len(rnxmeta['TypeObs']):
                    
                    print('Error! Failed to parse RINEX observables!')
                    print('Check that RINEX file is correctly formatted!')
                    print('Check that RINEX version is v2!')
                    print('/n')
                    return False
                
                TypeObs = rnxmeta['TypeObs']
                NumObsv = rnxmeta['NumObs']
                
                break # Break out of the for loop once RINEX header ends
        
        # Check if user selection of single/dual frequency is valid
        if (freqnum == 1):
            print('User selected single frequency processing.')
            print('\n')
            if ('L1' in TypeObs) and ('L2' in TypeObs):
                print('But RINEX observation contains L1 and L2 frequencies!')
                print('LEOGPS will proceed with processing L1-only data!')
                print('\n')
            if ('D1' not in TypeObs):
                print('Doppler not found in RINEX observation.')
                print('Pseudorange rates will have to be coarsely estimated.')
                print('LEO velocity will also have to be coarsely estimated.')
                print('\n')
                DopplerCheck = False # This will be used in posvel.py
            else:
                print('Doppler found in RINEX observation')
                print('LEO satellite velocity can also be estimated.')
                print('\n')
                DopplerCheck = True # This will be used in posvel.py
            
        elif (freqnum == 2):
            print('User selected dual frequency processing.')
            if ('L1' not in TypeObs) or ('L2' not in TypeObs):
                print('However, RINEX observation contains only 1 frequency!')
                print('Update config.txt freqnum parameter.')
                print('\n')
                return False
            if ('D1' not in TypeObs) or ('D2' not in TypeObs):
                print('Doppler not found in RINEX observation.')
                print('Pseudorange rates will have to be coarsely estimated.')
                print('LEO velocity will also have to be coarsely estimated.')
                print('\n')
                DopplerCheck = False # This will be used in posvel.py
            else:
                print('Doppler found in RINEX observation')
                print('LEO satellite velocity can also be estimated.')
                print('\n')
                DopplerCheck = True # This will be used in posvel.py
                
        else:
            print('Invalid input for freqnum in *config.txt*!')
            return False
        
        # Check for where the end of the RINEX header is located at.
        header_reached = False
        header_linenum = 0
        while header_reached == False:
            headerline = rnxlines[header_linenum]
            if 'END OF HEADER' in headerline:
                header_reached = True
            header_linenum += 1
        
        for data in rnxlines[header_linenum:]:
            
            # Replace the G prefix with a space.
            data = data.replace('G',' ') # GPS satellites
            dataspl = data.split()
            
            # READING OF DATES AND TIMES IN THIS SEGMENT.      
            # We use the number of decimal points '.' in one line,
            # To distinguish whether it is an observation line, or time entry.
            
            NumDecs = sum([m.count('.') for m in dataspl]) # Count decimals
        
            # The entry should be a date-time entry if there are > 7 elements
            if len(dataspl) >= 8 and NumDecs <= 2:
                
                if 'E' in data or 'B' in data or 'R' in data:
                    print('GLONASS, GALILEO or BEIDOU satellites observed!')
                    print('Unable to process, only GPS allowed in LEOGPS!')
                    return False
                
                epoch = datetime.datetime(int(dataspl[0]) + 2000,
                                          int(dataspl[1]),
                                          int(dataspl[2]),
                                          int(dataspl[3]),
                                          int(dataspl[4]),
                                          int(round(float(dataspl[5]))))
                
                # Now, let's clean up the previous epoch, in case there were
                # still float-zero values due to corrupted RINEX data.
                if epoch-rnxstep in rnxdata:
                    if obsready == False:
                        if SVNum in rnxdata[epoch-rnxstep]:
                            del rnxdata[epoch-rnxstep][SVNum]
                
                # Now, we check if observations are still within time frame.
                if epoch >= tstart and epoch <= tstop:
                    
                    # Initialise some flags and recording triggers.
                    reading_header = True # Mark that we're reading the epoch
                    trigger = True # Start recording RINEX observations
                    rnxdata[epoch] = {} # Save this time stamp in 'rnxdata'
                    
                    # Initialise some useful counters
                    satcount = 0 # To be used later...
                    obscount = 0 # To be used later...
                    
                    # Check if this entry's RINEX observation marker is OK
                    if dataspl[6] != '0':
                        
                        print('Warning! Potentially problematic epoch at:')
                        print(str(epoch) + '\n')
                        
                    # Setup the list of observable GPS satellites
                    SVList = [int(x) for x in dataspl[8:]]
                        
                    # Now populate this epoch with GPS observations
                    for SV in SVList:
                        rnxdata[epoch][SV] = {}
                        for obsv in TypeObs:
                            rnxdata[epoch][SV][obsv] = 0.0
                        
                    # Prepare an empty placeholder for observations
                    data_obs = []
                    
                else:
                    trigger = False # Stop recording RINEX observations
                
                continue # Do not remove this line!
            
            # Check if the epoch header has more than one line...?
            if reading_header == True:
                
                # Check if this is a continuation of the previous line.
                if '                                ' == data[:32]:
                    
                    # Setup the list of observable GPS satellites
                    SVList2 = [int(x) for x in dataspl]
                        
                    # Now populate this epoch with GPS observations
                    for SV in SVList2:
                        SVList.append(SV)
                        rnxdata[epoch][SV] = {}
                        for obsv in TypeObs:
                            rnxdata[epoch][SV][obsv] = 0.0
                
                    continue # Do not remove this line!
                    
            # READING OF OBSERVATION VALUES IN THIS SEGMENT.
            # Now, check if this is a valid observation line

            if trigger == True:
                
                reading_header = False # Switch off the epoch-header flag.
                SVNum = SVList[satcount] # What SV ID?
                obscount_now = 0 # Present number of observations
                
                # Parse through the data array and discard SNR / LLI info
                for obs in dataspl:
                     
                    if len(obs) > 2:
                        obscount += 1 # Cumulative number of observations
                        obscount_now += 1 # Present number of observations
                        decpoint = obs.find('.') # Find decimal point
                        obsv = float(obs[:decpoint+4]) # Observation value                      
                        data_obs.append(obsv)
                
                # RINEX observations will span multiple lines.
                # ... this portion of the code ensures the data integrity,
                # ... pertaining to each SV ID's number of observations,
                # ... if there are less than 'NumObsv' observations, then,
                # ... this GPS satellite's data is corrupt and will be erased.
                
                # Has observation count matched 'NumObsv' derived from header?
                if obscount < NumObsv: # Not finished reading observables.
                    obsready = False # Not ready to add data into 'rnxdata'.
                elif obscount == NumObsv:
                    obsready = True # Finished reading observables; record now.
                else: # We have read more than we should.
                    data_obs = data_obs[(obscount-obscount_now):]
                    obscount = obscount_now
                    SVList.remove(SVNum)
                    del rnxdata[epoch][SVNum] # Remove entry
                    SVNum = SVList[satcount] # Update the new SV ID
                    
                    # Now, we must re-check again if obscount matches NumObsv.
                    if obscount < NumObsv: # Not finished reading observables.
                        obsready = False # Not ready to add data.
                    elif obscount == NumObsv:
                        obsready = True # Finished reading observables.

                if obsready == True:
                    
                    # Now, populate rnxdata with non-zero GPS observations
                    for m in range(0,len(data_obs)):
                        
                        ObsT = TypeObs[m] # What observation?
                        ObsV = data_obs[m] # What value of observation?
                        rnxdata[epoch][SVNum][ObsT] = ObsV
                    
                    # If that was a bad satellite, filter it out!
                    if SVNum not in goodsats:
                        del rnxdata[epoch][SVNum]
                        
                    data_obs = [] # Reset the placeholder
                    obscount = 0 # Reset the observation counter
                    satcount += 1 # Update the i-th GPS satellite
                    
                    # If all satellites are accounted for,
                    # Variable 'trigger' can be switched off                    
                    if satcount == len(SVList):
                        trigger = False

    # Now, our final product is the dictionary of RINEX data: 'rnxdata'
    # We would like to parse through this dictionary and screen for cycle slips
    # However, we have to mark our carrier phase readings first with flags.
    rnxmark = phasep.phsmrk(rnxdata, rnxstep, goodsats, inps)

    # After flagging our data, we can check if Doppler values exist.
    if DopplerCheck == False:
        
        # If it doesn't, then estimate them.
        rnxmark = dopest.dopest(rnxmark,goodsats,tstart,tstop,rnxstep,inps)
        
        # Check if any errors occured during the Doppler estimation.
        if rnxmark == False:
            return False
    
    # Check if the user desires to perform hatch filtering.
    if hatchfilt == 'True':
        
        # The file containing phase processing routines is 'phasep.py'
        # Within it, contains several routines, we will call 3 main routines:
        # phasep.phsmrk(...) ---> marks data flags = start/end/none/slip
        # phasep.ph1fil(...) ---> single-frequency Hatch filtering
        # phasep.ph2fil(...) ---> dual-frequency divergence-free filtering
        
        # We will perform hatch filtering on the RINEX data now
        if freqnum == 1:
            rnxproc = phasep.ph1fil(rnxmark, rnxstep, goodsats, hatchlen)
        
        elif freqnum == 2:
            rnxproc = phasep.ph2fil(rnxmark, rnxstep, goodsats, hatchlen)
            
        return rnxproc
    
    else:
        return rnxmark
