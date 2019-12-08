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
## Extraction of position information of GPS satellites across time.         ##
## Derivation of velocity of GPS satellites using first order derivative.    ##
##                                                                           ##
## INPUTS:                                                                   ##
##                                                                           ##
## Final Product SP3 file of GPS ephemeris, and CLK_30s file of GPS.         ##
## (Auto-download from IGS CDDIS if online, and if file is not present).     ##
##                                                                           ##
## OUTPUT:                                                                   ##
##                                                                           ##
## GPS dictionary of the following nested key-value pairs:                   ##
## gpsdata = {epoch1:{1: {px:123, py:123, pz:123,                            ##
##                        vx:123, vy:123, vz:123,                            ##
##                        clkb:123, clkd:123},                               ##
##                    2: {px:123, py:123, pz:123,                            ##
##                        vx:123, vy:123, vz:123,                            ##
##                        clkb:123, clkd:123}, ...                           ##
##                         ... ... ... ... ... ...                           ##
##                    32:{px:123, py:123, pz:123,                            ##
##                        vx:123, vy:123, vz:123,                            ##
##                        clkb:123, clkd:123}} ...                           ##
##            epoch2:{1: {px:123, py:123, pz:123,                            ##
##                         ... ... ... ... ... ...}}}                        ##
##                                                                           ##
## REMARKS: Use only SP3 orbit format for GPS only (no multi-GNSS support)   ##
##                                                                           ##
## AUTHOR MODIFIED: 26-11-2019, by Samuel Y.W. Low                           ##
##                                                                           ##
###############################################################################
###############################################################################
'''

import os
import shutil
import datetime
import warnings
import subprocess
import numpy as np
import urllib.request

# IMPORT LOCAL LIBRARIES
from codes import pubplt

''' Now, this is the main routine that parses ephemeris and clock data '''

def gpsxtr(inps, tstart, tstop, tstep):
    
    # First, we would like to download GPS ephemeris and clock data    
    igsurl   = 'ftp://cddis.nasa.gov/gnss/products/' # IGS CDDIS URL
    
    # We will download them into our directories, below, later on.
    cwd = inps['cwd'] # Get current main working directory
    iwd = cwd + '\\input\\' # Get directory for ephemeris / clock files
    os.chdir(cwd) # Ensure the user is running in the main directory.
    
    # Then, we must retrieve all number of days of GPS and CLK data needed    
    days     = (tstop.date() - tstart.date()).days + 1 # Number of days
    filelist = [] # Stores all the required GPS / CLK ephemeris files
    
    # Now, check for desired ephemeris and clock files    
    for d in range(0,days):
        
        wd, wwww = gpsweekday( tstart + datetime.timedelta( days = d ) )
        name = 'igs' + wwww + wd # File name, without file extension
        ephurl = igsurl + wwww + '/' + name + '.sp3.Z' # IGS URL of ephemeris file
        clkurl = igsurl + wwww + '/' + name + '.clk_30s.Z' # IGS URL of clock file
        filelist.append(name) # Add this into the list of files for parsing
    
        # Check for SP3 ephemeris file, and download + unzip it if needed.
        if os.path.exists( iwd + name + '.sp3' ) != True:
            
            print('SP3 file for ' + name + ' not found! Attempt download now...')
            urllib.request.urlretrieve(ephurl, name + '.sp3.Z')
            print('Completed downloading the ephemeris file! Now unzipping...')
            subprocess.call(cwd+'\\gzip\\gzip.exe -d ' + name +'.sp3.Z')
            print('Files unzipped, moving them into the inputs folder.')
            shutil.move(cwd + '\\' + name + '.sp3', iwd + name + '.sp3')
            print('Unzipping completed! \n')
        
        else:
            
            print('SP3 file for ' + name + ' found! Proceeding to process!')
        
        # Check for CLK file, and download + unzip it if needed.
        if os.path.exists( iwd + name + '.clk_30s' ) != True:
            
            print('CLK file for ' + name + ' not found! Downloading now...')
            urllib.request.urlretrieve(clkurl, name + '.clk_30s.Z')
            print('Completed downloading the clock file! Now unzipping...')
            subprocess.call(cwd+'\\gzip\\gzip.exe -d ' + name +'.clk_30s.Z')
            print('Files unzipped, moving them into the inputs folder.')
            shutil.move(cwd + '\\' + name + '.clk_30s', iwd+name + '.clk_30s')
            print('Unzipping completed! \n')     
    
        else:
            
            print('CLK file for ' + name + ' found! Proceeding to process!')
    
    print('\n')
    gpsdata = {} # Main output of the program
    gpsdata['t_sp3'] = [] # Array of SP3 time values (datetime objects)
    gpsdata['t_clk'] = [] # Array of CLK time values (datetime objects)
    coords = ['x','y','z'] # List of coordinates, for use in interpolation
    goodsats = list(range(1,33)) # List of satellites without clock problems
    badsats = [] # List of satellites without clock problems
    
    for prn in goodsats: # Create an ephemeris dict for each GPS satellite
        gpsdata[prn] = {} # Initialise the main GPS ephemeris dictionary
        gpsdata[prn]['px'] = [] # Storage of X-coordinates position
        gpsdata[prn]['py'] = [] # Storage of Y-coordinates position
        gpsdata[prn]['pz'] = [] # Storage of Z-coordinates position
        gpsdata[prn]['vx'] = [] # Storage of X-coordinates velocity
        gpsdata[prn]['vy'] = [] # Storage of Y-coordinates velocity
        gpsdata[prn]['vz'] = [] # Storage of Z-coordinates velocity
        gpsdata[prn]['clkb'] = [] # Storage of clock biases
    
    for file in filelist:
        
        prnls = [] # List of PRN IDs read in each clock file
        
        ''' Start reading each downloaded GPS ephemeris file from IGS  '''
        
        f_gps = open(iwd + file + '.sp3', 'r') # Open up the SP3 file
        gps_record = False # Trigger for recording GPS data
        ti = datetime.timedelta(seconds=900) # Time step for SP3 file
        tc = datetime.timedelta(seconds=30) # Time step for CLK_30S file
        
        for line in f_gps:
            
            line = line.split() # Split up the strings
            
            # Check if the time reading is accurate...
            if line[0] == '*':
                
                timenow = datetime.datetime(int(line[1]),
                                            int(line[2]),
                                            int(line[3]),
                                            int(line[4]),
                                            int(line[5]),
                                            int(float(line[6])))
                
                if timenow not in gpsdata['t_sp3']:
                    if timenow <= tstop and timenow >= (tstart - ti):
                        gps_record = True # Trigger the recording off.
                        gpsdata['t_sp3'].append(timenow)
                    else:
                        gps_record = False # Trigger the recording off.
                        
            
            # Now, check for XYZ coordinates of each PRN ID
            if 'PG' in line[0] and gps_record == True:
                
                # Get the data we need below.
                prn = int(line[0][2:4]) # Which GPS satellite is this?
                x   = float(line[1]) # X-coordinate position
                y   = float(line[2]) # Y-coordinate position
                z   = float(line[3]) # Z-coordinate position
                
                # Assign the data into the gpsdata dictionary.
                gpsdata[prn]['px'].append(x) # X-coordinate position
                gpsdata[prn]['py'].append(y) # Y-coordinate position
                gpsdata[prn]['pz'].append(z) # Z-coordinate position
            
        f_gps.close()
        
        ''' Start reading each downloaded GPS clock (30s) file from IGS  '''
        
        f_clk = open(iwd + file + '.clk_30s', 'r') # Open up the SP3 file
        
        for line in f_clk:
            
            # First, parse out GPS satellites that have clock problems
            if 'PRN LIST' in line:
                
                line = line.replace('PRN LIST',' ') # Remove the header
                line = line.replace('G',' ') # Remove the GPS header string
                line = line.split() # Split up the strings
                line = [int(prn) for prn in line]
                
                # Now, let's parse through the bad satellites.
                prnls = prnls + line
                        
                continue
                
            # Now, let's save those clock biases and drifts.
            if 'AS G' in line:
                
                line = line.split()
                timenow = datetime.datetime(int(line[2]),
                                            int(line[3]),
                                            int(line[4]),
                                            int(line[5]),
                                            int(line[6]),
                                            int(float(line[7])))
                
                # Is the current clock readings in the user-defined time axis?
                if timenow <= tstop and timenow >= (tstart - tc):
                    
                    # If it is, then record the time now.
                    if timenow not in gpsdata['t_clk']:
                        gpsdata['t_clk'].append(timenow)
                    
                    # Then, record the clock biases for each PRN.
                    prn = int(line[1][1:])
                    gpsdata[prn]['clkb'].append(float(line[9]))
        
        f_clk.close()
        
        # Check the bad satellites with clock issues.
        for p in range(1,33):
            if p not in prnls and p not in badsats:
                badsats.append(p)
    
    # Sanity check for bad satellites and other for missing entries.
    for p in range(1,33):
        if len(gpsdata['t_clk']) != len(gpsdata[p]['clkb']):
            badsats.append(p)
            print('GPS Satellite ' + str(p) + ' with missing clock values.')
            print('Satellite will be removed from the scenario. \n')
        if p in gpsdata:
            if len(gpsdata[p]['px']) != len(gpsdata['t_sp3']):
                print('GPS Satellite ' + str(p) + ' XYZ not parsed correctly.')
            
    # Check the good satellites without clock issues.
    for p in badsats:
        if p in gpsdata:
            del gpsdata[p] # Remove the PRN from the GPS data dictionary
        if p in goodsats:
            goodsats.remove(p) # Remove the PRN from the list of good sats
    
    badsats.sort() # Sort in ascending order
    goodsats.sort() # Sort in ascending order
        
    # Sanity check: timings must match the clock data index-wise
    if len(gpsdata['t_clk']) != len(gpsdata[prn]['clkb']):
        print('Something went wrong during processing!')
        print('Length of time arrays do not match length of clock biases!')
        print('Check the format of the clock files!')
        return False
    
    ''' Begin preparing time variables for the interpolation process '''
    
    # Let us get the user-defined time axis for use in interpolation later.
    t_usr = tstep.seconds # User-defined step size in *config.txt*
    epoch = datetime.datetime(tstart.year,
                              tstart.month,
                              tstart.day)
    
    # Get the start time of the user-defined time axis, in seconds
    t_usr_i = (tstart - epoch).seconds
    
    # Get the end time of the user-defined time axis, in seconds
    t_usr_f = (tstop - epoch).days*86400 + (tstop - epoch).seconds
    
    # Finally, initiate the user-defined time axis.
    t_usr_ls = list(range(t_usr_i, t_usr_f, t_usr))
    gpsdata['t_usr'] = t_usr_ls
    
    # Let us get the time axis used in the original IGS ephemeris array
    t_eph = 900 # Step size in IGS ephemeris array, in seconds
    N_eph = len(gpsdata['t_sp3']) # No. of elements in ephemeris time array
    t_eph_i = (gpsdata['t_sp3'][0] - epoch).seconds # Start time
    t_eph_ls = list(range(t_eph_i, t_eph_i+(N_eph*t_eph), t_eph)) # t-axis
    
    # Let us get the time axis used in the original IGS clock array
    t_clk = 30 # Step size in IGS clock time step array, in seconds
    N_clk = len(gpsdata['t_clk']) # No. of elements in clk_30s time array
    t_clk_i = (gpsdata['t_clk'][0] - epoch).seconds # Start time of clock
    t_clk_ls = list(range(t_clk_i, t_clk_i+(N_clk*t_clk), t_clk)) # t-axis
    
    # Before we begin interpolation, check if the user wants to save plots.
    savefigs = inps['savefigs'] # User-defined option to save plots
    saverept = inps['savereport'] # User-defined option to save report
    
    # We also disable the warnings for poor polynomial fitting.
    warnings.simplefilter('ignore', np.RankWarning) # Ignore polyfit warnings
    
    print('Extracting ephemeris and clock information for GPS satellites \n')
    
    # For each PRN ID, we will begin to interpolate ephemeris and clocks
    for prn in goodsats: 
        
        ''' Main script for polynomial interpolation of GPS ephemeris '''
        
        for coord in coords:
            
            # Interpolate every 3 position points in X for a much tighter fit
            pos = gpsdata[prn]['p' + coord] # Array of coordinates
            pos2 = [] # Array of new coordinates
            vel2 = [] # Array of coordinate velocities
            
            # Perform polynomial fitting now through every 11 position points
            for p in range(0,len(pos)):
                
                if p <= len(pos) - 10:
                    coeff = np.polyfit(t_eph_ls[p:p+11],pos[p:p+11],19)
                else:
                    if len(t_eph_ls[p:]) > 2:
                        coeff = np.polyfit(t_eph_ls[p:],pos[p:],19)
                
                # Then, interpolate all values with t_usr_ls
                pt1 = t_eph_ls[p] # First timing in segmented t_eph
                pt2 = pt1 + t_eph # Second timing in segmented
                t_usr_i = [t for t in t_usr_ls if t < pt2 and t >= pt1]
                
                # Check if the time axis is an empty list:
                if t_usr_i != []:
                    pos2_np = np.polyval(coeff,t_usr_i)
                    pos2 = pos2 + pos2_np.tolist()
                
                    # Now we wish to get velocity values too.
                    t_usr_i_delta = np.array(t_usr_i) + 0.001 # Add 100ms
                    pos2_delta = np.polyval(coeff,t_usr_i_delta) # Delta-R
                    vel2 = vel2 + ((pos2_delta - pos2_np)*1000).tolist()
            
            gpsdata[prn]['p' + coord] = pos2 # m
            gpsdata[prn]['v' + coord] = vel2 # m/s
        
        ''' Main script for linear interpolation of GPS clock biases '''
                
        clkb2 = [] # New array of clock biases for this PRN satellite
        clkd2 = [] # New array of clock drifts for this PRN satellite
        
        # Parse through each time step in the clk_30s IGS file.
        for t in range(0,len(gpsdata['t_clk'])):
            
            pt1 = t_clk_ls[t] # Start of time segment
            pt2 = pt1 + t_clk # End of time segment
            
            # Linear interpolation of biases using drift as gradient.
            if t != len(gpsdata['t_clk']) - 1:

                drift = gpsdata[prn]['clkb'][t+1] - gpsdata[prn]['clkb'][t]
                drift = drift / t_clk # Rise over run -> gradient.
            
            # Extract the segment of user-defined steps within this 30s frame.
            t_usr_i = [t for t in t_usr_ls if t < pt2 and t >= pt1]
            
            # Perform the interpolation in this 30s segment.
            if len(t_usr_i) > 0:
                for k in t_usr_i:
                    clkb2.append(gpsdata[prn]['clkb'][t]+(drift*(k-pt1)) )
                    clkd2.append(drift)
        
        gpsdata[prn]['clkb'] = clkb2
        gpsdata[prn]['clkd'] = clkd2
        
        ''' Main script for plotting and saving figures of GPS satellites '''
        
        if savefigs == 'True':
            pubplt.gps_graphs(prn, t_usr_ls, gpsdata, inps)
    
    # From this on, all ephemeris and clock biases have been interpolated.
    gpsdata['t'] = [] # Update the primary time array of gpsdata dictionary
    for t in range(0,len(t_usr_ls)):
        gpsdata['t'].append( tstart + t*tstep )
    
    # Remove all the intermediate time arrays
    del gpsdata['t_usr']
    del gpsdata['t_sp3']
    del gpsdata['t_clk']
    
    # Save the report generated by parsing out GPS data, if user desires.
    if saverept == 'True':
        print('Saving graph plots and reports on GPS satellites \n')
        pubplt.gps_report(gpsdata, goodsats, inps)   
    
    # We now re-organise the gpsdata dictionary for ease of computation later.
    # gpsdata = {epoch1:{1:{'px':...,'py':...,'vz':...,'clkb':...,'clkd':...},
    #                    2:{...}, 3:{...}, ... 32:{...}}
    #            epochN:{1:{'px':...,'py':...,'vz':...,'clkb':...,'clkd':...},
    #                    2:{...}, 3:{...}, ... 32:{...}}}
    
    gpsdata_new = {} # Initialise a new structure for gpsdata
    epochs = gpsdata['t']
    for t in range(0,len(epochs)):
        epoch = epochs[t]
        gpsdata_new[ epoch ] = {}
        for p in goodsats:
            gpsdata_new[epoch][p] = {}
            gpsdata_new[epoch][p]['px'] = gpsdata[p]['px'][t]*1000 # (m)
            gpsdata_new[epoch][p]['py'] = gpsdata[p]['py'][t]*1000 # (m)
            gpsdata_new[epoch][p]['pz'] = gpsdata[p]['pz'][t]*1000 # (m)
            gpsdata_new[epoch][p]['vx'] = gpsdata[p]['vx'][t]*1000 # (m/s)
            gpsdata_new[epoch][p]['vy'] = gpsdata[p]['vy'][t]*1000 # (m/s)
            gpsdata_new[epoch][p]['vz'] = gpsdata[p]['vz'][t]*1000 # (m/s)
            gpsdata_new[epoch][p]['clkb'] = gpsdata[p]['clkb'][t] # (s)
            gpsdata_new[epoch][p]['clkd'] = gpsdata[p]['clkd'][t] # (s/s)
    
    
    print('Extraction of ephemeris and clock information completed! \n')
    return gpsdata_new, goodsats

'''' We define a function that returns the day-of-week and the GPS week. '''

def gpsweekday(t):
    
    wkday = (t.weekday() + 1) % 7 # Weekday from Python to GPST
    GPST_epoch = datetime.date(1980,1,6) # Date of GPST epoch
    user_epoch = t.date() # Get the date of the input time
    GPST_epoch_Monday = GPST_epoch - datetime.timedelta(GPST_epoch.weekday())
    user_epoch_Monday = user_epoch - datetime.timedelta(user_epoch.weekday())
    wwww = int(((user_epoch_Monday-GPST_epoch_Monday).days/7)-1) # GPS week
    
    return str(wkday), str(wwww)