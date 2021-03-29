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
## Extraction of position information of GPS satellites across time.         ##
## Derivation of velocity of GPS satellites using first order derivative.    ##
##                                                                           ##
## INPUTS:                                                                   ##
##                                                                           ##
## Final Product EPH file of GPS ephemeris, and CLK file of GPS.             ##
## (Auto-download from COD if online, and if file is not present).           ##
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
## AUTHOR MODIFIED: 12-01-2021, by Samuel Y.W. Low                           ##
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





# DELETE AFTER
import matplotlib.pyplot as plt






# IMPORT LOCAL LIBRARIES
from codes import pubplt


''' Now, this is the main routine that parses ephemeris and clock data '''

def gpsxtr(inps, tstart, tstop, tstep):
    
    warnings.simplefilter('ignore', np.RankWarning) # Ignore polyfit warnings
    
    ###########################################################################
    ###########################################################################
    
    ''' First, we initialise GPS ephemeris and clock download data. '''
    
    codurl   = 'ftp://ftp.aiub.unibe.ch/CODE/' # UBern COD URL
    
    # We will download them into our directories, below, later on.
    cwd = inps['cwd'] # Get current main working directory
    iwd = cwd + '\\input\\' # Get directory for ephemeris / clock files
    os.chdir(cwd) # Ensure the user is running in the main directory.
    
    # Then, we must retrieve all number of days of GPS and CLK data needed    
    days     = (tstop.date() - tstart.date()).days + 1 # Number of days
    filelist = [] # Stores all the required GPS / CLK ephemeris files
	
	# This is the command line call to use gzip from codes folder:
    gzip_call = '\\utils\\gzip\\gzip.exe -d '
    
    # Extract the current year.
    year = str(tstart.year)
    
    # Note, a bug will occur if processing is done across a new year... but
    # the chance of that happening is low, so I didn't fix it for now.
    
    ###########################################################################
    ###########################################################################
    
    ''' Now, we download GPS ephemeris and clock biases from IGS. '''
    
    # Now, check for desired clock files. If non-existent, download them.
    for d in range(-1,days+1):
        
        wd, wwww = gpsweekday( tstart + datetime.timedelta( days = d ) )
        name = 'COD' + wwww + wd # File name, without file extension
        clkurl = codurl + year + '/' + name + '.CLK.Z' # URL of clock file
        filelist.append(name) # Add this into the list of files for parsing
        
        # Check for CLK file, and download + unzip it if needed.
        if d in range(0,days) and os.path.exists(iwd+name+'.CLK') != True:
            
            print('CLK file for ' + name + ' not found! Downloading now at...')
            print(clkurl)
            urllib.request.urlretrieve(clkurl, name + '.CLK.Z')
            print('Completed downloading the clock file! Now unzipping...')
            subprocess.call(cwd + gzip_call + name + '.CLK.Z')
            print('Files unzipped, moving them into the inputs folder.')
            shutil.move(cwd + '\\' + name + '.CLK', iwd+name + '.CLK')
            print('Unzipping completed! \n')
            
        else:
            if d in range(0,days):
                print('CLK file for '+name+' found! Proceeding to process! \n')
    
    ###########################################################################
    ###########################################################################
    
    ''' In this segment, we extract GPS clock information from COD. '''
    
    # Now, we initialise the clock dictionary holding biases and drifts.
    
    gpsclks = {} # Initialise the dictionary that holds clock information.
    first_flag = True # Flag to mark the first reading of clock time
    
    for file in filelist[1:-1]:
        
        ''' Start reading each downloaded GPS clock (30s) file from IGS  '''
        
        tc = datetime.timedelta(seconds=30) # Step size of clock file.
        f_clk = open(iwd + file + '.CLK', 'r') # Open up the SP3 file
        
        # Read the clock file.
        for line in f_clk:
                
            # Now, let's save those clock biases and drifts.
            if 'AS G' in line:
                
                line = line.split()
                
                yyyy, mm, dd = int(line[2]), int(line[3]), int(line[4])
                hh,   mn, ss = int(line[5]), int(line[6]), int(float(line[7]))
                
                timenow = datetime.datetime(yyyy, mm, dd, hh, mn, ss)
                
                # Is the current clock readings in the user-defined time axis?
                if timenow <= tstop and timenow >= (tstart - tc):
                    
                    # Record the timing of the first epoch.
                    if first_flag == True:
                        first_time = timenow
                        first_flag = False
                    
                    # Then, record the clock biases for each SV.
                    SV      = int(line[1][1:])
                    clkbias = float(line[9])
                    
                    # Check if this SV has been recorded before.
                    if SV not in gpsclks:
                        gpsclks[SV] = {}
                    
                    # Check if this epoch has been recorded before.
                    if timenow not in gpsclks[SV]:
                        gpsclks[SV][timenow] = clkbias
                
                # Record the final time too.                    
                if timenow + tc > tstop and first_flag == False:
                    final_time = timenow
        
        # Close the current CLK file, open the next one if any.
        f_clk.close()
    
    ###########################################################################
    ###########################################################################
    
    ''' In this segment, we interpolate missing clock biases. '''
    
    # Now let's reconstruct the time axis across clock biases.
    # We create two time axes, one using datetime objects, and the other using
    # seconds in integers. This may seem duplicated, but we need the integer
    # array in order to do interpolation, whereas the datetime array is used
    # for keeping track of which missing values are found in the time axis.
    # It is not optimal code, but preference is given to the readability of it.
    
    goodsats = [] # We initialise an array to hold satellites with clock info.
    clktime_d = first_time # Initialise the first time reading (datetime)
    clkaxis_d = [] # Axis for time (datetime) across clock bias values.
    clktime_s = 0 # Initialise the first time reading (seconds)
    clkaxis_s = [] # Axis for time (seconds) across clock bias values.
    
    while clktime_d <= final_time:
        clkaxis_d.append(clktime_d) # Time axis appending of datetime objects.
        clkaxis_s.append(clktime_s) # Time axis appending of integer objects.
        clktime_d = clktime_d + tc # Add time-delta of 30 seconds.
        clktime_s = clktime_s + 30 # Add integer 30 seconds.
    
    # Let's first report which GPS satellite clock biases are completely gone.
    for SV in range(1,33):
        if SV not in gpsclks:
            print('GPS SV ' + str(SV) + ' clock biases completely missing.')
            print('SV ' + str(SV) + ' will not be used in the solver. \n')
    
    # We begin parsing through the GPS clocks dictionary to check for missing
    # clock values, and to interpolate them if possible.
    for SV in gpsclks:
        if SV not in goodsats:
            goodsats.append(SV)
            
        # Read the clock bias values for each SV as an array.
        clkbiases = list(gpsclks[SV].values()) # Clock bias values
        clkaxis_di = list(gpsclks[SV].keys()) # Recorded datetime objects
        
        # We need a time axis in floats, not datetimes, so we can do the
        # interpolation. Since the full time axis 'clkaxis_d' should be sorted
        # thanks to IGS, we can simplify our search of the sorted list by
        # matching the indices, instead of using brute force or binary search.
        clkaxis_si = [30.0 * clkaxis_d.index(t) for t in clkaxis_di]
        
        # Any SVs that were not recorded, will be not be entered in 'goodsats'.
        # However, there could be missing clock values scattered across time,
        # for a particular SV value. We should interpolate those values here.
        
        if len(clkaxis_si) != len(clkaxis_s):
            
            # Perform a least squares regression best fit line.
            bestfit = np.polyfit(clkaxis_si, clkbiases, 1)
            
            # Then, fill in all missing clock values.
            # Time axis 'clkaxis_d' contains all ideal epochs, without missing
            # values. We will compare its elements to the elements of the
            # epochs in the present 'timeaxisd', and if there are any missing
            # epochs, then we perform extrapolation or interpolation.
            
            for t in range(0,len(clkaxis_d)):
                
                # Now, we check which epochs in the time axis is missing.
                if clkaxis_d[t] not in clkaxis_di:
                    
                    # Extrapolate the clock bias.
                    clkbias_extr = float(np.polyval(bestfit, clkaxis_s[t]))
                    
                    # Update this interpolated value in the 'gpsclks'
                    gpsclks[SV][clkaxis_d[t]] = clkbias_extr
    
    # Sort the list of good satellites found.
    goodsats.sort() # Sort in ascending order
    
    ###########################################################################
    ###########################################################################
    
    ''' In this segment, we initialise the structure of our final output. '''
   
    # We initialise the structure of the final output dictionary.
    # gpsdata = {epoch1: { SV1: { 'px': ..., 'vz': ..., 'clkb':... } ...} ...}
    # Unlike gpsephm and gpsclks, the final gpsdata is indexed across time.
    gpsdata = {} # Final output object.
    ti = tstart # Initialise the time object.
    
    # Begin initialising the 'gpsdata' dictionary object;
    # (the final output of this program).
    while ti <= tstop:
        gpsdata[ti] = {}
        for p in goodsats:
            gpsdata[ti][p] = {}
            gpsdata[ti][p]['px']   = 0.0 # Position X of LEO at t = ti
            gpsdata[ti][p]['py']   = 0.0 # Position Y of LEO at t = ti
            gpsdata[ti][p]['pz']   = 0.0 # Position Z of LEO at t = ti
            gpsdata[ti][p]['vx']   = 0.0 # Velocity X of LEO at t = ti
            gpsdata[ti][p]['vy']   = 0.0 # Velocity Y of LEO at t = ti
            gpsdata[ti][p]['vz']   = 0.0 # Velocity Z of LEO at t = ti
            gpsdata[ti][p]['clkb'] = 0.0 # Clock bias of LEO at t = ti
            gpsdata[ti][p]['clkd'] = 0.0 # Clock drift of LEO at t = ti
        ti = ti + tstep # Update the time step.
    
    # From 'gpsdata' we can get the user-specified time axis in date-time.
    tstep_ss  = tstep.seconds
    t_usr_dt = sorted(list(gpsdata.keys()))
    t_usr_ss = np.array(sorted([tstep_ss * t for t in range(0,len(t_usr_dt))]))
    
    ###########################################################################
    ###########################################################################
    
    ''' In this segment, we extract GPS precise ephemerides from COD. '''
    
    # Now, check for desired ephemeris files. If non-existent, download them.
    for d in range(-1,days+1):
        
        wd, wwww = gpsweekday( tstart + datetime.timedelta( days = d ) )
        name = 'COD' + wwww + wd # File name, without file extension
        ephurl = codurl + year + '/' + name + '.EPH.Z' # URL of ephemeris file
    
        # Check for SP3 ephemeris file, and download + unzip it if needed.
        if os.path.exists(iwd+name+'.EPH') != True:
            
            print('EPH file for '+ name +' not found! Attempt download now...')
            urllib.request.urlretrieve(ephurl, name + '.EPH.Z')
            print('Completed downloading the ephemeris file! Now unzipping...')
            subprocess.call(cwd + gzip_call + name + '.EPH.Z')
            print('Files unzipped, moving them into the inputs folder.')
            shutil.move(cwd + '\\' + name + '.EPH', iwd + name + '.EPH')
            print('Unzipping completed! \n')
        
        else:
            
            print('EPH file for ' + name + ' found! Proceeding to process! \n')

    # We download the GPS precise ephemerides one day before and after,
    # to prevent the occurrence of Runge's phenomenon by adding buffer in the
    # extrapolation process.
    
    # This will be where ephemeris data will be parsed into.
    # gpsphm = {SV: {epoch1: {'px':xxx, 'py':yyy ... 'vy':vyy, 'vz',vzz}}}
    gpsephm = {} # Ephemeris dictionary.
    for SV in goodsats:
        gpsephm[SV] = {}
    
    # All epochs after 'gps_tstart' will be recorded from IGS ephemeris file.
    gps_tstart = tstart - datetime.timedelta(seconds = 7200) # Minus 02:00:00
    
    # All epochs will be recorded until 'gps_tstop', inclusive.
    gps_tstop  = tstop + datetime.timedelta(seconds = 7200) # Add 02:00:00
    
    # Time step in the IGS precise ephemeris file.
    gps_tstep  = datetime.timedelta(seconds = 900)
    
    # Flag to mark the first reading of ephemeris time
    first_flag = True
    
    # Now we parse through all the downloaded IGS files.
    for file in filelist:
        
        ''' Start reading each downloaded GPS ephemeris file from IGS  '''
        
        f_gps = open(iwd + file + '.EPH', 'r') # Open up the SP3 file
        gps_record = False # Trigger for recording GPS data
        
        for line in f_gps:
            
            # Split up the strings
            line = line.split()
            
            # Skip blank lines
            if len(line) <= 1:
                continue
            
            # Check if the time reading is accurate...
            if line[0] == '*':
                
                # Read the time reading now, save it as a datetime object.
                timenow = datetime.datetime(int(line[1]),
                                            int(line[2]),
                                            int(line[3]),
                                            int(line[4]),
                                            int(line[5]),
                                            int(float(line[6])))
                    
                # ... only if the epoch falls within our desired range.
                if timenow <= gps_tstop and timenow > gps_tstart - gps_tstep:
                    gps_record = True # Trigger the recording on.
                    
                    # Record the first ever epoch.
                    if first_flag == True:
                        first_time = timenow
                        first_flag = False
                        
                else:
                    gps_record = False # Trigger the recording off.
                        
            # Now, check for XYZ coordinates of each SV ID
            if 'PG' in line[0] and gps_record == True:
                
                # Get the data we need below.
                SV  = int(line[0][2:4]) # Which GPS satellite is this?
                x   = float(line[1]) # X-coordinate position
                y   = float(line[2]) # Y-coordinate position
                z   = float(line[3]) # Z-coordinate position
                
                # Check if this SV has clean clock values.
                if SV in goodsats:
                
                    # Check if this is the first epoch.
                    if timenow not in gpsephm[SV]:
                        gpsephm[SV][timenow] = {}
                    
                    # Assign the data into the gpsephm dictionary.
                    gpsephm[SV][timenow] = {} # Initialise...
                    gpsephm[SV][timenow]['px'] = x # Position X (ECEF)
                    gpsephm[SV][timenow]['py'] = y # Position Y (ECEF)
                    gpsephm[SV][timenow]['pz'] = z # Position Z (ECEF)
                    
                    # Initialise the velocity values too, for estimation.
                    gpsephm[SV][timenow]['vx'] = 0.0 # Velocity X (ECEF)
                    gpsephm[SV][timenow]['vy'] = 0.0 # Velocity Y (ECEF)
                    gpsephm[SV][timenow]['vz'] = 0.0 # Velocity Z (ECEF)
                
        f_gps.close()
    
    ###########################################################################
    ###########################################################################
    
    ''' In this segment, we interpolate the GPS precise ephemeris. '''
    
    print('We now interpolate the GPS precise ephemeris. \n')
    
    # We need to add about two hours of buffer before and after the validity
    # period for interpolation. See research paper "Polynomial interpolation
    # of GPS satellite coordinates" by Milan Horemuz (Feb 2006).
    # Thus, we use a sliding window interpolation, with fit interval of 4h,
    # and a validity window in the middle of 15 minutes only, leaving the 
    # rest of the 4h as interpolation buffer to prevent Runge's phenomenon.
    
    # Our first task is to generate the IGS ephemeris time axis.
    gpsephm_epoch_dt = [] # List of datetime objects.
    gpsephm_epoch_ss = [] # List of objects in integer seconds.
    gpsephm_dt = first_time # First epoch, as a datetime object.
    gpsephm_ss = 0 # First epoch but in integer seconds.
    
    # Fill the list.
    while gpsephm_dt <= gps_tstop:
        gpsephm_epoch_dt.append(gpsephm_dt)
        gpsephm_epoch_ss.append(gpsephm_ss)
        gpsephm_dt = gpsephm_dt + gps_tstep
        gpsephm_ss = gpsephm_ss + gps_tstep.seconds
    
    # We also initialise the starting time in GPS ephemeris interpolation.
    # This block of code is meant to create a time axis in integer seconds.
    t_offset_eph = (tstart-first_time).days*86400 + (tstart-first_time).seconds
    
    # We offset the interpolant time axis with the time axis offset.
    t_usr_eph = sorted(t_usr_ss + t_offset_eph)
    
    # Also, it would be wise to check if the IGS precise ephemerides had any
    # entries with missing GPS epochs. We simply check the lengths of arrays.
    for p in gpsephm:
        if len(gpsephm_epoch_dt) != len(list(gpsephm[p].keys())):
            print('Warning! IGS ephemerides missing epochs for SV ' + str(p))
            print('Error will occur during interpolation process!')
    
    # The interpolation now starts when iterable 't' > gpsephm_epoch_dt[0],
    # the first element of the array, and it will end when 't' exceeds 'tstop',
    # the user defined ending of the time axis.
    
    # We loop through each fit interval first, interpolate it, and then move
    # our sliding window interpolator every 'window_len' of 2 hours.
    
    coords = ['x', 'y', 'z'] # Position coordinate keys in dictionary.
    windex = 0 # Window index keeps track of sliding window interpolant.
    
    # Sanity check, GPS orbits must be at least 4 hours long for interpolation.
    if len(gpsephm_epoch_dt) <= 17:
        print('Warning! GPS satellite interpolation is not possible!')
        print('Duration of time for orbit interpolation is too short!')
        print('Please use at least 2 hours of duration in scenario!')
        print('Returning an error... \n')
        return {}, []
    
    while gpsephm_epoch_dt[windex+16] != gpsephm_epoch_dt[-1]:
        
        # First let us get the fit and validity window (as datetime objects).
        fit_dt = gpsephm_epoch_dt[windex : windex + 17]
        val_dt = gpsephm_epoch_dt[windex + 7 : windex + 11]
        
        # Then we get the equivalent windows in terms of integer seconds.
        fit_ss = gpsephm_epoch_ss[windex : windex + 17]
        val_ss = gpsephm_epoch_ss[windex + 7 : windex + 11]
            
        # Now we generate the interpolated time axis.
        intp_dt = [t for t in t_usr_dt if t <= val_dt[-1] and t >= val_dt[0]]
        intp_ss = [t for t in t_usr_eph if t <= val_ss[-1] and t >= val_ss[0]]
            
        # Numpy-rize the arrays, and also create a + 1ms time axis.
        # 'intp_ss_delta' will be used for velocity estimation (1st order)
        intp_ss = np.array(intp_ss)
        intp_ss_delta = intp_ss + 0.00001 # Add one ms
        
        # In this sliding window filter, intepolate across SVs.
        for SV in goodsats:
            
            # For each axes in the coordinate frame...
            for c in coords:
                
                # Get corresponding fit and validity intervals for ephemerides.
                gpsephm_fit = [gpsephm[SV][t]['p'+c] for t in fit_dt]
                
                # Perform 16th Order Polyfit Interpolation.
                poly = np.polyfit( fit_ss, gpsephm_fit, 16 )
                
                # Now, we can perform the interpolation for positions.
                gpsephm_intp = np.polyval(poly, intp_ss)
                
                # We also exploit the interpolation for velocity estimation.
                gpsephm_intp_delta = np.polyval(poly, intp_ss_delta)
                gpsephm_intv = (gpsephm_intp_delta - gpsephm_intp) * 100000
                
                # Add the interpolated positions and velocities into 'gpsdata'.
                for k in range(0,len(gpsephm_intp)):
                    
                    gpsdata[intp_dt[k]][SV]['p'+c] = gpsephm_intp[k] # Pos
                    gpsdata[intp_dt[k]][SV]['v'+c] = gpsephm_intv[k] # Vel
        
        # Check if we have reached the end of the while loop.
        if gpsephm_epoch_dt[windex+16] == gpsephm_epoch_dt[-1]:
            break # End the while loop if we are.
        else:
            windex += 1 # Update the sliding window filter by 2 hours
    
    # # For debugging
    # xt = []
    # xx= []
    # yp = []
    # ypi = []
    # count = t_offset_eph
    # for t in gpsdata.keys():
    #     xt.append(count)
    #     ypi.append(gpsdata[t][30]['px'])
    #     count += tstep.seconds
    # count = 0
    # for t1 in gpsephm_epoch_dt:
    #     xx.append(count)
    #     yp.append(gpsephm[30][t1]['px'])
    #     count += gps_tstep.seconds
    # plt.figure(1)
    # plt.plot(xt,ypi)
    # plt.scatter(xx,yp)
    # plt.show()
    
    ###########################################################################
    ###########################################################################
    
    ''' Finally, we interpolate clock biases to the user's time axis. '''
    
    print('Interpolation of GPS precise ephemerides done!')
    print('Now interpolating GPS clock biases. \n')
    
    # We also initialise the starting time in GPS clock biases interpolation.
    # This block of code is meant to create a time axis in integer seconds.
    
    # From 'gpsdata' we can get the user-specified time axis in date-time.
    tstep_ss  = tstep.seconds
    t_usr_dt = sorted(list(gpsdata.keys()))
    t_usr_ss = np.array(sorted([tstep_ss * t for t in range(0,len(t_usr_dt))]))
    
    # We then adjust the time axis of the user's to interpolate for 'clkb'.
    first_time = sorted(list(gpsclks[goodsats[0]].keys()))[0]
    t_offset_clk = (tstart-first_time).days*86400 + (tstart-first_time).seconds
    t_usr_clk = t_usr_ss + t_offset_clk
    
    # Now, start going through each SV for clock bias interpolation.    
    for SV in goodsats:
        
        # Read the clock bias values for each SV as an array.
        clkaxis_di = list(gpsclks[SV].keys()) # Recorded datetime objects.
        clkaxis_si = [30 * t for t in range(0,len(clkaxis_di))]
        clkbiases = list(gpsclks[SV].values()) # Clock bias values.
        clkbiases = [x for _,x in sorted(zip(clkaxis_di, clkbiases))]
        clk_bound = 0 # Lower bound index for clock bias.
        lower_bound = clkaxis_si[ clk_bound ] # Nominal lower bound time (sec).
        
        for k in range(0,len(t_usr_clk)):
            tss = t_usr_clk[k] # Epoch in seconds, on user-defined time axis.
            tdt = t_usr_dt[k] # Epoch in datetime, on user-defined time axis.
            
            while clkaxis_si[ clk_bound ] + 30 <= tss:
                clk_bound += 1 # Check the next time-bound in gpsclks.
                lower_bound = clkaxis_si[ clk_bound ] # Update the bound.
            
            # Now we perform the interpolation. First, we need the clock drift.
            if clk_bound + 1 < len(clkbiases):
                
                # Get the clock drift as the gradient between two clock biases.
                clkdrift = (clkbiases[clk_bound+1] - clkbiases[clk_bound]) / 30
            
            # Get the desired first order Delta-T.
            delta_tt = tss - lower_bound
            
            # Interpolate using (y + delta_y) = x + (drift * delta_x)
            clkbias_interp = clkbiases[ clk_bound ] + clkdrift*( delta_tt )
            
            # Append interpolated clock biases into 'gpsdata' output.
            gpsdata[tdt][SV]['clkb'] = clkbias_interp
            gpsdata[tdt][SV]['clkd'] = clkdrift
    
    ###########################################################################
    ###########################################################################
    
    ''' In a final step, we plot the GPS ephemeris and clock biases. '''
    
    # First, check if the user wishes to plot GPS ephemeris and clock biases.
    savefigs = inps['savefigs'] # User-defined option to save plots
    saverept = inps['savereport'] # User-defined option to save report
    
    # If so, then continue to save the final output report on GPS ephemeris.
    if saverept == 'True':
        print('Saving output report on interpolated GPS ephemerides. \n')
        pubplt.gps_report(gpsdata, goodsats, inps)
    
    # ... as well as ephemeris and clock plots.
    if savefigs == 'True':
        print('Saving plots on GPS position, velocity and clock biases. \n')
        for SV in goodsats:
            pubplt.gps_graphs(SV, t_usr_dt, t_usr_ss, gpsdata, inps)
    
    ###########################################################################
    ###########################################################################
    
    return gpsdata, goodsats

'''' We define a function that returns the day-of-week and the GPS week. '''

def gpsweekday(t):
    
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