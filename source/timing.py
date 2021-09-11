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
##    This module comprises two functions: tcheck() and get_startstop().     ##
##                                                                           ##
##    tcheck() takes in the paths (pathlib.Path) of 02x RINEX v2 files,      ##
##    as well as the user-defined dictionary of inputs from config.txt.      ##
##    If both RINEX files and the user-defined timing parameters do not      ##
##    have an overlapping duration (i.e. each RINEX file corresponds to      ##
##    a different date or time), or if they have incoherent step sizes       ##
##    then this function will output False. Else, this function will out-    ##
##    put the scenario start time, scenario final time, the scenario step    ##
##    size, as well as the RINEX file step size.                             ##
##                                                                           ##
##    get_startstop() takes in a list, with elements of the list being       ##
##    each line of the RINEX file, as strings, and returns three objects;    ##
##    the observed start and stop times as datetime.datetime objects, and    ##
##    the RINEX file time step as a datetime.timedelta object.               ##
##                                                                           ##
##    The rest of LEOGPS will base its time on the output of tcheck().       ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 05-Jun-2021.                                             ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##    Documentation: https://leogps.readthedocs.io/en/latest/                ##
##                                                                           ##
###############################################################################
###############################################################################

# Import global libraries
import datetime

def tcheck(rnx1file, rnx2file, inps):
    
    ''' Reads both RINEX files and output the desired start and stop times.
    Sets the start, stop, and time steps for the entire LEOGPS scenario.
    
    Parameters
    ----------
    rnx1file : pathlib.Path
        Path of RINEX file of LEO-A
    rnx2file : pathlib.Path
        Path of RINEX file of LEO-B
    inps : dict
        A dictionary of inputs created by `inpxtr.inpxtr()`
    
    Returns
    -------
    tstart : datetime.datetime
        Scenario start time for processing
    tstop : datetime.datetime
        Scenario stop time for processing
    tstep : datetime.timedelta
        Scenario time step used in processing
    rnx1step : datetime.timedelta
        Observed time step in the RINEX file
    
    '''
    
    # Let's first import all user-defined inputs
    name1    = inps['name1'] # 4-letter ID of spacecraft A
    name2    = inps['name2'] # 4-letter ID of spacecraft B
    dtstart  = inps['dtstart'] # The user-defined start time
    dtstop   = inps['dtstop'] # The user-defined stop time
    userstep = datetime.timedelta(seconds = inps['timestep']) # Step size
    
    # INSERT RNXFILE EXTRACTION HERE
    
    print('The user-defined start time of simulations is: '+str(dtstart))
    print('The user-defined end time of simulations is: '+str(dtstop))
    print('\n')
    
    # Grab the start-stop times for RINEX file 1
    with open(rnx1file) as file1:
        
        rnx1lines = file1.readlines()
        rnx1start, rnx1stop, rnx1step = get_startstop(rnx1lines)

    print('The start of RINEX observables in '+name1+' is: '+str(rnx1start))
    print('The end of RINEX observables in '+name1+' is: '+str(rnx1stop))
    print('The step of RINEX observables in '+name1+' is: ' + str(rnx1step))
    print('\n')
    
    # Grab the start-stop times for RINEX file 2
    with open(rnx2file) as file2:
        
        rnx2lines = file2.readlines()
        rnx2start, rnx2stop, rnx2step = get_startstop(rnx2lines)
    
    print('The start of RINEX observables in '+name2+' is: ' + str(rnx2start))
    print('The end of RINEX observables in '+name2+' is: ' + str(rnx2stop))
    print('The step of RINEX observables in '+name2+' is: ' + str(rnx2step))
    print('\n')
    
    # Check if RINEX files 1 and 2 have the same time step
    
    if rnx1step != rnx2step:
        
        print('Error! Both RINEX files do not have the same time step!')
        return False
    
    # Check if RINEX files 1 and 2 have an overlapping intersection of time
    
    if (rnx2start > rnx1stop) or (rnx1start > rnx2stop):
        
        print('Error! Both RINEX files do not have a common time overlap!')
        rnxstart = datetime.datetime(1980,1,6,0,0,0) # Set GPST epoch default
        rnxstop  = datetime.datetime(1980,1,6,0,0,0) # Set GPST epoch default
        
        return False
    
    # If they do have overlapping durations, pick the intersection
    
    else:
        
        rnxstart = max([rnx1start,rnx2start])
        rnxstop  = min([rnx1stop,rnx2stop])
        
        tstart = max([rnxstart,dtstart])
        tstop  = min([rnxstop,dtstop])
        
        if dtstart > rnxstop:
            
            print('Error!')
            print('*config.txt dtstart* is later than last observable!')
            print('LEOGPS will use the first possible RINEX observable.')
            tstart = rnxstart # Ignore config.txt
        
        if dtstop < rnxstart:
            
            print('Error!')
            print('*config.txt dtstop* is earlier than first observable!')
            print('LEOGPS will use the last possible RINEX observable.')
            tstop = rnxstop # Ignore config.txt
        
        if dtstart < rnxstart:
            
            print('Warning!')
            print('*config.txt dtstart* is earlier than first observable!')
            print('LEOGPS will use the first possible RINEX observable start!')
            print('\n')
            
        if dtstop > rnxstop:
            
            print('Warning!')
            print('*config.txt dtstop* is later than last observable!')
            print('LEOGPS will use the latest possible RINEX observable stop!')
            print('\n')
        
        tstep = max([userstep,rnx1step])
        
        if userstep < rnx1step:
            print('Warning!')
            print('*config.txt timestep* is smaller than RINEX timestep!')
            print('LEOGPS will use the time step in the RINEX observations.')
            print('\n')
            tstep = rnx1step
            
        elif userstep.seconds % rnx1step.seconds != 0:
            print('Error!')
            print('*config.txt timestep* is not a multiple of RINEX timestep!')
            print('LEOGPS will use the time step in the RINEX observations.')
            print('\n')
            tstep = rnx1step
            
    print('The scenario will use the start time ' + str(tstart))
    print('The scenario will use the stop time ' + str(tstop))
    print('\n')
            
    return tstart, tstop, tstep, rnx1step



def get_startstop(rnxlines):
    
    ''' Gets the first and last epoch after reading one RINEX file. 
    
    Parameters
    ----------
    rnxlines : list
        List of strings as lines of the RINEX file
    
    Returns
    -------
    rnxstart : datetime.datetime
        Observed initial time stamp in RINEX file
    rnxstop : datetime.datetime
        Observed final time stamp in RINEX file
    rnxstep : datetime.timedelta
        Observed time step in RINEX file
    
    '''
    
    # Check for where the end of the RINEX header is located at.
    header_reached = False
    header_linenum = 0
    while header_reached == False:
        headerline = rnxlines[header_linenum]
        if 'END OF HEADER' in headerline:
            header_reached = True
        header_linenum += 1
    
    tcount = 0
    
    # Get the starting date-time
    for data in rnxlines[header_linenum:]:
        
        data = data.replace('G',' ') # GPS satellites
        dataspl = data.split()
        
        # The entry should be a date-time entry if there are > 7 elements
        if len(dataspl) >= 8:
            
            if tcount == 0:
                                
                # Get the starting epoch of RINEX file
                rnxstart = datetime.datetime(int(dataspl[0]) + 2000,
                                             int(dataspl[1]),
                                             int(dataspl[2]),
                                             int(dataspl[3]),
                                             int(dataspl[4]),
                                             int(round(float(dataspl[5]))))
                tcount += 1
            
            elif tcount == 1:
                
                # Get the second epoch of RINEX file
                rnxstart2 = datetime.datetime(int(dataspl[0]) + 2000,
                                              int(dataspl[1]),
                                              int(dataspl[2]),
                                              int(dataspl[3]),
                                              int(dataspl[4]),
                                              int(round(float(dataspl[5]))))
                
                rnxstep = rnxstart2 - rnxstart # Get time delta
                break
    
    # Get the last date-time
    for data in reversed(rnxlines):
        
        data = data.replace('G',' ') # GPS satellites
        dataspl = data.split()
        
        # The entry should be a date-time entry if there are > 7 elements
        if len(dataspl) >= 8:
            
            rnxstop = datetime.datetime(int(dataspl[0]) + 2000,
                                         int(dataspl[1]),
                                         int(dataspl[2]),
                                         int(dataspl[3]),
                                         int(dataspl[4]),
                                         int(round(float(dataspl[5]))))
            break
        
    return rnxstart, rnxstop, rnxstep
