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
##    Extraction of input information from 'config.txt', by searching        ##
##    for the config.txt file, and then parses the user-defined inputs.      ##
##    Outputs a dictionary holding user-defined inputs from 'config.txt',    ##
##    which can then be conveniently called in other python scripts as       ##
##    function 'inpxtr.inpxtr()'                                             ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 07-Jun-2021.                                             ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##    Documentation: https://leogps.readthedocs.io/en/latest/                ##
##                                                                           ##
###############################################################################
###############################################################################

# Import global libraries
import datetime
from os.path import dirname, abspath, join

# We need a function that converts dates and times into GPS formats.

def inptim(t):
    '''Extract time parameters from **config.txt**, and translates them into 
    five timing-related parameters to be used for processing later on. Called
    internally only in `inpxtr.inpxtr()`.
    
    Parameters
    ----------
    t : datetime.datetime
        An epoch datetime object.
    
    Returns
    -------
    yyyy : int
        4-digit Gregorian year
    yy : str
        2-digit Gregorian year (for RINEX file name)
    doy : str
        3-digit day-of-year (for RINEX file name)
    wkday : int
        1-digit day-of-week (Sunday = 0, Saturday = 6)
    wwww : int
        GPS week number (generally 4-digits)
    
    '''
    # t is a datetime.datetime object
    
    yyyy = t.year # The four digit representation of Gregorian year
    yy = str(yyyy)[2:] # The two digit representation of Gregorian year
    doy = t - datetime.datetime(yyyy,1,1) # Subtract today from first day
    doy = str( int(doy.days) + 1 ) # Which day of the year is it (0 to 365)?
    wkday = (t.weekday() + 1) % 7 # Weekday from Python to GPST (Sunday = 0)
        
    # Now, we get the GPS week for GPST
    GPST_epoch        = datetime.date(1980,1,6)
    user_epoch        = t.date()
    GPST_epoch_Monday = GPST_epoch - datetime.timedelta(GPST_epoch.weekday())
    user_epoch_Monday = user_epoch - datetime.timedelta(user_epoch.weekday())
    
    # Get the GPS week number
    wwww = int( ( ( user_epoch_Monday - GPST_epoch_Monday).days / 7 ) - 1 )
    
    return yyyy, yy, doy, wkday, wwww
    

# We have to populate our inputs in the Python script now.

def inpxtr():
    '''Extract all parameters from **config.txt**. No input arguments.
    
    Returns
    -------
    inputdict : dict
        A dictionary of key-value pairs comprising of key-values found in the
        **config.txt** file. For example, the time step of the scenario would
        be the key-value pair `{ 'timestep' : 30 }`
        
    '''
    
    # First, some housekeeping to get the required files you need.
    
    cwd = dirname(dirname(abspath(__file__))) # Current working directory (str)
    iwd = join(cwd, 'config', 'config.txt') # Inputs files directory (str)
    
    # Some characteristic configuration parameters are accounted for here.
    
    integers = [ 'freq','timestep','hatchlength','cycsliplen' ]
    
    inputfile = open(iwd,'r')
    inputdict = {} # Create a dictionary to store all the input 
    
    for line in inputfile:
        
        # Check for input entry with an 'I'
        
        if line[0] == 'I':
            
            # First string splitting and formatting
            
            line_inp = line[3:].split()
            
            # Check for values that should be expressed in integers
            
            if line_inp[0] in integers:
                
                inputdict[ line_inp[0] ] = int(line_inp[1])
            
            # Check for values that should be expressed in date and time
            
            elif line_inp[0] == 'dtstart':
                
                dtlist = line_inp[1].split('-')
                
                dtstart = datetime.datetime(int(dtlist[0]),
                                            int(dtlist[1]),
                                            int(dtlist[2]),
                                            int(dtlist[3]),
                                            int(dtlist[4]),
                                            int(dtlist[5]))
                
                # Now, we input all starting epoch time parameters!
                
                inputdict[ 'dtstart' ]       = dtstart # Datetime object
                inputdict[ 'dtstart_yyyy' ]  = inptim(dtstart)[0] # int
                inputdict[ 'dtstart_yy' ]    = inptim(dtstart)[1] # str
                inputdict[ 'dtstart_doy' ]   = inptim(dtstart)[2] # str
                inputdict[ 'dtstart_wkday' ] = inptim(dtstart)[3] # int
                inputdict[ 'dtstart_wwww' ]  = inptim(dtstart)[4] # int

                
            # Check for values that should be expressed in date and time
            
            elif line_inp[0] == 'dtstop':
                
                dtlist = line_inp[1].split('-')
                
                dtstop = datetime.datetime(int(dtlist[0]),
                                           int(dtlist[1]),
                                           int(dtlist[2]),
                                           int(dtlist[3]),
                                           int(dtlist[4]),
                                           int(dtlist[5]))
                
                # Now, we input all ending epoch time parameters!
                
                inputdict[ 'dtstop' ]       = dtstop # Datetime object
                inputdict[ 'dtstop_yyyy' ]  = inptim(dtstop)[0] # int
                inputdict[ 'dtstop_yy' ]    = inptim(dtstop)[1] # str
                inputdict[ 'dtstop_doy' ]   = inptim(dtstop)[2] # str
                inputdict[ 'dtstop_wkday' ] = inptim(dtstop)[3] # int
                inputdict[ 'dtstop_wwww' ]  = inptim(dtstop)[4] # int    
            
            # Else, just append input parameters into input dictionary
            else:
                
                inputdict[ line_inp[0] ] = line_inp[1]
    
    # Close the file when done
    inputfile.close()
    
    # Add the current working directory into inpxtr for future reference
    inputdict['cwd'] = cwd
    
    # Return the inputs for future use
    return inputdict