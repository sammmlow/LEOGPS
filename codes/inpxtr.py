#!/usr/bin/env python3
'''
###############################################################################
###############################################################################
##                                                                           ##
##     _    ___  ___   ___ ___ ___                                           ##
##    | |  | __ /   \ / __| _ | __|                                          ##
##    | |__| __  ( ) | (_ |  _|__ \                                          ##
##    |____|___ \___/ \___|_| \___/                                          ##
##                                    v 1.1 (Stable)                         ##
##                                                                           ##
## FILE DESCRIPTION:                                                         ##
##                                                                           ##
## Extraction of input information from a 'config.txt' file.                 ##
##                                                                           ##
## INPUTS:                                                                   ##
##                                                                           ##
## Searches for the config.txt file, and then parses the user-defined inputs ##
##                                                                           ##
## OUTPUT:                                                                   ##
##                                                                           ##
## A variable holding all user-defined inputs from 'config.txt', which can   ##
## then be conveniently called in other python scripts as function 'inpxtr()'##
##                                                                           ##
## AUTHOR MODIFIED: 27-11-2019, by Samuel Y.W. Low                           ##
##                                                                           ##
###############################################################################
###############################################################################
'''

import datetime
from os.path import dirname, abspath

''' We need a function that converts dates and times into GPS formats '''

def inptim(t):
    
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
    

''' We have to populate our inputs in the Python script now '''

def inpxtr():
    
    ''' First, some housekeeping to get the required files you need '''

    cwd = dirname(dirname(abspath(__file__))) # Current working directory
    iwd = cwd + '\config\config.txt' # Inputs files
    
    ''' Some characteristic configuration parameters are accounted for here '''
    
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
                
            # Check for the array of bad GPS satellites
            elif line_inp[0] == 'badsats':
                
                badsats = []
                badsatlist = line_inp[1].split(',')
                
                for elem in badsatlist:
                    badsats.append(int(elem))
                    
                inputdict['badsats'] = badsats
            
            # Else, just append input parameters into input dictionary
            else:
                
                inputdict[ line_inp[0] ] = line_inp[1]
    
    # Close the file when done
    inputfile.close()
    
    # Add the current working directory into inpxtr for future reference
    inputdict['cwd'] = cwd
    
    # Return the inputs for future use
    return inputdict