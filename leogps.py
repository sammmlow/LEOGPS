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
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 03-12-2019.                                              ##
##                                                                           ##
###############################################################################
###############################################################################
'''

# Import local libraries
from codes import inpxtr
from codes import rnpath
from codes import timing
from codes import gpsxtr
from codes import rinxtr
from codes import posvel
from codes import ambest
from codes import pubplt

def main():
    
    # Let's start importing user-defined parameters!
    inps = inpxtr.inpxtr()
    
    # Get the RINEX file paths.
    rinex1f, rinex2f = rnpath.rnpath( inps )
    
    # Now, let's import our scenario timing parameters.
    timecheck = timing.tcheck( rinex1f, rinex2f, inps )
    
    # Check if the time axis has any conflicts?
    if timecheck == False:
        return None
    
    # Otherwise, proceed to get timings!
    tstart, tstop, tstep, rnxstep = timecheck
    
    # Grab XYZ and clock data from good GPS satellites.
    gpsdata, goodsats = gpsxtr.gpsxtr( inps, tstart, tstop, tstep )
    
    # Extract the RINEX observables.
    rinex1 = rinxtr.rinxtr(rinex1f, inps, goodsats, tstart, tstop, rnxstep)
    rinex2 = rinxtr.rinxtr(rinex2f, inps, goodsats, tstart, tstop, rnxstep)

    # Check if the RINEX files are corrupted?
    if rinex1 == False or rinex2 == False:
        return None
    
    # # Initialise and setup the time axis for processing.
    ti, time, results = tstart, [], {}
    while ti < tstop:
        time.append(ti)
        ti = ti + tstep
    
    # Now, begin getting PVT and baseline lengths of both LEOs.
    for t in time:
        
        # What are the GPS observables and RINEX observables at t?
        gps = gpsdata[t]
        rx1 = rinex1[t]
        rx2 = rinex2[t]
        
        # Extract PVT and DOP values for LEO satellite 1, at time = t.
        pos1, vel1, dop1 = posvel.posvel(t, goodsats, gps, rx1, inps)
        
        # Extract PVT and DOP values for LEO satellite 2, at time = t.
        pos2, vel2, dop2 = posvel.posvel(t, goodsats, gps, rx2, inps)
        
        # Perform double-differencing to get baseline length.
        baseline = ambest.ambest(t, gps, rx1, rx2, pos1, pos2, inps)
        
        # Log the results into a dictionary.
        results[t] = [pos1, vel1, dop1, pos2, vel2, dop2, baseline]
    
    # Publish the results of LEOGPS.
    pubplt.leo_results( results, inps )
    
    return None

main()
