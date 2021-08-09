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
##    This module has a single function called run(), which is basically     ##
##    the main process workflow for LEOGPS to perform double differencing    ##
##    baseline estimation. The function in this module, leorun.run(), is     ##
##    called in the GUI (leogui.py) when "Run LEOGPS" is clicked on.         ##
##    This function requires no inputs, and returns None.                    ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 08-Aug-2021                                              ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##    Documentation: https://leogps.readthedocs.io/en/latest/                ##
##                                                                           ##
###############################################################################
###############################################################################

# Import local libraries
from source import inpxtr
from source import rnpath
from source import timing
from source import gpsxtr
from source import rinxtr
from source import posvel
from source import ambest
from source import frames
from source import pubplt

def run():
    '''Basically runs the primary workflow for LEOGPS' relative positioning.
    No input arguments needed, and returns None.
    
    '''

    # IMPORTING USER DEFINED PARAMETERS:
    # ==================================
    # -> The routine 'inpxtr.py' extracts out all user-specified parameters.
    # -> It returns a dictionary object comprising each parameter.
    
    inps = inpxtr.inpxtr()
    
    n1  = inps['name1'] # Name of LEO1
    n2  = inps['name2'] # Name of LEO2

    # FINDING THE RINEX FILE PATHS:
    # =============================
    # -> The routine 'rnpath.py' will output the path of the RINEX file.
    # -> This is based on the 4-letter spacecraft ID input in 'config.txt'.
    # -> Do not change any of the folder names in LEOGPS!
    
    rinex1f, rinex2f = rnpath.rnpath( inps )

    # FINDING ALL TIMING PARAMETERS:
    # ==============================
    # -> The routine 'timing.py', will take in three sets of timings.
    # -> First, the user-specified start-stop times in 'config.txt'.
    # -> Next, the two start-stop time lines based on both of the RINEX files.
    # -> Then, 'timing.py' will output the intersection of all three timelines.
    # -> The routine also error-checks the timesteps, based on the RINEX files.
    
    timecheck = timing.tcheck( rinex1f, rinex2f, inps )

    # CHECK FOR ANY TIMING CONFLICTS:
    # ===============================
    # -> If time conflicts arise, timing.tcheck(...) will return a False.
    # -> One possible conflict is if the user-defined start-stop times,
    # -> are completely out of the time line of the RINEX observation files.
    # -> The routine will also print out an error log for the user to read.
    
    if timecheck == False:
        return None

    # EXTRACTING LEOGPS SCENARIO START, STOP, STEP, AND RINEX STEP TIMES:
    # ===================================================================
    # -> If all is good, then 'timecheck' outputs start, stop and step times.
    # -> This will be the start, stop, and step used for the entire scenario.
    tstart, tstop, tstep, rnxstep = timecheck

    # EXTRACTING ALL GPS EPHEMERIS AND CLOCK BIAS DATA FROM CDDIS.
    # ============================================================
    # -> The routine 'gpsxtr.py' downloads GPS clock and ephemeris data.
    # -> The routine then interpolates positions and clock biases.
    # -> It estimates GPS satellite velocities using a first order derivative.
    # -> It also outputs GPS satellites without clock issues (goodsats).
    gpsdata, goodsats = gpsxtr.gpsxtr( inps, tstart, tstop, tstep )
    
    # RINEX OBSERVABLES ARE EXTRACTED HERE, IN A NESTED DICTIONARY OBJECT.
    # ====================================================================
    # -> RINEX observables are extracted, screened, and marked for cycle slips.
    # -> If no Doppler values found, they will be estimated in 'dopest.py'.
    # -> Hatch filtering is an option, if the user toggles it in 'config.txt'.
    # -> For more details on the RINEX parsing process, checkout 'rinxtr.py'.
    rinex1 = rinxtr.rinxtr(rinex1f, inps, goodsats, tstart, tstop, rnxstep)
    rinex2 = rinxtr.rinxtr(rinex2f, inps, goodsats, tstart, tstop, rnxstep)
    
    # CHECK RINEX FILES FOR CORRUPTED DATA.
    # =====================================
    if rinex1 == False or rinex2 == False:
        print('Possible corrupted or mis-formatted RINEX file(s) detected!')
        return None
    
    # READ THE DESIRED COORDINATE REFERENCE FRAMES FOR OUTPUT DATA.
    # =============================================================
    frameOrb = inps['frameOrb']
    frameForm = inps['frameForm']
    
    # INITIALISE THE ARRAY OF TIME STEPS USED FOR PROCESSING.
    # =======================================================
    ti, time, results = tstart, [], {}
    while ti <= tstop:
        time.append(ti)
        ti = ti + tstep
    
    # SOLVING FOR POSITIONS (ABSOLUTE AND RELATIVE).
    # ==============================================
    # -> This for loop iterates through all steps generated by the list above.
    # -> In doing so, it refers to RINEX observables generated by 'rinxtr.py'.
    # -> Then, the routine will trilaterate absolute positions, estimate
    # -> receiver biases, offset other errors such as Earth's rotation and
    # -> relativistic effects, and finally estimate the velocities via Doppler.
    # -> The routine 'ambest.py' will then calculate the precise baseline
    # -> using the double-differenced baseline solution, with float ambiguity.
    
    print('Solving for absolute and relative positions now.')
    
    for t in time:
        
        # What are the GPS ephemeris and RINEX observables at t?
        gps = gpsdata[t]
        rx1 = rinex1[t]
        rx2 = rinex2[t]
        
        # Extract PVT, DOP and clock bias values for LEO 1, at time = t.
        pos1, vel1, dop1, cb1 = posvel.posvel(t, goodsats, gps, rx1, inps, n1)
        
        # Extract PVT, DOP and clock bias values for LEO 2, at time = t.
        pos2, vel2, dop2, cb2 = posvel.posvel(t, goodsats, gps, rx2, inps, n2)
        
        # Perform double-differencing to get baseline vector in ITRF.
        baseline = ambest.ambest(t, gps, rx1, rx2, pos1, pos2, inps)
        
        # By default, since the GPS ephemeris provided by CODE is given in
        # ITRF, then the resultant frame is ITRF. If the user selects ICRF
        # then we will perform the coordinate frame computation below.
        # Note that we should account for the possibility that the velocity
        # vector is null, due to a possible lack of Doppler observations.
        if frameOrb == 'ICRF':
            pos1, vel1 = frames.itrf2icrf( t, pos1, vel1 )
            pos2, vel2 = frames.itrf2icrf( t, pos2, vel2 )
        
        # By default, the coordinate reference frame for the relative baseline
        # is also in ITRF. LEOGPS allows for the option of ICRF or the local
        # Euler-Hill frame, taking LEO 1 as the chief satellite.
        if frameForm == 'ICRF':
            baseline, null = frames.itrf2icrf( t, baseline ) # Ignore null
        elif frameForm == 'Hill':
            baseline, null = frames.itrf2icrf( t, baseline ) # Ignore null
            baseline       = frames.icrf2hill( baseline, pos1, vel1 )
        
        # Log the results into a dictionary.
        # Result size: 1x3   1x3   1x3  1x1   1x3   1x3   1x3  1x1       1x3
        results[t] = [pos1, vel1, dop1, cb1, pos2, vel2, dop2, cb2, baseline]
        
    # FINALLY, WE PUBLISH THE RESULTS IN THE OUTPUTS FOLDER.
    # ======================================================
    pubplt.leo_results( results, inps )
    
    return results, inps
