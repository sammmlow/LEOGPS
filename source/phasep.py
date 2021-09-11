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
##    Screens the phase information in RINEX observables (C1,P2,L1,L2).      ##
##    RINEX observables are stored as a dict returned by 'rinxtr.py'.        ##
##                                                                           ##
##    Inputs: Takes in the RINEX data dictionary formed preliminarily        ##
##    from rinxtr.rinxtr(), for marking and cycle slip detection.            ##
##                                                                           ##
##    Outputs: Phase-processed nested dictionary of RINEX observations.      ##
##    The format is essentially the same, except in the third sub-dict       ##
##    rnxdata[epoch][SV] has two new key-value pairs added. First, the       ##
##    L4 intermediate carrier phase value; and second, a cycle slip flag     ##
##    based on L4 observables computed.                                      ##
##                                                                           ##
##    If freqnum == 1 (single frequency):                                    ##
##        L4 = Code-phase Melbourne-Wubbena linear combination               ##
##                                                                           ##
##    If freqnum == 2 (dual frequency)                                       ##
##        L4 = Geometry-free linear combination                              ##
##                                                                           ##
##    Where L4 will be used as a derived carrier phase value for cycle       ##
##    slip detection by checking for phase values beyond a tolerance.        ##
##                                                                           ##
##    The format is formatted as (and will be the format of 'rinxtr.py)      ##
##                                                                           ##
##    Output = {epoch1 : {1 : {'C1':123, 'L1':123, ...                       ##
##                             'L4':321, 'flag':'none'} ...} ...             ##
##              epoch2 : {2 : {'C1':123, 'L1':123, ...                       ##
##                             'L4':321, 'flag':'slip'} ...} ...             ##
##                                                                           ##
##              ... ... ... ... ... ...                                      ##
##                                                                           ##
##              epochX : {1 : {'C1':123, 'L1':123, ...                       ##
##                             'L4':321, 'flag':'none'} ...}}                ##
##                                                                           ##
##    This programme is run as a subroutine in 'rinxtr.py' only.             ##
##    It does not take in any other observation file, except 'config.txt'    ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 07-Jun-2021.                                             ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##    Documentation: https://leogps.readthedocs.io/en/latest/                ##
##                                                                           ##
###############################################################################
###############################################################################

# Import global libraries
import copy
import warnings
import numpy as np

# The first function flags the carrier phase status at each epoch.
def phsmrk(rnxdata, rnxstep, goodsats, inps):
    '''Marks the **rnxdata** dict with carrier phase flags at each epoch.

    Parameters
    ----------
    rnxdata : dict
        A nested dictionary comprising code observations, carrier phase,
        and doppler values.
    rnxstep :datetime.timedelta
        Observed time step in RINEX file
    goodsats : list
        Sorted list of GPS satellites without outages by PRN IDs
    inps : dict
        A dictionary of inputs created by `inpxtr.inpxtr()`

    Returns
    -------
    rnxproc : dict
        A nested dictionary comprising code observations, carrier phase,
        doppler values, with an added L4 phase and a carrier phase flag.
    
    '''
    
    # Get the desired input parameters.
    freqnum = inps['freq']
    cycleslip_tolerance = float(inps['cycsliptol'])
    cycleslip_filtlength = inps['cycsliplen']
    
    # Ignore polyfit warnings
    warnings.simplefilter('ignore', np.RankWarning)
    
    # Wavelengths of L1 and L2
    WL1 = 0.190293672798
    WL2 = 0.244210213425
    
    print('Marking and screening observables for cycle slips.')
    
    rnxproc = copy.deepcopy(rnxdata) # Dictionary of processed RINEX data
    slipcount = 0 # Counting the number of cycle slips

    # For each particular SV ID
    for SV in goodsats:
        
        # Across all the epochs recorded
        for epoch in rnxproc:
            
            # Initialise a time variable for the previous and next step
            prev_epoch = epoch-rnxstep
            next_epoch = epoch+rnxstep
                        
            # Now, we check if the current SV is being observed, and flag it.
            if SV in rnxproc[epoch]:
                
                # Check if this is the first observation.
                if prev_epoch in rnxproc:
                    if SV not in rnxproc[prev_epoch]:
                        rnxproc[epoch][SV]['flag'] = 'start'
                else:
                    rnxproc[epoch][SV]['flag'] = 'start'
                
                # Check if this is the final observation.
                if next_epoch in rnxproc:
                    if SV not in rnxproc[next_epoch]:
                        rnxproc[epoch][SV]['flag'] = 'end'
                else:
                    rnxproc[epoch][SV]['flag'] = 'end'
                
                # Check if this is a lone observable
                if prev_epoch in rnxproc and next_epoch in rnxproc:
                    if SV not in rnxproc[prev_epoch]:
                        if SV not in rnxproc[next_epoch]:
                            rnxproc[epoch][SV]['flag'] = 'solo'
                    if SV in rnxproc[prev_epoch]:
                        if SV in rnxproc[next_epoch]:
                            rnxproc[epoch][SV]['flag'] = 'none'
                elif prev_epoch not in rnxproc and next_epoch in rnxproc:
                    if SV not in rnxproc[next_epoch]:
                        rnxproc[epoch][SV]['flag'] = 'solo'
                elif next_epoch not in rnxproc and prev_epoch in rnxproc:
                    if SV not in rnxproc[prev_epoch]:
                        rnxproc[epoch][SV]['flag'] = 'solo'
                
                # At this stage, we will calculate an intermediate value L4.
                # L4 will then be poly-fitted with adjacent values across time.
                # If L4 remains an outlier from the polynomial fit,
                # Then that epoch's L4 value identifies as a cycle slip.
                
                # Perform Melbourne-Wubbena linear combination as the L4.
                if freqnum == 1:
                    if 'P1' in rnxproc[epoch][SV]:
                        L1p = rnxproc[epoch][SV]['P1'] # P1 code observable
                    elif 'C1' in rnxproc[epoch][SV]:
                        L1p = rnxproc[epoch][SV]['C1'] # P1 code observable
                    else:
                        print('Problem found, no C1 or P1 observable in: ')
                        print(str(epoch) + ' for satellite SV ' + str(SV))
                        return None
                    L2p = rnxproc[epoch][SV]['L1'] # L1 phase observable
                    L2p = L2p * WL1 # Multiply phase by wavelength
                
                # Perform the geometry-free linear combination as the L4.
                elif freqnum == 2:
                    L1p = rnxproc[epoch][SV]['L1'] # L1 phase observable
                    L1p = L1p * WL1 # Multiply phase by wavelength
                    L2p = rnxproc[epoch][SV]['L2'] # L2 phase observable
                    L2p = L2p * WL2 # Multiply phase by wavelength
                    
                # Conclude the entry for the geometry-free LC and the flag
                L4 = L1p - L2p # Time-stamped geometry-free LC for dual freq
                rnxproc[epoch][SV]['L4'] = L4 # Throw in L4 value
    
    # Now we begin the process of an N-window polynomial fitting filter
    # We will use a quadratic filter for this approach
    
    N = cycleslip_filtlength # Length of sliding window filter for poly-fit
    
    # For each SV ID
    for SV in goodsats:
        
        k = 1 # Restart the time counter
        t = [] # Restart the time array
        L = [] # Restart the observation array
        
        # Then check through each epoch based on a SV ID
        for epoch in rnxproc:
            
            prev_epoch = epoch - rnxstep
            
            # If this SV exists
            if SV in rnxproc[epoch]:
                
                # Read the flag of this SV observable
                obs = rnxproc[epoch][SV]['L4']
                flag = rnxproc[epoch][SV]['flag']
                
                if prev_epoch in rnxproc:
                    if SV in rnxproc[prev_epoch]:
                        preflag = rnxproc[prev_epoch][SV]['flag']
                else:
                    preflag = 'none'
                
                # Now is the critical steps that lead to the poly-fit...
                # We want to ensure that L4 observations are recorded for
                # Flags: start, none, end...
                # We would never have a case where the current flag is slip
                
                # First, let's check if it is the last entry of this SV...
                if flag == 'end':
                    
                    trigger = True # Trigger to do polynomial fitting
                    k += 1 # Increase the time stamp by 1
                    t.append(k) # Add in the first element
                    L.append(obs) # Add in the L4 observation
                
                # Next, let's see if this is a stand-alone observable.
                elif flag == 'solo':
                    
                    trigger = False # Do not record the poly-fit results.
                    k = 1 # Restart the time counter
                    t = [] # Restart the time array
                    L = [] # Restart the observation array
                
                # If it isn't, then let's check if it is the starting entry
                elif flag == 'start' or preflag == 'slip':
                    
                    trigger = False # Reset the poly-fit trigger
                    k = 1 # Restart the time counter
                    t = [] # Restart the time array
                    L = [] # Restart the observation array
                    
                    t.append(k) # Add in the first element
                    L.append(obs) # Add in the L4 observation
                
                # Otherwise then it is just a normal entry...
                elif flag == 'none':
                    
                    k += 1 # Increase the time stamp by 1
                    t.append(k) # Add in the first element
                    L.append(obs) # Add in the L4 observation
                    
                elif flag == 'slip':
                    print('Something is wrong here...')
                    print('Current epoch is recorded as a cycle slip?')
                    print('Epoch and SV ID:')
                    print(str(epoch) + ' ' + str(SV))
                    
                else:
                    print('Something else is wrong here...')
                    print('Flag string does not match any known flags...')
                    print('Current flag is: ' + str(flag))
                
                # Be mindful that the tolerance depends on the time step
                # Due to the time-dependence on the ionospheric variation
                # Now, is where we decide if we wish to do poly-fitting
                
                if k == N+1: # It has reached the filter size limit
                    trigger = True
                
                # If the trigger is true...
                # Then this is where we screen for cycle slips
                if trigger == True:
                    
                    # Over here, we do the polynomial curve fitting
                    tn = np.array(t) # Numpy-rize the time array
                    Ln = np.array(L) # Numpy-rize the original L4 data
                    pn = np.polyfit(tn,Ln,2)
                    Lp = np.polyval(pn,tn) # Polynomial
                    Lres = Ln - Lp # Get the residual values
                    Lstd = np.std(Lres) # Get the standard deviation
                    
                    # If the observed minus computed (O-C) value:
                    # Exceeds 'cycsliptol' number of standard deviations,
                    # Then, declare it as a cycle slip and change the flag
                    
                    if np.abs(Lres[-1]) > cycleslip_tolerance * Lstd:
                        
                        # Declare this SV at this epoch as a cycle slip!
                        # However, there may be edge cases where it is both:
                        # Flag = 'slip' and flag = 'end'.
                        # In this case, 'slip' should take priority.
                        
                        rnxproc[epoch][SV]['flag'] = 'slip'
                        
                        if prev_epoch in rnxproc:
                            if SV in rnxproc[prev_epoch]:
                                rnxp = rnxproc[prev_epoch][SV]
                                if rnxp['flag'] == 'start':
                                    rnxp['flag'] = 'solo'
                                elif rnxp['flag'] == 'none':
                                    rnxp['flag'] = 'end'
                                
                            
                        if next_epoch in rnxproc:
                            if SV in rnxproc[next_epoch]:
                                rnxn = rnxproc[next_epoch][SV]
                                if rnxn['flag'] != 'end':
                                    rnxn['flag'] = 'start'
                                
                        slipcount += 1
            
            # Maintain the filter window at length = N
            # Discard the old value of L4 observation
            if k == N+1:
                k -= 1
                del t[-1]
                del L[0]
    
    # At this junction, cycle slips have already been flagged!
    # We will return this dictionary 'rnxproc'.
    # Typically, after this step, hatch filtering will commence.
    
    print('Phase screening of cycle slips done!')
    print('Total number of cycle slips: ' + str(slipcount) + '\n')
    
    return rnxproc

# The single-frequency code-carrier hatch filter. 
def hatch1(rf,ri,pf,pi,M,ltype):
    '''Single-frequency code-carrier hatch filter for a single data point. 
    
    Parameters
    ----------
    rf : float
        Pseudorange at epoch [i]
    ri : float
        Pseudorange at epoch [i-1]
    pf : float
        Carrier phase at epoch [i]
    pi : float
        Carrier phase at epoch [i-1]
    M : int
        Hatch filter length
    ltype : str
        Frequency string ('L1' or 'L2')

    Returns
    -------
    float
        Smoothed code-phase value.

    '''

    # r1,r2 are the pseudoranges of epoch (i) and (i-1) respectively
    # p1,p2 are the carrier phases of epoch (i) and (i-1) respectively
    # M is the filter length
    
    C = 299792458.0
    
    if ltype == 'L1':
        freq = 1575420000.0
    elif ltype == 'L2':
        freq = 1227600000.0
    else:
        print('Invalid carrier phase observation type!')
        return None
    
    wavelength = C/freq
    pf = pf*wavelength
    pi = pi*wavelength
    
    return (rf/M) + ((M-1)/M)*(ri+pf-pi)

# The dual-frequency divergence free code-carrier hatch filter.
def hatch2(rf,ri,pf1,pf2,pi1,pi2,M):
    '''Dual-frequency code-carrier hatch filter for a single data point. 
    
    Parameters
    ----------
    rf : float
        Pseudorange at epoch [i]
    ri : float
        Pseudorange at epoch [i-1]
    pf1 : float
        Carrier phase L1 at epoch [i]
    pf2 : float
        Carrier phase L2 at epoch [i]
    pi1 : float
        Carrier phase L1 at epoch [i-1]
    pi2 : float
        Carrier phase L2 at epoch [i-1]
    M : int
        Hatch filter length

    Returns
    -------
    float
        Smoothed code-phase value.

    '''
    
    # r1,r2 are the pseudoranges of epoch (i) and (i-1) respectively
    # p1,p2 are the carrier phases of epoch (i) and (i-1) respectively
    # M is the filter length
    
    C = 299792458.0
    freq1 = 1575420000.0
    freq2 = 1227600000.0
    
    wavelength1 = C/freq1
    wavelength2 = C/freq2
    pi1 = pi1*wavelength1
    pi2 = pi2*wavelength2
    pf1 = pf1*wavelength1
    pf2 = pf2*wavelength2
    
    # Divergence free phase value
    pdfi = pi1 + 3.091455560326319*(pi1-pi2)
    pdff = pf1 + 3.091455560326319*(pf1-pf2)
    
    return (rf/M) + ((M-1)/M)*(ri+pdff-pdfi)

# Code-carrier smoothing for single-frequency (standard Hatch filter).
def ph1fil(rnxdata, rnxstep, goodsats, hatchlen):
    '''Single-frequency hatch filter to be applied to the entire RINEX data
    dictionary, for carrier phase values without cycle slips. This function
    calls hatch1() repeatedly in the filter loop.
    
    Parameters
    ----------
    rnxdata : dict
        A nested dictionary comprising code observations, carrier phase,
        and doppler values.
    rnxstep :datetime.timedelta
        Observed time step in RINEX file
    goodsats : list
        Sorted list of GPS satellites without outages by PRN IDs
    hatchlen : int
        Length of the hatch filter
        
    Returns
    -------
    rnxdata : dict
        A nested dictionary comprising code observations, carrier phase,
        and doppler values, with the code values smoothed by Hatch filtering.
    
    '''
    
    print('Performing single-frequency L1 Hatch filtering for smoothing.')
        
    # Initialise a variable that checks if there was a cycle slip before
    
    slipcheck = False # False is default => no cycle slip in previous epoch
    
    # Now, we begin the hatch filtering process by:
    # First, parsing across time, for each SV.
    # We will smooth values between start/end.
    # Cycle slips will also interrupt the sliding of the window filter

    for SV in goodsats:
        
        k = 0 # Initialise a filter length k
        
        for epoch in rnxdata:
            
            if SV in rnxdata[epoch]:
                
                trigger = True # Allow for code-carrier Hatch filtering
                flag = rnxdata[epoch][SV]['flag'] # Check the flag
                
                # We also want to check if it is P1 or C1 entry
                if 'P1' in rnxdata[epoch][SV]:
                    obstype = 'P1'
                    Ltype = 'L1'
                elif 'C1' in rnxdata[epoch][SV]:
                    obstype = 'C1'
                    Ltype = 'L1'
                else:
                    print('Something is wrong. No valid code observables')
                    return None
            
            else: # Ignore this entry if there is no SV
                
                trigger = False 
            
            # Now, if the trigger is true, we will begin Hatch filtering
                
            if trigger == True:
                
                if flag == 'solo':
                    
                    trigger = False
                    
                elif flag == 'start':
                    
                    k = 0 # Initialise a filter length k
                
                elif flag == 'none' or flag == 'end':
                    
                    if slipcheck == False:
                    
                        if k > hatchlen:
                            k -= 1
                        
                        for n in range(0,k,-1): # Sliding window hatch filter
                            
                            # Get current and previous code-carrier values
                            R1 = rnxdata[epoch+((n)*rnxstep)][SV][obstype]
                            R2 = rnxdata[epoch+((n-1)*rnxstep)][SV][obstype]
                            P1 = rnxdata[epoch+((n)*rnxstep)][SV][Ltype]
                            P2 = rnxdata[epoch+((n-1)*rnxstep)][SV][Ltype]
                            kM = np.abs(k) + 1
                            
                            # Generate the new smoothed code value
                            RNew = hatch1(R1,R2,P1,P2,kM,Ltype)
                            rnxdata[epoch+((n)*rnxstep)][SV][obstype] = RNew
                    
                    else: # Reset the cycle slip trigger if it was True.
                        slipcheck = False
                        continue
                                            
                elif flag == 'slip':
                    
                    slipcheck = True # Cycle slip check triggered
                    k = 0 # Restart the filtering process
                
                else:
                    print('Something is wrong. There is an empty flag entry')
                    print(str(epoch) + ' @ SV ' + str(SV))
                    return None
    
    print('Completed single-frequency L1 Hatch filtering! \n')
    
    return rnxdata

# Code-carrier smoothing for dual frequency divergence free Hatch filter.

def ph2fil(rnxdata, rnxstep, goodsats, hatchlen):
    '''Dual-frequency hatch filter to be applied to the entire RINEX data
    dictionary, for carrier phase values without cycle slips. This function
    calls hatch2() repeatedly in the filter loop.
    
    Parameters
    ----------
    rnxdata : dict
        A nested dictionary comprising code observations, carrier phase,
        and doppler values.
    rnxstep :datetime.timedelta
        Observed time step in RINEX file
    goodsats : list
        Sorted list of GPS satellites without outages by PRN IDs
    hatchlen : int
        Length of the hatch filter
        
    Returns
    -------
    rnxdata : dict
        A nested dictionary comprising code observations, carrier phase,
        and doppler values, with the code values smoothed by Hatch filtering.
    
    '''
    
    print('Performing dual frequency divergence free Hatch filtering')
    
    # Initialise a variable that checks if there was a cycle slip before
    
    slipcheck = False # False is default => no cycle slip in previous epoch
    
    # Now, we begin the hatch filtering process by:
    # First, parsing across time, for each SV.
    # We will smooth values between start/end.
    # Cycle slips will also interrupt the sliding of the window filter

    for SV in goodsats:
        
        k = 0 # Initialise a filter length k
        
        for epoch in rnxdata:
            
            if SV in rnxdata[epoch]:
                
                trigger = True # Allow for code-carrier Hatch filtering
                flag = rnxdata[epoch][SV]['flag'] # Check the flag
                
                # We also want to check if it is P1 or C1 entry
                if 'P1' in rnxdata[epoch][SV]:
                    obs1 = 'P1'
                elif 'C1' in rnxdata[epoch][SV]:
                    obs1 = 'C1'
                else:
                    print('Something is wrong. No valid C1/P1 observables')
                    return None
            
            else: # Ignore this entry if there is no SV
                
                trigger = False 
            
            # Now, if the trigger is true, we will begin Hatch filtering
                
            if trigger == True:
                
                if flag == 'solo':
                    
                    trigger = False
                    
                elif flag == 'start':
                    
                    k = 0 # Initialise a filter length k
                
                elif flag == 'none' or flag == 'end':
                    
                    if slipcheck == False:
                    
                        if k > hatchlen:
                            k -= 1
                        
                        for n in range(0,k,-1): # Sliding window hatch filter
                            
                            # Get current and previous code-carrier values
                            RF = rnxdata[epoch+((n)*rnxstep)][SV][obs1]
                            RI = rnxdata[epoch+((n-1)*rnxstep)][SV][obs1]
                            
                            # Get L1 phase values now and before
                            PF1 = rnxdata[epoch+((n)*rnxstep)][SV]['L1']
                            PI1 = rnxdata[epoch+((n-1)*rnxstep)][SV]['L1']
                            
                            # Get L2 phase values now and before
                            PF2 = rnxdata[epoch+((n)*rnxstep)][SV]['L2']
                            PI2 = rnxdata[epoch+((n-1)*rnxstep)][SV]['L2']
                            
                            # Get length of filter
                            kM = np.abs(k) + 1
                            
                            # Generate the new smoothed code value
                            RNew = hatch2(RF,RI,PF1,PF2,PI1,PI2,kM)
                            rnxdata[epoch+((n)*rnxstep)][SV][obs1] = RNew
                    
                    else: # Reset the cycle slip trigger
                        slipcheck = False                    
                                            
                elif flag == 'slip':
                    
                    slipcheck = True # Cycle slip check triggered
                    k = 0 # Restart the filtering process
                
                else:
                    print('Something is wrong. There is an empty flag entry')
                    print(str(epoch) + ' @ SV ' + str(SV))
                    return None
    
    print('Completed dual-frequency L2 Hatch filtering! \n')
    
    return rnxdata