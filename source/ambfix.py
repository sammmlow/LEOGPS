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
##    This is the classical LAMBDA method that was originally authored by    ##
##    Teunissen, Jonge, and Tiberius (1993). The code was later written      ##
##    in MATLAB by Dr Sandra Verhagen and Dr Bofeng Li. It takes in a        ##
##    vector of float ambiguities to the integer least-squares problem,      ##
##    and covariance of the float ambiguities. It then runs the LAMBDA's     ##
##    ILS search-&-shrink and spits out the ambiguity integers. The other    ##
##    5 methods in original LAMBDA MATLAB code are not supported here        ##
##    (feel free to edit the code and implement it youself!). The default    ##
##    ncands = 2, as per original code. All support functions from the       ##
##    original MATLAB code (decorrel, ldldecom) have been nested within      ##
##    the main function as sub functions.                                    ##
##                                                                           ##
##    Inputs:                                                                ##
##    - ahat    : numpy array of float ambiguities                           ##
##    - Qahat   : numpy covariance matrix for float ambiguities              ##
##    - ncands  : number of candidates (optional parameter, default = 2)     ##
##                                                                           ##
##    Outputs:                                                               ##
##    - afixed  : Array of size (n x ncands) with the estimated integer      ##
##                candidates, sorted according to the corresponding squared  ##
##                norms, best candidate first.                               ##
##    - sqnorm  : Distance between integer candidate and float ambiguity     ##
##                vectors in the metric of the variance-covariance matrix.   ##
##                                                                           ##
##    Some notes on translating the Matlab LAMBDA to Python:                 ##
##     - Everything is identical EXCEPT MATLAB is ones-based indexing.       ##
##     - Python is zeros-based indexing, and range function does not         ##
##       include the upper limit index. Thus, only indices have changed.     ##
##     - Example in MATLAB: for i = 1:5 => {1,2,3,4,5}                       ##
##     - Equivalently in Python: for i in range(0,5) => {0,1,2,3,4}          ##
##     - Indices are thus updated accordingly.                               ##
##                                                                           ##
##    Original Developer: Professor Peter Teunissen (TU Delft)               ##
##    Original Authors: Dr Sandra Verhagen and Dr Bofeng Li (TU Delft)       ##
##                                                                           ##
##    Author (Modified) by Samuel Y. W. Low, with permissions.               ##
##    Last modified 08-Jun-2021.                                             ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##    Documentation: https://leogps.readthedocs.io/en/latest/                ##
##                                                                           ##
###############################################################################
###############################################################################

import numpy as np

def LAMBDA( ahat, Qahat, ncands = 2 ):
    '''Integer least-squares method with search-and-shrink for integer
    estimation based on the provided float ambiguity vector (Nx1) ahat and
    associated variance-covariance matrix (NxN) Qahat.
    
    Parameters
    ----------
    ahat : numpy.ndarray
        N-length array of float ambiguities
    Qahat : numpy.ndarray
        NxN covariance matrix of ambiguities
    ncands : int, optional
        Number of search candidates (default = 2)

    Returns
    -------
    afixed : 
        (N x ncands) Array of with estimated integer candidates, sorted 
        according to the corresponding squared norms, best candidate first.
    sqnorm : 
        (ncands x 1) Distance between integer candidate and float ambiguity 
        vectors in the metric of the variance-covariance matrix Qahat.
        
    '''
    
    ###########################################################################
    ###########################################################################
    
    # [afixed, sqnorm] = LAMBDA( ahat, Qahat, ncands )
    #
    # This is the main routine of the LAMBDA software package. By default the
    # ILS method will be used for integer estimation based on the provided
    # float ambiguity vector ahat and associated variance-covariance matrix
    # Qahat. In this Pythonic version (modified by Samuel Low, 2019), only
    # the ILS method is implemented. For other techniques: integer rounding,
    # bootstrapping or Partial Ambiguity Resolution (PAR), the user is free
    # to modify this code and adapt it to their own needs.
    #
    # NOTE 1: LAMBDA always first applies a decorrelation before the integer
    # estimation (for ILS this is required to guarantee an efficient search,
    # for rounding and bootstrapping it is required in order to get higher
    # success rates).
    #
    # INPUTS:
    #
    #     ahat: Float ambiguities (must be a column!)
    #    Qahat: Variance/covariance matrix of ambiguities
    #   ncands: number of candidates (optional parameter, default = 2)
    #
    # OUTPUTS:
    #
    #   afixed: Array of size (n x ncands) with the estimated integer
    #           candidates, sorted according to the corresponding squared
    #           norms, best candidate first.
    #   sqnorm: Distance between integer candidate and float ambiguity vectors 
    #           in the metric of the variance-covariance matrix Qahat.
    #           Only available for ILS.
    #
    # -------------------------------------------------------------------------
    # Release date  : 1-SEPT-2012                                          
    # Authors       : Bofeng LI and Sandra VERHAGEN                                             
    #  
    # GNSS Research Centre, Curtin University
    # Mathematical Geodesy and Positioning, Delft University of Technology                              
    # -------------------------------------------------------------------------
    #
    # REFERENCES: 
    #  1. LAMBDA Software Package: Matlab implementation, Version 3.0.
    #     Documentation provided with this software package.
    #  2. Teunissen P (1993) Least-squares estimation of the integer GPS
    #     ambiguities. In: Invited lecture, section IV theory and methodology,
    #     IAG General Meeting, Beijing, China
    #  3. Teunissen P (1995) The least-squares ambiguity decorrelation
    #     adjustment: a method for fast GPS ambiguity estitmation. J Geod
    #     70:651-7
    #  4. De Jonge P, Tiberius C (1996) The LAMBDA method of intger ambiguity 
    #     estimation:implementation aspects.
    #  5. Chang X ,Yang X, Zhou T (2005) MLAMBDA: a modified LAMBDA method for
    #     integer least-squares estimation
    
    ###########################################################################
    ###########################################################################
    
    # A function for obtaining the decimals only from float arrays.
    
    def floatrem( fltarray ):
        
        # This function is NECESSARY because of the differences between:
        # MATLAB's rem function
        # (computes the true mathematical remainder)
        # And Python's modulo % operator
        # (computes remainder complementary to the floor_divide function)
        
        fltarray = np.array(fltarray)
        fltarray = fltarray + 0.000001
        intarray = fltarray.astype(int)
        decarray = fltarray - intarray
        
        return decarray, intarray
    
    ###########################################################################
    ###########################################################################
    
    # A function to perform LtDL decomposition of the covariance matrix.
        
    def ldldecom( Qahat1 ):
        
        #  This routine finds the LtDL decomposition of a given variance or
        #  covariance matrix.
        # 
        #  Input arguments:
        #     Qahat: Symmetric n by n matrix to be factored
        # 
        #  Output arguments:
        #     L:     n by n factor matrix (strict lower triangular)
        #     D:     Diagonal n-vector
        #  ------------------------------------------------------------------
        #  File.....: ldldecom
        #  Date.....: 19-MAY-1999
        #  Author...: Peter Joosten
        #             Mathematical Geodesy and Positioning
        #             Delft University of Technology
        #  ------------------------------------------------------------------
        
        Qahat2 = Qahat1.copy() 
        
        # If we do not use copy, we will overwrite the original Qahat...
        # ... even the one outside the function! This doesn't occur in MATLAB.
        
        n = len(Qahat2)
        
        D = np.zeros((n))
        
        L = np.zeros((n,n))
        
        for i in range(n-1,-1,-1):
            
            D[i] = Qahat2[i][i]
            
            L[i,0:i+1] = Qahat2[i,0:i+1] / ((Qahat2[i][i])**0.5)

            for j in range(0,i):
                
                Qahat2[j,0:j+1] = Qahat2[j,0:j+1] - L[i,0:j+1]*L[i][j]
                
            L[i,0:i+1] = L[i,0:i+1] / L[i][i]
        
        return L,D
    
    ###########################################################################
    ###########################################################################
    
    # Decorrelation function for LAMBDA.
    
    def decorrel( ahat, Qahat ):
                
        # function [Qzhat,Z,L,D,zhat,iZt] = decorrel (Qahat,ahat)
        # DECORREL: Decorrelate a (co)variance matrix of ambiguities
        # 
        # [Qzhat,Z,L,D,zhat] = decorrel (Qahat,ahat)
        # 
        # This routine creates a decorrelated Q-matrix, by finding the
        # Z-matrix and performing the corresponding transformation.
        # 
        # The method is described in:
        # The routine is based on Fortran routines written by Paul de Jonge
        # and on Matlab-routines written by Kai Borre.
        # The resulting Z-matrix can be used as follows:
        # zhat = Zt * ahat; \hat(z) = Z' * \hat(a);
        # Q_\hat(z) = Z' * Q_\hat(a) * Z
        #
        # Input arguments:
        #   Qahat: Variance-covariance matrix of ambiguities (original)
        #   ahat:  Original ambiguities (optional)
        # 
        # Output arguments:
        #   Qzhat: Variance-covariance matrix of decorrelated ambiguities
        #   Z:     Z-transformation matrix
        #   L:     L matrix (from LtDL-decomposition of Qzhat)
        #   D:     D matrix (from LtDL-decomposition of Qzhat)
        #   zhat:  Transformed ambiguities (optional)
        #   iZt:   inv(Z')-transformation matrix
        #
        # ------------------------------------------------------------------
        # Function.: decorrel
        # Date.....: 19-MAY-1999 / modified 12-APRIL-2012
        # Author...: Peter Joosten / Sandra Verhagen
        #            Mathematical Geodesy and Positioning
        #            Delft University of Technology
        # Modified.: Samuel Low, July 2019, DSO National Laboratories
        # ------------------------------------------------------------------
    
        # Initialisations
        n = len(Qahat)
        iZt = np.identity(n)
        i1 = n - 1
        sw = True
        
        # LtDL decomposition
        L, D = ldldecom(Qahat)
        
        while sw == 1:
            
            i = n # Loop for column from n to 1
            sw = 0
            
            while sw == 0 and i > 1:
                
                i = i - 1 # The i-th column
                
                if i <= i1:
                    
                    for j in range(i,n):
                        
                        # We have to do some manual coding here, as python's
                        # rounding for .5's are different from MATLAB's
                        
                        mu = L[j,i-1] # Get the float mu
                        
                        mu_dec = mu%1 # Get the decimal float of mu
                        
                        if mu_dec == 0.5:
                            
                            mu += 0.01 # Just to make it round up properly.
                            
                        mu = round(mu)
                        
                        if mu != 0.0:
                            
                            L[j:n,i-1] = L[j:n,i-1] - mu * L[j:n,j]
                            
                            iZt[:,j] = iZt[:,j] + mu * iZt[:,i-1]
                            
                delta = D[i-1] + (L[i,i-1]**2) * D[i]
                
                if delta < D[i]:
                    
                    lam = D[i] * L[i,i-1] / delta
                    
                    eta = D[i-1] / delta
                    
                    D[i-1] = eta * D[i]
                    
                    D[i] = delta
                    
                    mult1 = np.array([-1*L[i,i-1], 1])
                    
                    mult2 = np.array([eta,lam])
                    
                    mult3 = np.stack((mult1,mult2))
                    
                    L[i-1:i+1,0:i-1] = np.matmul(mult3,L[i-1:i+1,0:i-1])
                    
                    L[i,i-1] = lam
                    
                    # Flip rows i and i+1
                    
                    L[i+1:n,i-1:i+1] = np.flip(L[i+1:n,i-1:i+1], axis=0)
                    
                    iZt[:,i-1:i+1] = np.flip(iZt[:,i-1:i+1], axis=0)
                    
                    i1 = i
                    
                    sw = 1
        
        iZt = iZt + 0.000001 # Resolves Python 3's rounding definition
        
        Z = np.round(np.linalg.inv(iZt.transpose()))
        
        Qzhat = np.matmul( Qahat, Z )
        
        Qzhat = np.matmul( Z.transpose(), Qzhat )
        
        zhat = np.matmul(Z.transpose(),ahat)
        
        iZt = np.round(iZt)
        
        return Qzhat, Z, L, D, zhat, iZt
    
    ###########################################################################
    ###########################################################################
    
    def ssearch( ahat, L, D, ncands):
        
        #------------------------------------------------------------------|
        #
        #   Integer ambiguity vector search via search-and-shrink technique.
        #
        #   INPUTS:
        #
        #       ahat  : Float ambiguities (should be decorrelated for
        #                computational efficiency)
        #       L,D   : LtDL-decomposition of the variance-covariance matrix 
        #                of the float ambiguities ahat
        #       ncands: Number of requested candidates
        #
        #   OUTPUTS:
        #
        #       afixed: estimated integers (n, x, ncands)
        #       sqnorm: corresponding squared norms (n-vector, ascending order)
        #
        #------------------------------------------------------------------|
        # Date    : 02-SEPT-2010                                           |
        # Author  : Bofeng LI                                              |
        #           GNSS Research Center, Department of Spatial Sciences   |
        #           Curtin University of Technology                        |
        # E-mail  : bofeng.li@curtin.edu.au                                | 
        #------------------------------------------------------------------|
        
        # First, check that float ambiguity and D have same length
        
        if len(ahat) != len(D):
            
            print('Error! Float ambiguity vector must be a column vector!')
            print('It must also have the same dimension as D')
            return None
        
        # Initialising outputs
        n = len(ahat)
        afixed = np.zeros((n, ncands))
        sqnorm = np.zeros(ncands)
        
        # Initializing the variables for searching
        Chi2 = 1.0e+18 # Start search with an infinite chi-square
        dist = np.zeros(n) # MATLAB distance function 
        endsearch = False # Search trigger
        count = 0 # Count the number of candidates

        acond = np.zeros(n)
        acond[n-1] = ahat[n-1]
        
        zcond = np.zeros(n)
        zcond[n-1] = np.round(acond[n-1]+0.000001)
        
        left = acond[n-1] - zcond[n-1]
        
        step = np.zeros(n)
        step[n-1] = np.sign(left)
        
        if step[n-1] == 0:
            step[n-1] = 1 # Give a positive step.
            
        imax = ncands - 1 # Initially, the maximum F(z) is at ncands
        
        S = np.zeros((n,n)) # Used to compute conditional ambiguities
        
        k = n
        
        # Now we start the main search loop.
        
        while endsearch == False:
            
            newdist = dist[k-1] + (left**2) / D[k-1]
                        
            if newdist < Chi2:
                
                if k != 1: # Case 1: move down
                    
                    k -= 1
                    
                    dist[k-1] = newdist
                    
                    S[k-1,0:k] = S[k,0:k] + (zcond[k] - acond[k])*L[k,0:k]
                    
                    acond[k-1] = ahat[k-1] + S[k-1,k-1]
                    
                    zcond[k-1] = np.round(acond[k-1]+0.000001)
                    
                    left = acond[k-1] - zcond[k-1]
                    
                    step[k-1] = np.sign(left)
                    
                    if step[k-1] == 0: # Very rarely would this happen...
                        
                        step[k-1] = 1 # ... but just in case, you know.
                        
                else: # Case 2: store the found candidate and try the next.
                    
                    if count < (ncands - 1):
                        
                        # Store the 1st ncands-1 initial points as candidates
                        
                        count += 1
                        
                        afixed[:,count-1] = zcond[0:n];
                        
                        sqnorm[count-1] = newdist # Store F(zcond)
                        
                    else:
                        
                        afixed[:,imax] = zcond[0:n]
                        
                        sqnorm[imax] = newdist
                        
                        Chi2 = max(sqnorm)
                        
                        imax = np.argmax(sqnorm) # No need to add '-1' to imax
                        
                    zcond[0] = zcond[0] + step[0]
                    
                    left = acond[0] - zcond[0]
                    
                    step[0] = -1*step[0] - np.sign(step[0])
                    
            else: # Case 3: exit or move up
                
                if k == n:
                    
                    endsearch = True
                    
                else:
                    
                    k += 1 # Move up
                    
                    zcond[k-1] = zcond[k-1] + step[k-1]
                    
                    left = acond[k-1] - zcond[k-1]
                    
                    step[k-1] = -1*step[k-1] - np.sign(step[k-1])
                
        order = np.argsort(sqnorm) # Get an array of INDICES for a sort.
        
        sqnormf = np.sort(sqnorm) # Get an array of ACTUAL SORTS for sqnorm.
        
        afixedf = np.copy(afixed)
        
        for k in range(0,len(order)):
            
            afixedf[:,k] = afixed[:,order[k]]
            
        return afixedf, sqnormf
    
    ###########################################################################
    ###########################################################################
    
    # Test inputs: Is the Q-matrix symmetric?
    if np.array_equal(Qahat,Qahat.transpose()) == False:
        print('Variance-covariance matrix is not symmetric!')
        return None
    
    # Test inputs: Is the Q-matrix positive-definite?
    if np.sum(np.linalg.eig(Qahat)[0] > 0.0) != len(Qahat):
        print('Variance-covariance matrix is not positive definite!')
        return None

    # Test inputs: Does Q-matrix and amb vector have identical dimensions?
    if len(ahat) != len(Qahat):
        print('Variance-covariance matrix and vector of ambiguities...')
        print('... do not have identical dimensions!')
        return None
    
    ###########################################################################
    ###########################################################################
    
    # Begin least-squares ambiguity decorrelation adjustment!
    
    # Remove integer numbers from float solution, so that all values are
    # between -1 and 1 (for computational convenience only)
    ahat, incr = floatrem( ahat )
    
    # Compute Z matrix based on the decomposition  Q=L^T*D*L;
    Qzhat, Z, L, D, zhat, iZt = decorrel( ahat, Qahat )
    
    # Integer ambiguity vector search via search-and-shrink
    zfixedff, sqnormff = ssearch( zhat, L, D, ncands )
    
    # Perform the back-transformation and add the increments
    afixed = np.matmul(iZt,zfixedff)
    repmat = np.repeat(np.array([incr]),ncands,axis=0)
    repmat = repmat.transpose()
    afixed = afixed + repmat
    afixed = afixed.transpose()
    
    ###########################################################################
    ###########################################################################
    
    # Returns best amb-fix, second best amb-fix, and the square norm.
    
    return afixed, sqnormff

    ###########################################################################
    ###########################################################################
