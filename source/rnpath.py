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
##    This module is a RINEX path finder. If an observation file is          ##
##    found, with the 4-letter spacecraft ID matching the user's inputs      ##
##    in config.txt, the filepath will be forwarded to 'rinxtr.py' and       ##
##    'timing.py' for data extraction. If observation files are found,       ##
##    but they are 'hatanaka' compressed, then this module runs the          ##
##    hatanaka decompression executable to decompress the file, from         ##
##    the (.D) to (.O) extension.                                            ##
##                                                                           ##
##    Yuki Hatanaka (hatagsi.go.jp) (GSI) wrote and maintains rnx2crx        ##
##    and crx2rnx, which allows the user to compress or decompress a         ##
##    RINEX observation file into smaller ASCII format. The Hatanaka-        ##
##    compressed ASCII format version of a RINEX observation file is         ##
##    frequently used in conjuction with the UNIX compress, zip, gzip        ##
##    or other general compression utilities to create a very small          ##
##    file for Internet transfer.                                            ##
##                                                                           ##
##    In LEOGPS v1.1, Martin Valgur contributed the Pythonic translation     ##
##    of the Hatanaka compression and GZIP libraries. This will replace      ##
##    the previous Windows-only executable files for Hatanaka compression,   ##
##    decompression, and unzipping, thereby making LEOGPS friendly on all    ##
##    other non-Windows platforms.                                           ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 14-Apr-2021 (Martin Valgur)                              ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##    Documentation: https://leogps.readthedocs.io/en/latest/                ##
##                                                                           ##
###############################################################################
###############################################################################

# Import global libraries
import hatanaka
from pathlib import Path

def rnpath(inps):
    '''Gets the paths to RINEX files, following AIUB CODE's naming convention.

    Parameters
    ----------
    inps : dict
        A dictionary of inputs created by `inpxtr.inpxtr()`

    Returns
    -------
    rnx1file : pathlib.Path
        Path of RINEX file of LEO-A
    rnx2file : pathlib.Path
        Path of RINEX file of LEO-B
    
    '''
    
    # Let us start by parsing out all relevant user configuration inputs.
    yy       = inps['dtstart_yy'] # The starting two-digit year number
    doy      = inps['dtstart_doy'] # The starting day of year in GPST
    cwd      = Path(inps['cwd']) # Current working directory
    iwd      = cwd / 'input'
    name1    = inps['name1'] # 4-letter ID of the first spacecraft
    name2    = inps['name2'] # 4-letter ID of the second spacecraft

    # Now, we get all necessary file paths, in the input folder.
    rnx1file = iwd / (name1 + doy + '0.' + yy + 'O') # pathlib Path class
    rnx2file = iwd / (name2 + doy + '0.' + yy + 'O') # pathlib Path class
    crx1file = iwd / (name1 + doy + '0.' + yy + 'D') # pathlib Path class
    crx2file = iwd / (name2 + doy + '0.' + yy + 'D') # pathlib Path class
    
    # Check if there is a need for hatanaka decompression for LEO 1.
    if rnx1file.exists():
        print('Decompressed RINEX obs file observed for ' + name1 + '\n')
    elif crx1file.exists():
        print('Hatanaka compressed file observed for ' + name1)
        hatanaka.decompress_on_disk(crx1file)
        if rnx1file.exists():
            print('Hatanaka decompression successful for LEO1! \n')
        else:
            print('Decompression failed. Did you rename any folders? \n')
    else:
        print('Neither compressed nor decompressed RINEX files were found! \n')

    # Check if there is a need for hatanaka decompression for LEO 2.
    if rnx2file.exists():
        print('Decompressed RINEX obs file observed for ' + name2 + '\n')
    elif crx2file.exists():
        print('Hatanaka compressed file observed for ' + name2)
        hatanaka.decompress_on_disk(crx2file)
        if rnx2file.exists():
            print('Hatanaka decompression successful for LEO2! \n')
        else:
            print('Decompression failed. Did you rename any folders? \n')
    else:
        print('Neither compressed nor decompressed RINEX files were found! \n')    

    if rnx1file.exists() and rnx2file.exists():
        return rnx1file, rnx2file
    else:
        print('Error, somehow RINEX observation files still not found!')
        print('Did you rename any folders accidentally? \n')
        return False