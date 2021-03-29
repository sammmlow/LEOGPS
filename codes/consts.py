#!/usr/bin/env python
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
## No inputs or outpots. Just a List of useful constants.                    ##
## Adapted from (and credits to) Bernese GNSS v5.2                           ##
##                                                                           ##
## AUTHOR MODIFIED: 07-05-2019, by Samuel Y.W. Low                           ##
##                                                                           ##
###############################################################################
###############################################################################
'''

C      = 299792458.0      # VELOCITY OF LIGHT                M/SEC
FREQ1  = 1575420000.0     # L1-CARRIER FREQUENCY   GPS       1/SEC
FREQ2  = 1227600000.0     # L2-CARRIER FREQUENCY   GPS       1/SEC
FREQ5  = 1176450000.0     # L5-CARRIER FREQUENCY   GPS       1/SEC
FREQP  = 10230000.0       # P-CODE     FREQUENCY   GPS       1/SEC
FREQG1 = 1602000000.0     # L1-CARRIER FREQUENCY   GLONASS   1/SEC
FREQG2 = 1246000000.0     # L2-CARRIER FREQUENCY   GLONASS   1/SEC
DFRQG1 = 562500.0         # L1-CARRIER FREQ. DIFF. GLONASS   1/SEC
DFRQG2 = 437500.0         # L2-CARRIER FREQ. DIFF. GLONASS   1/SEC
FREQGP = 5110000.0        # P-CODE     FREQUENCY   GLONASS   1/SEC
FRQE1  = 1575420000.0     # L1-CARRIER FREQUENCY   GALILEO   1/SEC
FRQE5  = 1191795000.0     # L5-CARRIER FREQUENCY   GALILEO   1/SEC
FRQE5a = 1176450000.0     # L5a-CARRIER FREQUENCY  GALILEO   1/SEC
FRQE5b = 1207140000.0     # L5b-CARRIER FREQUENCY  GALILEO   1/SEC
FRQE6  = 1278750000.0     # L6-CARRIER FREQUENCY   GALILEO   1/SEC
FRQS1  = 1575420000.0     # L1-CARRIER FREQUENCY   SBAS      1/SEC
FRQS5  = 1176450000.0     # L5-CARRIER FREQUENCY   SBAS      1/SEC
FRQC1  = 1589740000.0     # L1-CARRIER FREQUENCY   COMPASS   1/SEC
FRQC2  = 1561098000.0     # L2-CARRIER FREQUENCY   COMPASS   1/SEC
FRQC5b = 1207140000.0     # L5b-CARRIER FREQUENCY  COMPASS   1/SEC
FRQC6  = 1268520000.0     # L6-CARRIER FREQUENCY   COMPASS   1/SEC
FRQJ1  = 1575420000.0     # L1-CARRIER FREQUENCY   QZSS      1/SEC
FRQJ2  = 1227600000.0     # L2-CARRIER FREQUENCY   QZSS      1/SEC
FRQJ5  = 1176450000.0     # L5-CARRIER FREQUENCY   QZSS      1/SEC
FRQJ6  = 1278750000.0     # L6-CARRIER FREQUENCY   QZSS      1/SEC
GM     = 398.6004415e12   # GRAVITY CONSTANT*EARTH MASS      M**3/SEC**2
GMS    = 1.3271250e20     # GRAVITY CONSTANT*SOLAR MASS      M**3/SEC**2
GMM    = 4.9027890e12     # GRAVITY CONSTANT*LUNAR MASS      M**3/SEC**2
AU     = 149597870691     # ASTRONOMICAL UNIT                M
AE     = 6378137.0        # EQUATORIAL RADIUS OF EARTH       M
CONRE  = 6371000.0        # MEAN RADIUS OF THE EARTH         M
J2     = 1.0826359e-3     # DYNAMICAL FORM-FACTOR IERS(2003) 1
FACTEC = 40.3e16          # IONOSPHERIC FACTOR               M/SEC**2/TECU
P0     = -0.94e-7         # NOMINAL RAD.PR. ACCELERAT.       M/SEC**2
OMEGA  = 7292115.1467e-11 # ANGULAR VELOCITY OF EARTH        RAD/USEC
EPHUTC = 55.0             # EPH. TIME (ET) MINUS UTC         SEC
WGTPHA = 1.0              # WEIGHT FOR PHASE OBSERVATIONS    1
WGTCOD = 1.0e-4           # WEIGHT FOR CODE OBSERVATIONS     1
HREF   = 0.0              # REFERENCE HEIGHT FOR METEO MODEL M
PREF   = 1013.25          # PRESSURE AT HREF                 MBAR
TREF   = 18.0             # TEMPERATURE AT HREF              DEG. CELSIUS
HUMREF = 50.0             # HUMIDITY AT HREF                 %
ERR    = 7.2921150e-5     # EARTH INERTIAL ROTATION RATE     RAD/SEC