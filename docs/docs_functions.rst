..
   ###########################################################################
   ###########################################################################
   ##                                                                       ##
   ##     _    ___  ___   ___ ___ ___                                       ##
   ##    | |  | __ /   \ / __| _ | __|                                      ##
   ##    | |__| __  ( ) | (_ |  _|__ \                                      ##
   ##    |____|___ \___/ \___|_| \___/                                      ##
   ##                                    v 1.2 (Stable)                     ##
   ##                                                                       ##
   ###########################################################################
   ###########################################################################

.. image:: /_static/leogps_logo.png

|

API Reference 
=============

The order of functions in this API reference goes according to the chronological order of which they are called in the native LEOGPS processing work flow (see previous page). Functions that are currently not in use in the current native work flow are listed at the end of this page.

leogui.py
---------

The first item to be called when running **leogps.py** is the root GUI object (tkinter.Tk() object). This object interfaces directly with **config.txt** to save or load your input parameters.

.. automodule:: leogui
   :members:

leorun.py
---------

The **leorun.run()** function is called once the user clicks on **Run LEOGPS**.

.. automodule:: leorun
   :members:

inpxtr.py
---------

The first thing that is run in the native work flow is to extract all the inputs from **config.txt**. This is done through the **inpxtr** module which comprises two functions. First, **inpxtr.inptim()**, extracts the dates and times from the config file and outputs them in GPS time format. Second, **inpxtr.inpxtr()** extracts all the other processing parameters.

.. automodule:: inpxtr
   :members:

rnpath.py
---------

The second step is to search for the RINEX file paths based on the 4-letter IDs of both spacecraft, and to perform Hatanaka decompression if necessary.

.. automodule:: rnpath
   :members:

timing.py
---------

The third step involves deconflicting all the timing parameters. The routine **timing.py** will take in three sets of timings: the user-specified start-stop times in the GUI or in **config.txt**, and the start-stop times from both of the RINEX observation files. This module will then output the intersection of all three timelines. This module will also check if the time steps of the RINEX files and time step specified by the user are compatible.

.. automodule:: timing
   :members:

gpsxtr.py
---------

The fourth step is to extract out the interpolated GPS satellite ephemeris and clock data.

The primary output of **gpsxtr.gpsxtr()** is the **gpsdata** dictionary, which is a three-layer nested dictionary comprising all of the Lagrange-interpolated ephemeris and clock data. The first layer keys are the epochs (*datetime.datetime*). The second layer keys are the PRN IDs of GPS satellites. The third layer keys are the XYZ coordinates of the position vectors, XYZ coordinates of the velocity vectors, the clock bias, and clock drift. Third layer values are all floats.

All GPS coordinates are in ITRF by default. No coordinate reference frame transformations were performed between the raw and the interpolated GPS ephemeris. 

.. automodule:: gpsxtr
   :members:

An example structure of this dictionary is given below::

	gpsdata = {epoch1 : {1 : {px:0,   py:0,   pz:0,
	                          vx:0,   vy:0,   vz:0,
	                          clkb:0, clkd:0      },
						 
	                     2 : {px:0,   py:0,   pz:0,
	                          vx:0,   vy:0,   vz:0,
	                          clkb:0, clkd:0      }, ...
							  
	                         ... ... ... ... ... ... ...
							 
	                    32 : {px:0,   py:0,   pz:0,    
	                          vx:0,   vy:0,   vz:0,    
	                          clkb:0, clkd:0      }},
							 
	           epoch2 : {1 : {px:0,   py:0,   pz:0,    
	                         ... ... ... ... ... ...} ...} ...}

The secondary output of **gpsxtr.gpsxtr()** is the **goodsats** list, which is a sorted list of GPS satellite PRN IDs whose observables had no outages.

An auxiliary function **gpsxtr.gpsweekday()** exists but is not documented here. It converts a datetime object into two string variables: the GPS day-of-week and the GPS week number. This is a similar function to **inpxtr.inptim()**.

rinxtr.py
---------

The fifth step involves the extraction of RINEX observations C1/P1, P2, L1, L2, D1, D2, with carrier phase marking, and cycle slip detection performed using the **phasep.py** module.

An additional 'flag' key will be added to the RINEX observables to mark them. These are the possible values belonging to the 'flag' keys.

- **"start" :** start of carrier phase observed, from an N-th GPS satellite, in a sequence.
- **"end" :** last carrier phase observed, from some N-th GPS satellite, in a sequence.
- **"solo" :** single carrier phase observation by some N-th GPS satellite, no sequence.
- **"none" :** carrier phase observation within a sequence, from some N-th GPS satellite.
- **"slip" :** cycle slip flag for carrier phase observed from some N-th GPS satellite.

An additional carrier phase term, L4, will also be added to the RINEX observables. This helps LEOGPS to perform cycle slip detection. If the user opts for single frequency processing, then L4 is the Melbourne-Wubbena linear combination. If the user opts for dual frequency processing, then L4 is the geometry-free linear combination.

If Doppler observables D1/D2 are not found in the RINEX observations, then Doppler values will be estimated through the **dopest.py** module, by estimating a first-order derivative of the L1/L2 phase values numerically using polynomial fitting.

.. automodule:: rinxtr
   :members:

The output RINEX observations are structured as a dictionary of epochs, each epoch with a sub-dictionary of GPS satellites based on IDs (1 to 32), and each ID with a sub dictionary of the various observations (C1/P1, P2, L1, L2, D1, D2). Example output structure::

	rnxproc = {epoch1 : {1 : {'C1':123, 'L1':123, ...
	                          'L4':321, 'flag':'none'} ...} ...
	           epoch2 : {2 : {'C1':123, 'L1':123, ...
	                          'L4':321, 'flag':'slip'} ...} ...
							 
			   ...
			  
	           epochX : {1 : {'C1':123, 'L1':123, ...
	                          'L4':321, 'flag':'none'} ...}}

If the user opts for code-carrier smoothing Hatch filter (either through the GUI or in the config file), then hatch filtering will be called in **rinxtr.rinxtr()** using the **phasep.py** module.

You may also change the length of the cycle slip filter, and the filter tolerance in terms of the number of standard deviations through the GUI or manually in the config file.

Do ensure that RINEX observation files follow 4-letter ID naming convention followed by the DOY + 0, with the file extension .YYO.

phasep.py
---------

The **phasep.py** module augments the RINEX data dictionary parsed out by **rinxtr.py** module. It comprises the cycle slip detection and marking function **phsmrk()**, as well as the hatch filtering algorithm **ph1fil()** for L1 observables, and **ph2fil()** for L1 + L2 observables. Within the hatch filtering loop, each code-phase data point at each time step is computed by the **hatch1()** or **hatch2()** functions.

.. automodule:: phasep
   :members:

.. note:: The RINEX data dictionary returned by **phasep.phsmrk(), phasep.ph1fil(), phasep.ph2fil()** all share the same nested dictionary key-value pairs as the output of the **rinxtr.rinxtr()** function.

dopest.py
---------

If Doppler observables D1/D2 are not found in the RINEX observations, then Doppler values will be estimated through the **dopest.py** module, by estimating a first-order derivative of the L1/L2 phase values numerically using polynomial fitting.

.. automodule:: dopest
   :members:

posvel.py
---------

The sixth step is to perform single point positioning (SPP), using the code pseudorange equations, solved via weighted least squares epoch-wise. Thus, the primary function **posvel.posvel()** is called once for each satellite and for each epoch.

In the code phase ranging equation, GPS satellite clock biases are offset from the output of the **gpsxtr.py** module. Ionospheric delays are handled too. In the L1 case, the GRAPHIC linear combination is used. In the L2 case, the ionosphere-free linear combination is used. Other effects such as the signal time-of-flight, the effects of Earth rotation, relativistic Shapiro effect, relativistic clock delays and clock advances are also offset. Functions to compute the relativistic effects are given in the next section, under **einstn.py**.

Doppler-based estimation of velocities will also be performed if Doppler data is available in the RINEX data dictionary parsed out by **rinxtr.py** module. By default, if Doppler data is missing, the **dopest.py** module would have worked its magic to estimate the Doppler values.

.. automodule:: posvel
   :members:

.. note:: Tropospheric effects are **not** handled in **posvel.py** as LEOGPS was built for spaceborne receivers and not for terrestrial receivers. However, if you wish to include terrestrial receivers, you can include the tropospheric modelled offsets (i.e. Saastamoinen, Hopfield, or Differential Refraction models etc) to the observed range variable *gpsrng_obsv* from lines 258 to 297.

In the main work flow **leorun.run()**, this function will be called once in each epoch for each of the two satellites. Both SPP results will be used in the carrier phase ambiguity estimation.

einstn.py
---------

The 'Einstein' module, comprises two main functions: one to compute the clock advance and one to compute the Shapiro path delay. In both functions, the output is converted to the equivalent path-length in meters.

.. automodule:: einstn
   :members:

.. note:: The rate of advance of two identical clocks, one in the LEO satellite and the other on the GPS satellite, will differ due to differences in the gravitational potential and to the relative speed between them.

.. note:: **Shapiro Delay:** Due to the space time curvature produced by the gravitational field, the Euclidean range travelled by the signal, which is computed by **posvel.py** must be corrected by the extra distance travelled. Typically, Shapiro effects corrupt the range model with about ~2cm ranging error.

This module will be called in **posvel.py** during the setup of pseudorange model in the iterative least squares solution of single-point positioning.

azimel.py
---------

While elevation-dependent or azimuth-dependent weighting is not built into the iterative least squares processing of LEOGPS for single point positioning, this module exists (but is not used) if the user wishes to compute azimuths and elevations from the LEO to GPS satellites anyway.

.. automodule:: azimel
   :members:

ambest.py
---------

This is the seventh processing step, which is the carrier phase integer or float ambiguity resolution step. This module contains functions that support epoch-wise processing for integer ambiguity resolution step.

The chief function in the module is the **'ambest()'** function, which outputs the precise relative baseline vector between the two spacecraft. This processing is done snapshot-wise, and thus has to be called for each epoch of carrier phase observations. All other supporting functions are called within **'ambest()'**. At the user-level, it is advised to modify contents only within 'ambest' unless the user wishes to modify the core float ambiguity resolution algorithm.

.. automodule:: ambest
   :members:

An option exists to perform integer fixing (see the `fix` argument above) using Peter Teunissen's LAMBDA method. A Pythonic translation of his original LAMBDA Integer-Least-Squares (ILS) Search-and-Shrink algorithm has been provided in the **ambfix.py** module. 

.. note:: If the user wishes to set their own custom zero difference covariance matrix, the user can input this in the optional `covZD` argument. If the user does not specify the `covZD` argument, then `covZD` by default will revert to an identity matrix scaled by the `sigma` argument above. Thus, running LAMBDA (by setting `fix = True` when calling **ambest.py**), but not setting a custom `covZD` argument, only has the equivalent effect of integer rounding.

ambfix.py
---------

This module holds the classical LAMBDA method that was originally authored by Teunissen, Jonge, and Tiberius (1993). The code was later written in MATLAB by Dr Sandra Verhagen and Dr Bofeng Li. It takes in a vector of float ambiguities to the integer least-squares (ILS) problem, and covariance of the float ambiguities. It then runs the LAMBDA's ILS search-&-shrink and spits out the ambiguity integers. The other 5 methods in original LAMBDA MATLAB code are not supported here (feel free to edit the code and implement it youself). The default `ncands = 2`, as per original code. All supporting functions from the original MATLAB code (`decorrel`, `ldldecom`, `ssearch`) have been nested within the main function as sub functions.

.. automodule:: ambfix
   :members:

.. note:: LAMBDA always first applies a decorrelation before the integer estimation. For ILS this is required to guarantee an efficient search. For rounding and bootstrapping it is required in order to get higher success rates (although rounding and bootstrapping is not included in LEOGPS).

frames.py
---------

This is the eighth step in the LEOGPS native processing work flow. This step performs the conversion of the coordinate reference frames between ITRF and ICRF, via the IAU1976 Theory of Precession and IAU1980 Theory of Nutation. For the visualisation of the formation geometry, it is recommended that the user select the Hill frame as the relative coordinate frame. By default, the reference frame in the downloaded ephemeris files in AIUB CODE's FTP is the ITRF.

.. automodule:: frames
   :members:

pubplt.py
---------

In the final stage, after all processing is done, the **pubplt.py** module publishes the information into output files in the `outputs` folder, found in the LEOGPS root directory.

Specifically, there are three functions in this module: a function to save as a plot graph the interpolated GPS ephemeris and clock biases; a function to save as a text report the interpolated GPS ephemeris and clock biases; and a function to save the final ephemeris and precise baselines estimated of both LEO-A and LEO-B.

.. automodule:: pubplt
   :members:

consts.py
---------

The following is a list of common constants used throughout LEOGPS, extracted from the University of Bern, Center for Orbit Determination in Europe (CODE)::

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

This API reference was automatically generated using Sphinx' Autodoc feature, using the `NumPy docstring format <https://numpydoc.readthedocs.io/en/latest/format.html>`_, and last updated on 6th June 2021.
