..
   ###########################################################################
   ###########################################################################
   ##                                                                       ##
   ##     _    ___  ___   ___ ___ ___                                       ##
   ##    | |  | __ /   \ / __| _ | __|                                      ##
   ##    | |__| __  ( ) | (_ |  _|__ \                                      ##
   ##    |____|___ \___/ \___|_| \___/                                      ##
   ##                                    v 1.3 (Stable)                     ##
   ##                                                                       ##
   ###########################################################################
   ###########################################################################

.. image:: /_static/leogps_logo.png

|

First Steps
===========

First, launch LEOGPS by running **'leogps.py'** in a Python IDE or your terminal, and you will see a user-interface (UI) with input parameters as shown below.

.. figure:: /_figures/gui-v1-3.jpg
   :align: center
   
   Graphical User Interface (GUI) of LEOGPS

On the left face of the GUI, you would find text entries where you would input your scenario and spacecraft parameters. On the right face of the GUI, you would find an embedded plot object (matplotlib) which would display the relative orbit trajectory after resolving for the relative positions (you can select the frame, but the default should always be the local Hill Frame with satellite LEO-A as the chief or origin).

Let us walk through the individual input entries and options on the left side of the GUI. This will be where LEOGPS reads and loads the `config.txt` file into the user interface.

####

**01. Spacecraft A 4-letter ID (LEOA)**

LEOGPS searches for the Spacecraft A's RINEX file using the spacecraft 4-letter ID. You can create one arbitrarily (letters only), just make sure it matches the first 4-letters of the file. For example, the 4-letter ID of GRACE A was set to GRCA in the file name, followed by the 3-digit GPS day-of-the-year, a 1-digit zero, and a file extension of the 2-digit year plus an 'O' for an observation file, or 'D' for a Hatanaka-compressed observation file (GRCA2080.10D for example).

####

**02. Spacecraft B 4-letter ID (LEOB)**

LEOGPS searches for the Spacecraft B's RINEX file using the spacecraft 4-letter ID. You can create one arbitrarily (letters only), just make sure it matches the first 4-letters of the file. For example, the 4-letter ID of GRACE B was set to GRCB in the file name, followed by the 3-digit GPS day-of-the-year, a 1-digit zero, and a file extension of the 2-digit year plus an 'O' for an observation file, or 'D' for a Hatanaka-compressed observation file (GRCB2080.10D for example).

####

**03. Epoch start (YYYY-MM-DD-HH-MN-SS)**

Specify the starting epoch for processing in LEOGPS (not necessarily the epoch where the first RINEX entry is found). Follow the date-time format in the heading above, starting with YYYY. Note that the two RINEX files must have an overlapping time frame in order to compute valid differential GPS results. Note also that the time scale is in GPST and not UTC!

####

**04. Epoch final (YYYY-MM-DD-HH-MN-SS)**

Specify the ending epoch for processing in LEOGPS (not necessarily the epoch where the last RINEX entry is found). Follow the date-time format in the heading above, starting with YYYY. Note that the two RINEX files must have an overlapping time frame in order to compute valid differential GPS results. Note also that the time scale is in GPST and not UTC!

####

**05. Timestep in seconds (i.e. 30)**

Time step used in processing of RINEX data (integer seconds only).

####

**06. Enable L1 or L1/L2 processing?**

LEOGPS can read L1 or L2 signals from GPS observations in the RINEX file. The L1 signal is the oldest legacy GPS signal, comprising two parts: the Coarse/Acquisition Code (C1) and the Precision Code (P1). The P-code is reserved for military use, while the C/A is open to the public. The L1 signal uses the frequency 1575.42 MHz. The L2 uses the frequency 1227.60 MHz, and was implemented after the L1 and must be used along with L1 frequencies. It also has a military code and a civilian use code. Enabling both L1/L2 allows LEOGPS to perform single point positioning computations using the ionosphere-free linear combination, where a differential signal between L1 and L2 is used to cancel out correlated ionospheric delay errors.

####

**07. Enable code-carrier Hatch filter?**

LEOGPS allows for noisy (but unambiguous) code measurements to be smoothed using the precise (but ambiguous) carrier phase measurements. Such an algorithm combines measurements of both from the RINEX file and the filter that implements the algorithm is known as a Hatch filter. Note that enabling Hatch filtering significantly increases the processing time and memory usage.

####

**08. Set Hatch filter window length**

LEOGPS implements the Hatch filter as a sliding window filter. Thus, you can set the length of the sliding window here as an integer number of samples. For example, if your time step selected was 10 seconds, and the window length input in this option was set to 30, then your effective sliding window duration is 300 seconds.

####

**09. Number of sigmas in cycle slip detection**

LEOGPS implements a very basic carrier phase cycle slip detection algorithm by performing an interpolation of combined carrier phase data, and observing if there are any single points of data that exceed "X" sigmas of the interpolated carrier phase. "X" is the input specified in this option here. For L1-only processing, LEOGPS uses the Melbourne-Wubbena linear combination as the carrier phase combination. For L1/L2 processing, LEOGPS uses the geometry-free linear combination as the carrier phase combination. For carrier phase observations that exceed the "X" number of sigmas, LEOGPS will mark that point with a cycle slip flag print a warning message to the user in the terminal. However, as of Version 1.3, LEOGPS does not attempt to repair or reject that particular observation.

####

**10. Set cycle slip detection filter length**

LEOGPS performs interpolation across a sliding window of samples in the cycle slip detection process. The user specifies how many samples N are used for carrier phase interpolation and screening. The larger the value N, the longer the processing but the lesser the likelihood of flagging a false cycle slip.

####

**11. Set the GPS antenna X offset (DISABLED)**

As the GPS antenna is unlikely to be positioned exactly in the center of mass of the spacecraft, the navigation solutions likely center about the GPS antenna phase center rather than the body center of mass. This option allows the user to specify the body-frame offset. **Note: this option is disabled for now as LEOGPS has not been configured to read spacecraft attitude files, needed for orbit-to-body-frame transformations.**

####

**12. Set the GPS antenna Y offset (DISABLED)**

As the GPS antenna is unlikely to be positioned exactly in the center of mass of the spacecraft, the navigation solutions likely center about the GPS antenna phase center rather than the body center of mass. This option allows the user to specify the body-frame offset. **Note: this option is disabled for now as LEOGPS has not been configured to read spacecraft attitude files, needed for orbit-to-body-frame transformations.**

####

**13. Set the GPS antenna Z offset (DISABLED)**

As the GPS antenna is unlikely to be positioned exactly in the center of mass of the spacecraft, the navigation solutions likely center about the GPS antenna phase center rather than the body center of mass. This option allows the user to specify the body-frame offset. **Note: this option is disabled for now as LEOGPS has not been configured to read spacecraft attitude files, needed for orbit-to-body-frame transformations.**

####

**14. Orbit frame (not shown in the plot)**

LEOGPS will output a text file reporting the single point positions and velocities of the individual spacecraft, as well as the relative baseline vectors. This option allows the user to toggle the output reference frame for the single point positions and velocities. The ITRF-ICRF (ECEF-ECI) conversion uses the IAU1976 Theory of Precession and IAU1980 Theory of Nutation in the celestial-to-terrestrial (and vice versa) conversion. Note that the plotter displays only the relative baselines and not the single point positions.

####

**15. Relative orbit frame (default Hill)**

LEOGPS will plot in the GUI (the big plot box on the right) as well as output in a text file the relative baseline vectors. In practice, the Euler-Hill frame is typically used (and is the default selection) but the user can also toggle the relative baseline data and plots to be displayed in ITRF or ICRF too.

####

We may now proceed to run the default LEOGPS scenario that comes with the build, for the GRACE formation flying satellite mission. On the UI, click **Load Config**. This loads the configuration inputs from **config.txt** into the blank entries explained ebove. If you did not make any changes to **config.txt**, the inputs should correspond to the GRACE formation flying scenario on 27-07-2010, as the default example. 

Whenever saving (**Save Config**) or running LEOGPS (**Run LEOGPS**), all inputs are saved in the **../LEOGPS/config/config.txt/** file.

Now, we are ready to run LEOGPS. Once we do, several things will happen (next page).
