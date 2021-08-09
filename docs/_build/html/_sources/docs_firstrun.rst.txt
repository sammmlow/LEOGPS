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

First Scenario Run
==================

First, for your own custom scenarios, set the correct 4-letter ID and date-times of your spacecraft scenario. When saving or running LEOGPS, these inputs for your scenario and spacecraft will be saved in the **../LEOGPS/config/config.txt/** file which the GUI interfaces with. Note that configuration input strings follow a strict format.

Second, to use your own custom RINEX observation files, you should provide this in the **/../input/** folder, with two RINEX (v2.xx) observation files, one for each LEO satellite, by pasting it in the **../LEOGPS/input/** folder.

The file naming convention for the RINEX observation files is an 8-character string comprising the 4-letter ID of the spacecraft (you can create one arbitrarily), the 3-digit GPS day-of-the-year, followed by a single zero digit. The file extension is the last two digits of the year, followed by the letter 'O'. If the RINEX observation is Hatanaka-compressed, then the last letter of the file extension is a 'D'.

For example, the GRACE satellites' decompressed observation file names are::

   GRCA2080.10O
   GRCB2080.10O

The 4-letter ID of the two spacecrafts in the LEOGPS GUI must match the first 4-letters of the RINEX files provided!

.. note:: For help in understanding the RINEX format better, LEOGPS comes packaged with the RINEX v2.10 documentation in the **reference** directory.

For precise ephemeris and clock files, as long as the user has an active internet connection, you do not need to provide your own ephemeris and clock files, as LEOGPS will automatically source it from the University of Bern, CODE. However, in the event where the FTP is down, but the user has their own precise ephemeris and clock files, you can transfer these files into the **/../input/** folder and rename these files to follow the CODE naming convention so that LEOGPS detects these files and can use them offline without access to the FTP.

The naming convention for precise clock and ephemeris files are also 8-character strings, starting with a 'COD', followed by the 4-digit GPS week number, followed by the 1-digit GPS day-of-week number. The file extension for precise ephemeris files is a '.EPH' string, and for clock files it is a '.CLK' string.

For example, an ephemeris (.EPH) and clock file (.CLK) from CODE, on the second day of the Week (Tuesday), and on GPS Week Number 1594, would be named::

   COD15942.EPH
   COD15942.CLK

.. note:: For a reference on geodetic calendar formats, refer to the `this calendar <https://geodesy.noaa.gov/CORS/Gpscal.shtml>`_.

This concludes the basic setup of the LEOGPS default build.

