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

First Scenario Run
==================

Now, your LEOGPS scenario should have begun running.

First, let us discuss about the RINEX observation files. Notice that the RINEX observation files have file names that match the 4-letter ID and date-times in your LEOGPS scenario. For example, the GRACE satellites' decompressed observation file names are::

   GRCA2080.10O
   GRCB2080.10O

Notice also that in the **../LEOGPS/input/** directory, RINEX files ending with **D** in the file extension will automatically be Hatanaka-decompressed (file extension changes from 'D' to 'O'). LEOGPS accepts Hatanaka compressed and uncompressed RINEX files. The native LEOGPS v1.3 comes with four Hatanaka-compressed RINEX observation files for the GRACE A and B formation (a joint formation flying project by NASA and DLR), and the Lumelite A and B formation (by NUS STAR).

To use your own custom RINEX observation files for your own spacecraft, you should provide the RINEX files following the same file naming convention in the **/../LEOGPS/input/** folder, with RINEX version 2.XX only, one for each LEO satellite. The file naming convention for RINEX observation files is an 8-character string comprising the 4-letter ID of the spacecraft (you can create one arbitrarily), the 3-digit GPS day-of-the-year, followed by a single zero digit. The file extension is the last two digits of the year, followed by the letter 'O'. If the RINEX observation is Hatanaka-compressed, then the last letter of the file extension is a 'D'.

.. note:: For help in understanding the contents of RINEX formatted files better, LEOGPS comes packaged with the RINEX v2.10 documentation in the **reference** directory.

Second, for supporting files, notice also that LEOGPS will automatically download daily precise (final) ephemeris and clock files from the `University of Bern's CODE FTP <ftp://ftp.aiub.unibe.ch/CODE/>`_, and unzip them automatically. Ephemeris files downloaded will span one day before and one day after the scenario date(s) so that the GPS satellite orbit interpolation can be done beyond the edges of the scenario time-line to prevent edge effects (i.e. ringing, poor polynomial fits etc).

As long as the user has an active internet connection, you do not need to provide your own ephemeris and clock files, as LEOGPS will automatically source it from the University of Bern, CODE. However, in the event where the FTP is down, but the user has their own precise ephemeris and clock files, you can transfer these files into the **/../input/** folder and rename these files to follow the CODE naming convention so that LEOGPS detects these files and can use them offline without access to the FTP.

The naming convention for precise clock and ephemeris files are also 8-character strings, starting with a 'COD', followed by the 4-digit GPS week number, followed by the 1-digit GPS day-of-week number. The file extension for precise ephemeris files is a '.EPH' string, and for clock files it is a '.CLK' string.

For example, an ephemeris (.EPH) and clock file (.CLK) from CODE, on the second day of the Week (Tuesday), and on GPS Week Number 1594, would be named::

   COD15942.EPH
   COD15942.CLK

.. note:: For a reference on geodetic calendar formats, refer to the `this calendar <https://geodesy.noaa.gov/CORS/Gpscal.shtml>`_.

Third, the interpolated data of GPS ephemeris and clocks are saved and plotted in the **../LEOGPS/output/** folders. There would be a text file report of interpolated GPS ephemeris, as well as graphical plots of GPS satellites of PRN IDs 1 to 32 (or fewer, due to spacecraft outages, manoeuvres, servicing etc).

Fourth and finally, the desired main output of LEOGPS - the single point positioning (SPP) and carrier phase differential GPS (CDGPS) solutions, are saved in **'LEOGPS_Results.txt'**. You can open this text file using Excel or any string formatting language you would like to use. Note that for the individual spacecraft ephemeris, the outputs are in the coordinate reference frame chosen in the GUI under the input label **Orbit frame**; for the relative baseline vector (relative position vector) the outputs are in the coordinate reference frame chosen in the GUI under the input label **Relative orbit frame**. It is recommended (and of more use) to display the relative orbit in the Euler-Hill reference frame, which takes Spacecraft A as the chief or reference spacecraft.

In the native LEOGPS build, the ground truths for the GRACE formation are also provided in the **reference** directory. These truths are the precise orbit determination solutions, with relative position magnitudes validated by the GRACE formation's actual K-Band ranging radar data.

This concludes the explanation of the basic setup of LEOGPS. For more information on relative baseline determination, refer to the tutorials on differential GPS in the next page(s).
