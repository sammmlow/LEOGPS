![LEOGPS - GPS-Aided Relative Satellite Navigation in Python](https://raw.githubusercontent.com/sammmlow/LEOGPS/master/docs/LEOGPS_Logo.png)

LEOGPS is an open-source Python package for absolute and relative navigation.

Absolute positioning is performed by trilaterating GPS pseudorange measurements, and Doppler (pseudorange rate) measurements. If Doppler values are missing in the RINEX observation file, LEOGPS will attempt to estimate them. The relative navigation between LEO satellites are performed using a double-differencing of carrier phase values, with a simple rounding of the float ambiguities.

To use LEOGPS, the user first inputs configuration parameters in the 'config.txt' file. Then, the user pastes two RINEX (v2.xx) observation files, one for each LEO satellite. Next, the user simply has to run 'leogps.py' with an internet connection online. LEOGPS will process the raw GPS measurements to produce a report comprising:

- The absolute positions and absolute velocities of both LEOs.
- Precise (centimeter-level) baseline estimation (relative position vector).
- Dilution of precision values.

LEOGPS is useful for formation flying satellite missions such as GRACE A and B, and can also be adapted for the rapid prototyping of navigation algorithms in Python, or for testing out integer ambiguity resolution techniques. A Pythonic translation of Professor Peter Teunissen's LAMBDA method has also been adapted, in the 'file ambfix.py'

To install, run, and develop your own programs with LEOGPS, please see the detailed documentation below:
<Insert URL here>
