![LEOGPS - GPS-Aided Relative Satellite Navigation in Python](https://raw.githubusercontent.com/sammmlow/LEOGPS/master/gui/logo.png)

### What is LEOGPS?
LEOGPS is an open-source Python package for absolute and relative satellite navigation.

### How do I use LEOGPS?

First, clone this repository. Then, that's it! The user first runs the application by running 'leogps.py', the main python file. The app interfaces with the inputs configuration file, allowing the user to save or load parameters in the 'config.txt' file, through the GUI.

Then, the user simply provides the program with two RINEX (v2.xx) observation files, one for each LEO satellite, by pasting it in the 'inputs' folder. By default, LEOGPS should be provided with two example RINEX files comprising GPS pseudorange and carrier phase measurements from two satellites in formation - GRACE A and B. Ground truths are also provided in the 'outputs' directory.

Next, the user simply has to run LEOGPS via the GUI, once the user saves the desired configuration options. Keep an internet connection online, as LEOGPS will attempt to download precise ephemerides and clock bias data from the University of Bern's COD FTP (unless you already have those files). Note that in the previous versions, LEOGPS used International GNSS Service (IGS) products. We have now switched to Bern's COD in v0.3.

The app should look something like the image below.

![LEOGPS - Graphical User Interface](https://raw.githubusercontent.com/sammmlow/LEOGPS/master/gui/gui.jpg)

LEOGPS will process the raw GPS measurements to produce a report comprising:

- The absolute positions and absolute velocities of both LEOs.
- Precise (centimeter-level) baseline estimation (relative position vector).
- Dilution of precision values.

LEOGPS will also (optionally, depending on your choice in the GUI) output plots or reports on the interpolated GPS satellite ephemeris and clock biases. Currently ephemeris errors are still present and we are working to resolve them.

### What is the method of navigation used behind LEOGPS?

Absolute positioning and velocity estimation is performed by trilaterating GPS pseudorange measurements, and Doppler (pseudorange rate) measurements respectively. Single-frequency measurements employ the GRAPHIC linear combination (T. P. Yunck, 1993), while dual frequency measurements employ the ionosphere-free linear combination (Misra & Enge, 2006). If Doppler values are missing in the RINEX observation file, LEOGPS will also attempt to estimate them using a gradient estimation via first order Taylor expansion approximation of the carrier phase values (in v0.2 onwards). However, the velocity estimation results so far are not promising.

The relative navigation between LEO satellites are performed using a double-differencing of carrier phase values, and using the float ambiguities directly. An optional module exists ('ambfix.py') in the codes folder for Dr Peter Teunissen's LAMBDA method for estimating integer ambiguities given the float ambiguities and covariance matrices.

### Why did you create LEOGPS?

LEOGPS is a personal pet project of mine, useful as a modular software used to study relative satellite navigation techniques for distributed space systems with a GPS receiver, such as formation missions like GRACE A and B. LEOGPS' sole purpose is as an adaptable platform for the rapid prototyping of navigation algorithms in Python, or for testing out integer ambiguity resolution techniques. A Pythonic translation of Professor Peter Teunissen's LAMBDA method has also been adapted, in the 'file ambfix.py'

### Requirements for LEOGPS:

Requires Python version 3.6.5 (Anaconda with Spyder recommended).

Core libraries necessary: NumPy (v1.14 and above) and matplotlib

Standard Python libaries: os, copy, math, datetime, decimal, shutil, subprocess, warnings, urllib.request

Libraries for GUI: PIL, tkinter

The above are all standard Python libraries that should come with a typical Anaconda installation.

This README was last updated on: 12-Jan-2021.

### Contact:

If you have any queries feel free to reach out to me at:
sammmlow@gmail.com
