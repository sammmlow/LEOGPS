![LEOGPS - GPS-Aided Relative Satellite Navigation in Python](https://raw.githubusercontent.com/sammmlow/LEOGPS/master/gui/logo.png)

### LEOGPS
LEOGPS is an open-source Python package for absolute and relative satellite navigation, with specific intent for relative navigation between formation flying satellites. It currently supports only observations from the GPS constellation, up to L1/L2 frequency, with input files as RINEX v2.XX. It uses ephemeris and clock files from the University of Bern, CODE.

### How do I use LEOGPS?

First, for the full documentation, please refer to it here:

First, clone this repository. Then, that's it! The user first runs the application by running 'leogps.py', the main python file. The app interfaces with your inputs, and saves it in a configuration file, allowing the user to save or load parameters through the GUI.

Then, the user simply provides the program with two RINEX (v2.xx) observation files, one for each LEO satellite, by pasting it in the 'inputs' folder. By default, LEOGPS should be provided with two example RINEX files comprising GPS pseudorange and carrier phase measurements from two satellites in formation - GRACE A and B. Ground truths are also provided in the 'outputs' directory.

Next, the user can run LEOGPS. The UI of the app should pop out and look something like below:

![LEOGPS - Graphical User Interface](https://raw.githubusercontent.com/sammmlow/LEOGPS/master/gui/gui_v1.jpg)

Keep an internet connection online, as LEOGPS will attempt to download precise ephemerides and clock bias data from the University of Bern's COD FTP (unless you already have those files). Note that in the previous versions, LEOGPS used International GNSS Service (IGS) products. We have now switched to Bern's COD in v0.3 onwards.

LEOGPS will process the raw GPS measurements to produce a report comprising:

- The single point positions and velocities of both LEOs.
- Precise baseline vectors between the two LEOs.
- Dilution of precision values.

LEOGPS will also (optionally, depending on your choice in the GUI) output plots or reports on the interpolated GPS satellite ephemeris and clock biases. Currently ephemeris errors are still present and we are working to resolve them.

### Requirements for LEOGPS:

Tested on Python version 3.6.5 (Anaconda with Spyder).

Core libraries necessary: NumPy (v1.14 and above) and matplotlib

Standard Python libaries: os, copy, math, datetime, decimal, shutil, subprocess, warnings, urllib.request

Libraries for GUI: PIL, tkinter

The above are all standard Python libraries that should come with a typical Anaconda installation.

This README was last updated on: 29-March-2021.

### Contact:

If you have any queries feel free to reach out to me at:
sammmlow@gmail.com
