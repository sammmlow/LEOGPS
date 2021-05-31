#!/usr/bin/env python3

###############################################################################
###############################################################################
##                                                                           ##
##     _    ___  ___   ___ ___ ___                                           ##
##    | |  | __ /   \ / __| _ | __|                                          ##
##    | |__| __| ( ) | (_ |  _|__ \                                          ##
##    |____|___ \___/ \___|_| \___/                                          ##
##                                    v 1.2 (Stable)                         ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 30-May-2021.                                             ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##    Documentation: https://leogps.readthedocs.io/en/latest/                ##
##                                                                           ##
###############################################################################
###############################################################################

''' Module docstring for the main LEOGPS file

This is the main script that calls the GUI (tkinter) to run.

'''

# Import our GUI libraries.
import tkinter

# Import local libraries.
from source import leogui

# Initialise the GUI.
root = tkinter.Tk()
root_gui = leogui.run_gui( root )
root.mainloop()

# To run LEOGPS without GUI, comment out lines 27-29.
# Set your config.txt file manually.
# Uncomment below and run this script.
# from source import leorun
# leorun.run()