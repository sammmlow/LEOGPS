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
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 18-Jul-2022 (By Andrew G. T. Ng).                                             ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##    Documentation: https://leogps.readthedocs.io/en/latest/                ##
##                                                                           ##
###############################################################################
###############################################################################

# Import global libraries
import tkinter

# Import local libraries.
from source import leogui

# Initialise the GUI.
root = tkinter.Tk()
root_gui = leogui.RunGUI( root )
root.mainloop()

# To run LEOGPS manually without GUI, comment out lines 27-29. Then, edit 
# your config.txt file manually. Finally, uncomment lines 34-35 and run.

# from source import leorun
# leorun.run()