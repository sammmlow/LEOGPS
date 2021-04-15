#!/usr/bin/env python3
'''
###############################################################################
###############################################################################
##                                                                           ##
##     _    ___  ___   ___ ___ ___                                           ##
##    | |  | __ /   \ / __| _ | __|                                          ##
##    | |__| __  ( ) | (_ |  _|__ \                                          ##
##    |____|___ \___/ \___|_| \___/                                          ##
##                                    v 1.1 (Stable)                         ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 15-Apr-2021.                                             ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##                                                                           ##
###############################################################################
###############################################################################
'''

# Import our GUI libraries.
import tkinter

# Import local libraries.
from codes import leogui

# Initialise the GUI.
root = tkinter.Tk()
root_gui = leogui.run_gui( root )
root.mainloop()

# To run LEOGPS without GUI, comment out lines 27-29.
# Set your config.txt file manually.
# Uncomment below and run this script.
# from codes import leorun
# leorun.run()