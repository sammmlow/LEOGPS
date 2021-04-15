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
##    FILE DESCRIPTION:                                                      ##
##                                                                           ##
##    This file contains the GUI class which based on Python tkinter.        ##
##    The class will be called in the main LEOGPS python file.               ##
##    (No inputs and outputs, this file only holds the GUI class object.     ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 13-03-2020.                                              ##
##                                                                           ##
###############################################################################
###############################################################################
'''

# IMPORT PUBLIC LIBRARIES
import tkinter as tk
from PIL import Image, ImageTk
from os.path import dirname, abspath

# IMPORT THE LOCAL LEOGPS LIBRARY
from codes import leorun

class run_gui:

    def __init__(self, master):
        
        # Create the main frame and window.
        master.title('LEOGPS v1.0')
        master.geometry('1024x896')
        
        # Initialise the basic text labels (found in the configuration file):
        self.txt1  = 'Input the 4-letter ID of the 1st spacecraft (i.e. LEOA).'
        self.txt2  = 'Input the 4-letter ID of the 2nd spacecraft (i.e. LEOB).'
        self.txt3  = 'Input the starting epoch in YYYY-MM-DD-HH-MN-SS.'
        self.txt4  = 'Input the ending epoch in YYYY-MM-DD-HH-MN-SS.'
        self.txt5  = 'Input the timestep in seconds (i.e. 30).'
        self.txt6  = 'Input single or dual frequency processing (i.e. 1 or 2).'
        self.txt7  = 'Offset Earths rotation during signal time-of-flight?'
        self.txt8  = 'Offset ephemeris errors due to relativistic effects?'
        self.txt9  = 'Save graph plots for GPS ephemeris and clocks?'
        self.txt10 = 'Save the report for GPS ephemeris and clocks?'
        self.txt11 = 'Code-carrier smoothing via Hatch filtering? (True/False)'
        self.txt12 = 'Window length (no.of observations) of the hatch filter?'
        self.txt13 = 'Standard deviation tolerance for cycle slip detection?'
        self.txt14 = 'Length of sliding window cycle slip detection filter?'
        self.txt15 = 'Set the GPS antenna X offset (m) for vehicle.'
        self.txt16 = 'Set the GPS antenna Y offset (m) for vehicle.'
        self.txt17 = 'Set the GPS antenna Z offset (m) for vehicle.'        
        
        # Initialise tkinter variables for the entries corresponding to above.
        self.var01 = tk.StringVar() # 4-letter ID of LEO A
        self.var02 = tk.StringVar() # 4-letter ID of LEO B
        self.var03 = tk.StringVar() # Start epoch datetime
        self.var04 = tk.StringVar() # Ending epoch datetime
        self.var05 = tk.IntVar()    # Time step per epoch
        self.var06 = tk.IntVar()    # Frequency number
        self.var07 = tk.IntVar()    # Earth rotation offset flag
        self.var08 = tk.IntVar()    # Relativistic effects flag
        self.var09 = tk.IntVar()    # Save flag for GPS plots
        self.var10 = tk.IntVar()    # Save flag for GPS report
        self.var11 = tk.IntVar()    # Enable flag for hatch filter
        self.var12 = tk.IntVar()    # Hatch filter window length
        self.var13 = tk.DoubleVar() # Cycle slip tolerance STDs
        self.var14 = tk.IntVar()    # Cycle slip detection window
        self.var15 = tk.DoubleVar() # Antenna offset X
        self.var16 = tk.DoubleVar() # Antenna offset Y
        self.var17 = tk.DoubleVar() # Antenna offset Z
        
        # Define the path to the LEOGPS logo file.
        leogps_logo = dirname(dirname(abspath(__file__)))
        leogps_logo = leogps_logo + '\gui\logo.jpg'
        
        # Configure the background image and load the logo.
        image = Image.open( leogps_logo )
        photo = ImageTk.PhotoImage(image)
        self.logo = tk.Label(image=photo)
        self.logo.image = photo
        self.logo.grid(row=0, column=0, padx=20, pady=20)
        
        # Add a button to read default entries from 'config.txt'.
        self.cfgR = tk.Button(master, text='Load Config', command=self.cfg_R)
        self.cfgR.grid(row=0, column=1, padx=20, pady=2)
        self.cfgR.configure(bg="light blue")
        
        # Add a button to save entries into 'config.txt'.
        self.cfgW = tk.Button(master, text='Save Config', command=self.cfg_W)
        self.cfgW.grid(row=0, column=2, padx=20, pady=2)
        self.cfgW.configure(bg="light blue")
        
        # Add a button to run LEOGPS.
        self.cfgW = tk.Button(master, text='Run LEOGPS', command=self.run)
        self.cfgW.grid(row=0, column=3, padx=20, pady=2)
        self.cfgW.configure(bg="light blue")
        
        # Input the 4-letter ID of the first spacecraft (i.e. LEOA).
        self.label01 = tk.Label(master, text=self.txt1 )
        self.label01.grid(row=1, column=0, padx=20, pady=2, sticky='w')
        self.entry01 = tk.Entry(master, width=5, textvariable=self.var01)
        self.entry01.grid(row=1, column=1, padx=5, pady=2, sticky='w')
        self.errtx01 = tk.Label(master, text='', fg='red' )
        self.errtx01.grid(row=1, column=3, padx=5, pady=2, sticky='w')
        
        # Input the 4-letter ID of the second spacecraft (i.e. LEOB).
        self.label02 = tk.Label(master, text=self.txt2 )
        self.label02.grid(row=2, column=0, padx=20, pady=2, sticky='w')
        self.entry02 = tk.Entry(master, width=5, textvariable=self.var02)
        self.entry02.grid(row=2, column=1, padx=5, pady=2, sticky='w')
        self.errtx02 = tk.Label(master, text='', fg='red' )
        self.errtx02.grid(row=2, column=3, padx=5, pady=2, sticky='w')
        
        # Input the starting epoch in YYYY-MM-DD-HH-MN-SS.
        self.label03 = tk.Label(master, text=self.txt3 )
        self.label03.grid(row=3, column=0, padx=20, pady=2, sticky='w')
        self.entry03 = tk.Entry(master, width=20, textvariable=self.var03)
        self.entry03.grid(row=3, column=1, padx=5, pady=2, sticky='w',
                          columnspan=2)
        self.errtx03 = tk.Label(master, text='', fg='red' )
        self.errtx03.grid(row=3, column=3, padx=5, pady=2, sticky='w')
        
        # Input the ending epoch in YYYY-MM-DD-HH-MN-SS.
        self.label04 = tk.Label(master, text=self.txt4 )
        self.label04.grid(row=4, column=0, padx=20, pady=2, sticky='w')
        self.entry04 = tk.Entry(master, width=20, textvariable=self.var04)
        self.entry04.grid(row=4, column=1, padx=5, pady=2, sticky='w',
                          columnspan=2)
        self.errtx04 = tk.Label(master, text='', fg='red' )
        self.errtx04.grid(row=4, column=3, padx=5, pady=2, sticky='w')
        
        # Input the timestep in seconds (i.e. 30).
        self.label05 = tk.Label(master, text=self.txt5 )
        self.label05.grid(row=5, column=0, padx=20, pady=2, sticky='w')
        self.entry05 = tk.Entry(master, width=4, textvariable=self.var05)
        self.entry05.grid(row=5, column=1, padx=5, pady=2, sticky='w')
        self.errtx05 = tk.Label(master, text='', fg='red' )
        self.errtx05.grid(row=5, column=3, padx=5, pady=2, sticky='w')
        
        # Input single or dual frequency processing (i.e. 1 or 2).
        self.label06 = tk.Label(master, text=self.txt6 )
        self.label06.grid(row=6, column=0, padx=20, pady=2, sticky='w')
        self.entry06a = tk.Radiobutton(master, text='L1',
                                       variable=self.var06, value=1)
        self.entry06a.grid(row=6, column=1, padx=5, pady=2, sticky='w')
        self.entry06b = tk.Radiobutton(master, text='L1+2',
                                       variable=self.var06, value=2)
        self.entry06b.grid(row=6, column=2, padx=5, pady=2, sticky='w')
        self.errtx06 = tk.Label(master, text='', fg='red' )
        self.errtx06.grid(row=6, column=3, padx=5, pady=2, sticky='w')
        
        # Offset Earths rotation during signal TOF? (True/False)
        self.label07 = tk.Label(master, text=self.txt7 )
        self.label07.grid(row=7, column=0, padx=20, pady=2, sticky='w')
        self.entry07 = tk.Checkbutton(master, text='True/False',
                                      variable=self.var07)
        self.entry07.grid(row=7, column=1, padx=5, pady=2, sticky='w')
        self.errtx07 = tk.Label(master, text='', fg='red' )
        self.errtx07.grid(row=7, column=3, padx=5, pady=2, sticky='w')
        
        # Account for relativistic effects? (True/False)
        self.label08 = tk.Label(master, text=self.txt8 )
        self.label08.grid(row=8, column=0, padx=20, pady=2, sticky='w')
        self.entry08 = tk.Checkbutton(master, text='True/False',
                                      variable=self.var08)
        self.entry08.grid(row=8, column=1, padx=5, pady=2, sticky='w')
        self.errtx08 = tk.Label(master, text='', fg='red' )
        self.errtx08.grid(row=8, column=3, padx=5, pady=2, sticky='w')
        
        # Save plots for GPS ephemeris and clocks? (True/False)
        self.label09 = tk.Label(master, text=self.txt9 )
        self.label09.grid(row=9, column=0, padx=20, pady=2, sticky='w')
        self.entry09 = tk.Checkbutton(master, text='True/False',
                                      variable=self.var09)
        self.entry09.grid(row=9, column=1, padx=5, pady=2, sticky='w')
        self.errtx09 = tk.Label(master, text='', fg='red' )
        self.errtx09.grid(row=9, column=3, padx=5, pady=2, sticky='w')
        
        # Save the report for GPS ephemeris and clocks? (True/False)
        self.label10 = tk.Label(master, text=self.txt10 )
        self.label10.grid(row=10, column=0, padx=20, pady=2, sticky='w')
        self.entry10 = tk.Checkbutton(master, text='True/False',
                                      variable=self.var10)
        self.entry10.grid(row=10, column=1, padx=5, pady=2, sticky='w')
        self.errtx10 = tk.Label(master, text='', fg='red' )
        self.errtx10.grid(row=10, column=3, padx=5, pady=2, sticky='w')
        
        # Code-carrier smoothing via Hatch filtering? (True/False)
        self.label11 = tk.Label(master, text=self.txt11 )
        self.label11.grid(row=11, column=0, padx=20, pady=2, sticky='w')
        self.entry11 = tk.Checkbutton(master, text='True/False',
                                      variable=self.var11)
        self.entry11.grid(row=11, column=1, padx=5, pady=2, sticky='w')
        self.errtx11 = tk.Label(master, text='', fg='red' )
        self.errtx11.grid(row=11, column=3, padx=5, pady=2, sticky='w')
        
        # Window length (number of observations) of the hatch filter?
        self.label12 = tk.Label(master, text=self.txt12 )
        self.label12.grid(row=12, column=0, padx=20, pady=2, sticky='w')
        self.entry12 = tk.Entry(master, width=4, textvariable=self.var12)
        self.entry12.grid(row=12, column=1, padx=5, pady=2, sticky='w')
        self.errtx12 = tk.Label(master, text='', fg='red' )
        self.errtx12.grid(row=12, column=3, padx=5, pady=2, sticky='w')
        
        # Standard deviation tolerance for cycle slip detection?
        self.label13 = tk.Label(master, text=self.txt13 )
        self.label13.grid(row=13, column=0, padx=20, pady=2, sticky='w')
        self.entry13 = tk.Entry(master, width=4, textvariable=self.var13)
        self.entry13.grid(row=13, column=1, padx=5, pady=2, sticky='w')
        self.errtx13 = tk.Label(master, text='', fg='red' )
        self.errtx13.grid(row=13, column=3, padx=5, pady=2, sticky='w')
        
        # Length of the sliding window cycle slip detection filter?
        self.label14 = tk.Label(master, text=self.txt14 )
        self.label14.grid(row=14, column=0, padx=20, pady=2, sticky='w')
        self.entry14 = tk.Entry(master, width=4, textvariable=self.var14)
        self.entry14.grid(row=14, column=1, padx=5, pady=2, sticky='w')
        self.errtx14 = tk.Label(master, text='', fg='red' )
        self.errtx14.grid(row=14, column=3, padx=5, pady=2, sticky='w')
        
        # Set the GPS antenna offset in X-direction (m) for vehicle.
        self.label15 = tk.Label(master, text=self.txt15 )
        self.label15.grid(row=15, column=0, padx=20, pady=2, sticky='w')
        self.entry15 = tk.Entry(master, width=10, textvariable=self.var15)
        self.entry15.grid(row=15, column=1, padx=5, pady=2, sticky='w')
        self.errtx15 = tk.Label(master, text='', fg='red' )
        self.errtx15.grid(row=15, column=3, padx=5, pady=2, sticky='w')
        
        # Set the GPS antenna offset in Y-direction (m) for vehicle.
        self.label16 = tk.Label(master, text=self.txt16 )
        self.label16.grid(row=16, column=0, padx=20, pady=2, sticky='w')
        self.entry16 = tk.Entry(master, width=10, textvariable=self.var16)
        self.entry16.grid(row=16, column=1, padx=5, pady=2, sticky='w')
        self.errtx16 = tk.Label(master, text='', fg='red' )
        self.errtx16.grid(row=16, column=3, padx=5, pady=2, sticky='w')
        
        # Set the GPS antenna offset in Z-direction (m) for vehicle.
        self.label17 = tk.Label(master, text=self.txt17 )
        self.label17.grid(row=17, column=0, padx=20, pady=2, sticky='w')
        self.entry17 = tk.Entry(master, width=10, textvariable=self.var17)
        self.entry17.grid(row=17, column=1, padx=5, pady=2, sticky='w')
        self.errtx17 = tk.Label(master, text='', fg='red' )
        self.errtx17.grid(row=17, column=3, padx=5, pady=2, sticky='w')
    
    # Method to load default values from the configuration file.
    def cfg_R(self):
        
        cwd = dirname(dirname(abspath(__file__))) # Current working directory
        iwd = cwd + '\config\config.txt' # Inputs files
        inputfile = open(iwd,'r') # Open the config.txt file
        inps = {} # Create a dictionary to store all the input 
        integers = [ 'freq','timestep','hatchlength','cycsliplen' ]
        floats = ['cycsliptol','antoffsetX','antoffsetY','antoffsetZ']
        
        # Now we parse through the config.txt file.
        for line in inputfile:
            
            # Check for input entry with an 'I', then split and format.
            if line[0] == 'I':
                line_inp = line[3:].split()
                
                # Now, let's try to parse parameters meant to be integers.
                if line_inp[0] in integers:
                    
                    try:
                        inps[ line_inp[0] ] = int(line_inp[1])
                    except ValueError:
                        print('Warning, error when reading '+ line_inp[0] +'!')
                        inps[ line_inp[0] ] = line_inp[1]
                
                # then we parse parameters meant to be floats.
                elif line_inp[0] in floats:
                    
                    try:
                        inps[ line_inp[0] ] = float(line_inp[1])
                    except ValueError:
                        print('Warning, error when reading '+ line_inp[0] +'!')
                        inps[ line_inp[0] ] = line_inp[1]
                        
                # For all other parameters, just log them down as they are.
                else:
                    inps[ line_inp[0] ] = line_inp[1]
        
        # Close the file when done
        inputfile.close()
        
        # Check for 4-letter ID of LEO-A
        if len(inps['name1']) != 4:
            self.errtx01.configure(text='Invalid! Check for 4-letter ID.')
        else:
            self.var01.set(inps['name1'])
            self.errtx01.configure(text='')
        
        # Check for 4-letter ID of LEO-B
        if len(inps['name2']) != 4:
            self.errtx02.configure(text='Invalid! Check for 4-letter ID.')
        else:
            self.var02.set(inps['name2'])
            self.errtx02.configure(text='')
        
        # Check for the starting epoch format
        if inps['dtstart'].count('-') == 5:
            self.var03.set(inps['dtstart'])
            self.errtx03.configure(text='')
        else:
            self.errtx03.configure(text='Invalid format, check config.txt!')
            
        # Check for the ending epoch format
        if inps['dtstop'].count('-') == 5:
            self.var04.set(inps['dtstop'])
            self.errtx04.configure(text='')
        else:
            self.errtx04.configure(text='Invalid format, check config.txt!')
        
        # Check for the time step format
        if type(inps['timestep']) == int:
            self.var05.set(inps['timestep'])
            self.errtx05.configure(text='')
        else:
            self.errtx05.configure(text='Error! Timestep is a non-integer!')
        
        # Check for the frequency number.
        if inps['freq'] == 1:
            self.entry06a.select()
            self.entry06b.deselect()
            self.errtx06.configure(text='')
        elif inps['freq'] == 2:
            self.entry06a.deselect()
            self.entry06b.select()
            self.errtx06.configure(text='')
        else:
            self.entry06a.deselect()
            self.entry06b.deselect()
            self.errtx06.configure(text='Invalid frequency! Check config.txt!')
        
        # Check box for Earth's rotation.
        if inps['earthrotation'] == 'True':
            self.entry07.select()
            self.errtx07.configure(text='')
        elif inps['earthrotation'] == 'False':
            self.entry07.deselect()
            self.errtx07.configure(text='')
        else:
            self.entry07.deselect()
            self.errtx07.configure(text='Error! Input must be True or False!')
        
        # Check box for accounting for relativistic effects.
        if inps['relativity'] == 'True':
            self.entry08.select()
            self.errtx08.configure(text='')
        elif inps['relativity'] == 'False':
            self.entry08.deselect()
            self.errtx08.configure(text='')
        else:
            self.entry08.deselect()
            self.errtx08.configure(text='Error! Input must be True or False!')
        
        # Check box for user to save output graph plots.
        if inps['savefigs'] == 'True':
            self.entry09.select()
            self.errtx09.configure(text='')
        elif inps['savefigs'] == 'False':
            self.entry09.deselect()
            self.errtx09.configure(text='')
        else:
            self.entry09.deselect()
            self.errtx09.configure(text='Error! Input must be True or False!')
        
        # Check box for user to save text report for GPS ephemeris.
        if inps['savereport'] == 'True':
            self.entry10.select()
            self.errtx10.configure(text='')
        elif inps['savereport'] == 'False':
            self.entry10.deselect()
            self.errtx10.configure(text='')
        else:
            self.entry10.deselect()
            self.errtx10.configure(text='Error! Input must be True or False!')
        
        # Check box for enabling hatch filtering.
        if inps['hatchfilter'] == 'True':
            self.entry11.select()
            self.errtx11.configure(text='')
        elif inps['hatchfilter'] == 'False':
            self.entry11.deselect()
            self.errtx11.configure(text='')
        else:
            self.entry11.deselect()
            self.errtx11.configure(text='Error! Input must be True or False!')
        
        # Entry for the user to define the length of the hatch filter.
        if type(inps['hatchlength']) == int:
            self.var12.set(inps['hatchlength'])
            if inps['hatchlength'] > 900:
                self.errtx12.configure(text='Warning! Hatch length is long!')
            elif inps['hatchlength'] <= 2:
                self.errtx12.configure(text='Warning! Hatch length is short!')
            else:
                self.errtx12.configure(text='')
        else:
            self.errtx12.configure(text='Error! Hatch length is non-integer!')

        # Entry for the number of standard deviations for cycle slip tolerance.
        if type(inps['cycsliptol']) == float:
            self.var13.set(inps['cycsliptol'])
            self.errtx13.configure(text='')
        else:
            self.errtx13.configure(text='Error! Tolerance must be float!')
            
        # Entry for the length of cycle slip detection filter window.
        if type(inps['cycsliplen']) == int:
            self.var14.set(inps['cycsliplen'])
            if inps['cycsliplen'] > 90:
                self.errtx14.configure(text='Warning! Filter window too long!')
            elif inps['cycsliplen'] <= 2:
                self.errtx14.configure(text='Warning! Filter window is short!')
            else:
                self.errtx14.configure(text='')
        else:
            self.errtx14.configure(text='Error! Filter length invalid!')

        # Antenna offset in X direction.
        if type(inps['antoffsetX']) == float:
            self.var15.set(inps['antoffsetX'])
            self.errtx15.configure(text='')
        else:
            self.errtx15.configure(text='Error! Antenna offset must be float!')
        
        # Antenna offset in Y direction.
        if type(inps['antoffsetY']) == float:
            self.var16.set(inps['antoffsetY'])
            self.errtx16.configure(text='')
        else:
            self.errtx16.configure(text='Error! Antenna offset must be float!')
        
        # Antenna offset in Z direction.
        if type(inps['antoffsetZ']) == float:
            self.var17.set(inps['antoffsetZ'])
            self.errtx17.configure(text='')
        else:
            self.errtx17.configure(text='Error! Antenna offset must be float!')
        
        return None
    
    # Method for writing the entries into the config.txt file.
    def cfg_W(self):
        
        cwd = dirname(dirname(abspath(__file__))) # Current working directory
        iwd = cwd + '\config\config.txt' # Inputs files
        input_r = open(iwd,'r') # Open the config.txt file
        record = [] # Array to record the strings
        error_flag = True # Check if user-defined variables is correct.
        
        # Variables to be recorded based on tkinter entries.
        var_arr = [self.var01, self.var02, self.var03, self.var04, self.var05,
                   self.var06, self.var07, self.var08, self.var09, self.var10,
                   self.var11, self.var12, self.var13, self.var14, self.var15,
                   self.var16, self.var17]
        
        # Key values (to be referred).
        key_arr = ['name1', 'name2', 'dtstart', 'dtstop', 'timestep', 'freq', 
                   'earthrotation', 'relativity', 'savefigs', 'savereport',
                   'hatchfilter', 'hatchlength', 'cycsliptol', 'cycsliplen',
                   'antoffsetX', 'antoffsetY', 'antoffsetZ']
        
        t_f_arr = ['earthrotation', 'relativity', 'savefigs', 'savereport',
                   'hatchfilter']
        
        # Check for 4-letter ID of LEO-A
        if len(self.var01.get()) != 4:
            self.errtx01.configure(text='Invalid! Check for 4-letter ID.')
            error_flag = False
        else:
            self.errtx01.configure(text='')
        
        # Check for 4-letter ID of LEO-B
        if len(self.var02.get()) != 4:
            self.errtx02.configure(text='Invalid! Check for 4-letter ID.')
            error_flag = False
        else:
            self.errtx02.configure(text='')
        
        # Check for the starting epoch format
        if self.var03.get().count('-') == 5:
            self.errtx03.configure(text='')
        else:
            self.errtx03.configure(text='Invalid format, check config.txt!')
            error_flag = False
            
        # Check for the ending epoch format
        if self.var04.get().count('-') == 5:
            self.errtx04.configure(text='')
        else:
            self.errtx04.configure(text='Invalid format, check config.txt!')
            error_flag = False
        
        # Check for the time step format
        if type(self.var05.get()) == int:
            self.errtx05.configure(text='')
        else:
            self.errtx05.configure(text='Error! Timestep is a non-integer!')
            error_flag = False
        
        # Check for the frequency number.
        if self.var06.get() == 1 or self.var06.get() == 2:
            self.errtx06.configure(text='')
        else:
            self.entry06a.deselect()
            self.entry06b.deselect()
            self.errtx06.configure(text='Invalid frequency! Check config.txt!')
            error_flag = False
        
        # Check box for Earth's rotation.
        if self.var07.get() == 1 or self.var07.get() == 0:
            self.errtx07.configure(text='')
        else:
            self.errtx07.configure(text='Error! Input must be True or False!')
            error_flag = False
        
        # Check box for accounting for relativistic effects.
        if self.var08.get() == 1 or self.var08.get() == 0:
            self.errtx08.configure(text='')
        else:
            self.errtx08.configure(text='Error! Input must be True or False!')
            error_flag = False
        
        # Check box for user to save output graph plots.
        if self.var09.get() == 1 or self.var09.get() == 0:
            self.errtx09.configure(text='')
        else:
            self.errtx09.configure(text='Error! Input must be True or False!')
            error_flag = False
        
        # Check box for user to save text report for GPS ephemeris.
        if self.var10.get() == 1 or self.var10.get() == 0:
            self.errtx10.configure(text='')
        else:
            self.errtx10.configure(text='Error! Input must be True or False!')
            error_flag = False
        
        # Check box for enabling hatch filtering.
        if self.var11.get() == 1 or self.var11.get() == 0:
            self.errtx11.configure(text='')
        else:
            self.errtx11.configure(text='Error! Input must be True or False!')
            error_flag = False
        
        # Entry for the user to define the length of the hatch filter.
        if type(self.var12.get()) == int:
            if self.var12.get() > 900:
                self.errtx12.configure(text='Warning! Hatch length is long!')
            elif self.var12.get() <= 2:
                self.errtx12.configure(text='Warning! Hatch length is short!')
            else:
                self.errtx12.configure(text='')
        else:
            self.errtx12.configure(text='Error! Hatch length is non-integer!')
            error_flag = False

        # Entry for the number of standard deviations for cycle slip tolerance.
        if type(self.var13.get()) == float:
            self.errtx13.configure(text='')
        else:
            self.errtx13.configure(text='Error! Tolerance must be float!')
            error_flag = False
            
        # Entry for the length of cycle slip detection filter window.
        if type(self.var14.get()) == int:
            if self.var14.get() > 90:
                self.errtx14.configure(text='Warning! Filter window too long!')
            elif self.var14.get() <= 2:
                self.errtx14.configure(text='Warning! Filter window is short!')
            else:
                self.errtx14.configure(text='')
        else:
            self.errtx14.configure(text='Error! Filter length invalid!')
            error_flag = False

        # Antenna offset in X direction.
        if type(self.var15.get()) == float:
            self.errtx15.configure(text='')
        else:
            self.errtx15.configure(text='Error! Antenna offset must be float!')
        
        # Antenna offset in Y direction.
        if type(self.var16.get()) == float:
            self.errtx16.configure(text='')
        else:
            self.errtx16.configure(text='Error! Antenna offset must be float!')
            error_flag = False
        
        # Antenna offset in Z direction.
        if type(self.var17.get()) == float:
            self.errtx17.configure(text='')
        else:
            self.errtx17.configure(text='Error! Antenna offset must be float!')
            error_flag = False
            
        # If error_flag == 1, then no errors found, proceed with overwriting.
        if error_flag == True:
            
            # Now we parse through the config.txt file.
            for line in input_r:
                
                if line[0] == 'I':
                    words = line.split() # Split string into list of words.
                    key   = words[1] # Get the key from config.txt
                    value = words[2] # Get the value from config.txt
                    
                    # Get the updated value based on the index from key list.
                    value_new = str(var_arr[ key_arr.index(key) ].get())
                    
                    if key in t_f_arr:
                        if value_new == '1':
                            value_new = 'True'
                        if value_new == '0':
                            value_new = 'False'
                    
                    line_new  = line.replace(value, value_new)
                
                else:
                    line_new = line
                
                # Now, record the entries.
                record.append(line_new)
            
            # Close the file when done
            input_r.close()
            
            # Now, we open and overwrite the config.txt file.
            input_w = open(iwd,'w') # Open the config.txt file
            
            for text in record:
                input_w.write(text)
            
            input_w.close()
        
        return None
    
    def run(self):
        
        try:
            self.cfg_W()
            leorun.run()
        except Exception as excpt:
            print('Error in running!')
            print(excpt)
            pass
        
        return None
