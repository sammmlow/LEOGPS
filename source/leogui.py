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
##    This file contains the GUI class which based on Python tkinter.        ##
##    The class will be called in the main LEOGPS python file.               ##
##    (No inputs and outputs, just the GUI class object).                    ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 09-Sep-2021                                              ##
##    Website: https://github.com/sammmlow/LEOGPS                            ##
##    Documentation: https://leogps.readthedocs.io/en/latest/                ##
##                                                                           ##
###############################################################################
###############################################################################

# Import global libraries
import datetime
import tkinter as tk
from PIL import Image, ImageTk
from os.path import dirname, abspath, join
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

# Import local libraries
from source import leorun

class RunGUI:
    
    '''This class represents the entire LEOGPS GUI, as a TKinter object.
    The constructor takes in a single tkinter.Tk() object as the root GUI.
    All buttons in the GUI are linked to the methods described below.
    
    Methods
    -------
    cfg_R( self )
        This method does two things. First, this method checks that all inputs
        in config.txt are correct. Second, it copies the input parameters into
        the GUI (as TKinter variables).
    cfg_W( self )
        This method does two things. First, this method checks that all inputs
        in the GUI are correct. Second, it copies the GUI parameters into the
        config.txt file, overwriting it.
    clr( self )
        Clears all existing relative orbit plots in the LEOGPS GUI.
    run( self )
        Run the LEOGPS program using the leorun.py script and plots the
        relative trajectory.
    '''

    def __init__(self, master):
        
        '''
        Loads the GUI class object, taking a tkinter.Tk() object as input.
        
        Example initialisation:
        >> root = tkinter.Tk()
        >> root_gui = run_gui( root )
        >> root.mainloop()
        '''
        
        # Create the main frame and window.
        master.title('LEOGPS v1.3 - Relative Satellite Navigation in Python')
        screen_height = master.winfo_screenheight()
        screen_width = master.winfo_screenwidth()
        master_geometry  = str(int(screen_width*0.85)) + 'x'
        master_geometry += str(int(screen_height*0.85))
        master.geometry(master_geometry)
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###         Initialisation of text labels and variables           ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        # Initialise the basic text labels (found in the configuration file):
        self.txt01 = 'Spacecraft A 4-letter ID (LEOA)'
        self.txt02 = 'Spacecraft B 4-letter ID (LEOB)'
        self.txt03 = 'Epoch start (YYYY-MM-DD-HH-MN-SS)'
        self.txt04 = 'Epoch final (YYYY-MM-DD-HH-MN-SS)'
        self.txt05 = 'Timestep in seconds (i.e. 30)'
        self.txt06 = 'Enable L1 or L1/L2 processing?'
        self.txt07 = 'Enable code-carrier Hatch filter?'
        self.txt08 = 'Set Hatch filter window length'
        self.txt09 = 'No. of sigmas in cycle slip detection'
        self.txt10 = 'Set cycle slip detection filter length'
        self.txt11 = 'GPS antenna Body-Frame X offset (m)'
        self.txt12 = 'GPS antenna Body-Frame Y offset (m)'
        self.txt13 = 'GPS antenna Body-Frame Z offset (m)'
        self.txt14 = 'Absolute orbit frame (not in plot)' 
        self.txt15 = 'Relative orbit frame (default Hill)'
        
        # Initialise tkinter variables for the entries corresponding to above.
        self.var01 = tk.StringVar() # 4-letter ID of LEO A
        self.var02 = tk.StringVar() # 4-letter ID of LEO B
        self.var03 = tk.StringVar() # Start epoch datetime
        self.var04 = tk.StringVar() # Final epoch datetime
        self.var05 = tk.IntVar()    # Time step per epoch
        self.var06 = tk.IntVar()    # Frequency number
        self.var07 = tk.IntVar()    # Enable flag for hatch filter
        self.var08 = tk.IntVar()    # Hatch filter window length
        self.var09 = tk.DoubleVar() # Cycle slip tolerance STDs
        self.var10 = tk.IntVar()    # Cycle slip detection window
        self.var11 = tk.DoubleVar() # Antenna offset X
        self.var12 = tk.DoubleVar() # Antenna offset Y
        self.var13 = tk.DoubleVar() # Antenna offset Z
        self.var14 = tk.IntVar()    # Frequency number
        self.var15 = tk.IntVar()    # Frequency number
        
        # Initialise an error flag to stop LEOGPS when there are bugs.
        self.error_flag = False
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###             Configure the software logo display               ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        # Define the path to the LEOGPS logo file.
        leogps_logo = dirname(dirname(abspath(__file__)))
        leogps_logo = leogps_logo + '\docs\_static\leogps_logo.png'
        
        # Get the users current screen height.
        screen_height = master.winfo_screenheight()
        
        # Configure the background image and load the logo.
        image = Image.open( leogps_logo )
        photo = ImageTk.PhotoImage(image)
        image_h = photo.height()
        image_w = photo.width()
        logo_scale = image_w / image_h
        logo_height = int(screen_height/6)
        logo_width = int(logo_height*logo_scale)
        image_resize = image.resize(( logo_width, logo_height ))
        image_logo = ImageTk.PhotoImage(image_resize)
        self.logo = tk.Label(image=image_logo)
        self.logo.image = image_logo
        self.logo.grid(row=0, column=0, padx=20, pady=20,
                       sticky='w', columnspan=4)
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###         Add basic buttons for LOAD, SAVE, CLEAR, RUN          ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        # Add a button to read default entries from 'config.txt'.
        self.cfgR = tk.Button(master, text='Load Config', command=self.cfg_R)
        self.cfgR.grid(row=0, column=5, padx=20, pady=2)
        self.cfgR.configure(bg="light blue")
        
        # Add a button to save entries into 'config.txt'.
        self.cfgW = tk.Button(master, text='Save Config', command=self.cfg_W)
        self.cfgW.grid(row=0, column=6, padx=20, pady=2)
        self.cfgW.configure(bg="light blue")
        
        # Add a button to clear the relative orbit plots.
        self.clrBtn = tk.Button(master, text='Clear Plots', command=self.clr)
        self.clrBtn.grid(row=0, column=7, padx=20, pady=5)
        self.clrBtn.configure(bg="light blue")
        
        # Add a button to run LEOGPS.
        self.runBtn = tk.Button(master, text='Run LEOGPS', command=self.run)
        self.runBtn.grid(row=0, column=8, padx=20, pady=2)
        self.runBtn.configure(bg="light blue")
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###               Add labels and main data entries                ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        # Input the 4-letter ID of the first spacecraft (i.e. LEOA).
        self.label01 = tk.Label(master, text=self.txt01 )
        self.label01.grid(row=1, column=0, padx=40, pady=2, sticky='w')
        self.entry01 = tk.Entry(master, width=5, textvariable=self.var01)
        self.entry01.grid(row=1, column=1, padx=5, pady=2, sticky='w')
        self.errtx01 = tk.Label(master, text='', fg='red' )
        self.errtx01.grid(row=1, column=4, padx=5, pady=2, sticky='w')
        
        # Input the 4-letter ID of the second spacecraft (i.e. LEOB).
        self.label02 = tk.Label(master, text=self.txt02 )
        self.label02.grid(row=2, column=0, padx=40, pady=2, sticky='w')
        self.entry02 = tk.Entry(master, width=5, textvariable=self.var02)
        self.entry02.grid(row=2, column=1, padx=5, pady=2, sticky='w')
        self.errtx02 = tk.Label(master, text='', fg='red' )
        self.errtx02.grid(row=2, column=4, padx=5, pady=2, sticky='w')
        
        # Input the starting epoch in YYYY-MM-DD-HH-MN-SS.
        self.label03 = tk.Label(master, text=self.txt03 )
        self.label03.grid(row=3, column=0, padx=40, pady=2, sticky='w')
        self.entry03 = tk.Entry(master, width=20, textvariable=self.var03)
        self.entry03.grid(row=3, column=1, padx=5, pady=2, sticky='w',
                          columnspan=3)
        self.errtx03 = tk.Label(master, text='', fg='red' )
        self.errtx03.grid(row=3, column=4, padx=5, pady=2, sticky='w')
        
        # Input the ending epoch in YYYY-MM-DD-HH-MN-SS.
        self.label04 = tk.Label(master, text=self.txt04 )
        self.label04.grid(row=4, column=0, padx=40, pady=2, sticky='w')
        self.entry04 = tk.Entry(master, width=20, textvariable=self.var04)
        self.entry04.grid(row=4, column=1, padx=5, pady=2, sticky='w',
                          columnspan=3)
        self.errtx04 = tk.Label(master, text='', fg='red' )
        self.errtx04.grid(row=4, column=4, padx=5, pady=2, sticky='w')
        
        # Input the timestep in seconds (i.e. 30).
        self.label05 = tk.Label(master, text=self.txt05 )
        self.label05.grid(row=5, column=0, padx=40, pady=2, sticky='w')
        self.entry05 = tk.Entry(master, width=10, textvariable=self.var05)
        self.entry05.grid(row=5, column=1, padx=5, pady=2, sticky='w',
                          columnspan=2)
        self.errtx05 = tk.Label(master, text='', fg='red' )
        self.errtx05.grid(row=5, column=4, padx=5, pady=2, sticky='w')
        
        # Input single or dual frequency processing (i.e. 1 or 2).
        self.label06 = tk.Label(master, text=self.txt06 )
        self.label06.grid(row=6, column=0, padx=40, pady=2, sticky='w')
        self.entry06a = tk.Radiobutton(master, text='L1',
                                       variable=self.var06, value=1)
        self.entry06a.grid(row=6, column=1, padx=5, pady=2, sticky='w')
        self.entry06b = tk.Radiobutton(master, text='L1+2',
                                       variable=self.var06, value=2)
        self.entry06b.grid(row=6, column=2, padx=5, pady=2, sticky='w')
        self.errtx06 = tk.Label(master, text='', fg='red' )
        self.errtx06.grid(row=6, column=4, padx=5, pady=2, sticky='w')
        
        # Code-carrier smoothing via Hatch filtering? (True/False)
        self.label07 = tk.Label(master, text=self.txt07 )
        self.label07.grid(row=7, column=0, padx=40, pady=2, sticky='w')
        self.entry07 = tk.Checkbutton(master, text='T/F',
                                      variable=self.var07,
                                      command=self._callback_hatch)
        self.entry07.grid(row=7, column=1, padx=5, pady=2, sticky='w')
        self.errtx07 = tk.Label(master, text='', fg='red' )
        self.errtx07.grid(row=7, column=4, padx=5, pady=2, sticky='w')
        
        # Window length (number of observations) of the hatch filter?
        self.label08 = tk.Label(master, text=self.txt08 )
        self.label08.grid(row=8, column=0, padx=40, pady=2, sticky='w')
        self.entry08 = tk.Entry(master, width=10, textvariable=self.var08)
        self.entry08.grid(row=8, column=1, padx=5, pady=2, sticky='w',
                          columnspan=2)
        self.errtx08 = tk.Label(master, text='', fg='red' )
        self.errtx08.grid(row=8, column=4, padx=5, pady=2, sticky='w')
        
        # Standard deviation tolerance for cycle slip detection?
        self.label09 = tk.Label(master, text=self.txt09 )
        self.label09.grid(row=9, column=0, padx=40, pady=2, sticky='w')
        self.entry09 = tk.Entry(master, width=10, textvariable=self.var09)
        self.entry09.grid(row=9, column=1, padx=5, pady=2, sticky='w',
                          columnspan=2)
        self.errtx09 = tk.Label(master, text='', fg='red' )
        self.errtx09.grid(row=9, column=4, padx=5, pady=2, sticky='w')
        
        # Length of the sliding window cycle slip detection filter?
        self.label10 = tk.Label(master, text=self.txt10 )
        self.label10.grid(row=10, column=0, padx=40, pady=2, sticky='w')
        self.entry10 = tk.Entry(master, width=10, textvariable=self.var10)
        self.entry10.grid(row=10, column=1, padx=5, pady=2, sticky='w',
                          columnspan=2)
        self.errtx10 = tk.Label(master, text='', fg='red' )
        self.errtx10.grid(row=10, column=4, padx=5, pady=2, sticky='w')
        
        # Set the GPS antenna offset in X-direction (m) for vehicle.
        self.label11 = tk.Label(master, text=self.txt11 )
        self.label11.grid(row=11, column=0, padx=40, pady=2, sticky='w')
        self.entry11 = tk.Entry(master, width=10, textvariable=self.var11,
                                state=tk.DISABLED)
        self.entry11.grid(row=11, column=1, padx=5, pady=2, sticky='w',
                          columnspan=2)
        self.errtx11 = tk.Label(master, text='', fg='red' )
        self.errtx11.grid(row=11, column=4, padx=5, pady=2, sticky='w')
        
        # Set the GPS antenna offset in Y-direction (m) for vehicle.
        self.label12 = tk.Label(master, text=self.txt12 )
        self.label12.grid(row=12, column=0, padx=40, pady=2, sticky='w')
        self.entry12 = tk.Entry(master, width=10, textvariable=self.var12,
                                state=tk.DISABLED)
        self.entry12.grid(row=12, column=1, padx=5, pady=2, sticky='w',
                          columnspan=2)
        self.errtx12 = tk.Label(master, text='', fg='red' )
        self.errtx12.grid(row=12, column=4, padx=5, pady=2, sticky='w')
        
        # Set the GPS antenna offset in Z-direction (m) for vehicle.
        self.label13 = tk.Label(master, text=self.txt13 )
        self.label13.grid(row=13, column=0, padx=40, pady=2, sticky='w')
        self.entry13 = tk.Entry(master, width=10, textvariable=self.var13,
                                state=tk.DISABLED)
        self.entry13.grid(row=13, column=1, padx=5, pady=2, sticky='w',
                          columnspan=2)
        self.errtx13 = tk.Label(master, text='', fg='red' )
        self.errtx13.grid(row=13, column=4, padx=5, pady=2, sticky='w')
        
        # Select the orbit ephemeris coordinate frame.
        self.label14 = tk.Label(master, text=self.txt14 )
        self.label14.grid(row=14, column=0, padx=40, pady=2, sticky='w')
        self.entry14a = tk.Radiobutton(master, text='ITRF',
                                       variable=self.var14, value=1)
        self.entry14a.grid(row=14, column=1, padx=5, pady=2, sticky='w')
        self.entry14b = tk.Radiobutton(master, text='ICRF',
                                       variable=self.var14, value=2)
        self.entry14b.grid(row=14, column=2, padx=5, pady=2, sticky='w')
        self.errtx14 = tk.Label(master, text='', fg='red' )
        self.errtx14.grid(row=14, column=4, padx=5, pady=2, sticky='w')
        
        # Select the relative baseline coordinate frame.
        self.label15 = tk.Label(master, text=self.txt15 )
        self.label15.grid(row=15, column=0, padx=40, pady=2, sticky='w')
        self.entry15a = tk.Radiobutton(master, text='ITRF',
                                       variable=self.var15, value=1)
        self.entry15a.grid(row=15, column=1, padx=5, pady=2, sticky='w')
        self.entry15b = tk.Radiobutton(master, text='ICRF',
                                       variable=self.var15, value=2)
        self.entry15b.grid(row=15, column=2, padx=5, pady=2, sticky='w')
        self.entry15c = tk.Radiobutton(master, text='Hill',
                                       variable=self.var15, value=3)
        self.entry15c.grid(row=15, column=3, padx=5, pady=2, sticky='w')
        self.errtx15 = tk.Label(master, text='', fg='red' )
        self.errtx15.grid(row=15, column=4, padx=5, pady=2, sticky='w')
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###              Configure the plot area in the GUI               ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        # Now, we add a sub-frame in the tkinter GUI so that we can embed the
        # the interactive matplotlib 3D plot of the relative orbit.
        self.toolbarFrame = tk.Frame(master)
        self.toolbarFrame.grid(row=1, column=5, padx=20, pady=10,
                               columnspan=4, rowspan=16)
        
        # Create the 3D axes matplotlib figure object, using the pack() method
        # of tkinter within the toolbarFrame object.
        self.orbFig = Figure(figsize=(5,4), dpi = master.winfo_fpixels('2.0c'),
                             linewidth=8, edgecolor="#DDDDDD")
        self.orbFig.set_tight_layout(True)
        self.orbPlot = FigureCanvasTkAgg(self.orbFig, self.toolbarFrame)
        self.orbPlot.get_tk_widget().pack(expand=True)
        
        # Note, the plotting should happen after the figure object is called.
        self.orbAxis = self.orbFig.add_subplot(projection='3d')
        self.orbAxis.set_xlabel('Hill Frame Cross-Track Axis (m)')
        self.orbAxis.set_ylabel('Hill Frame In-Track Axis (m)')
        self.orbAxis.set_zlabel('Hill Frame Radial Axis (m)')
        
        # At this point, you can insert plots if you want. For example,
        # self.orbAxis.scatter([1,2,3],[1,2,3],[1,2,3])
        
        self.orbPlot.draw()
        
        # Add the matplotlib navigation toolbar.
        self.toolbar = NavigationToolbar2Tk(self.orbPlot, self.toolbarFrame)
        self.toolbar.update()
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###    Finally, define string containers for error and warning    ###
        ###    messages to inform the user if input conditions violate    ###
        ###    formatting or physical principles. If the length of this   ###
        ###    string variable > 0, then it triggers an error message.    ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        self.error_msgprint = '' # Error message to print.
    
    #####################################################################
    #####################################################################
    ###                                                               ###
    ###    Define a call-back function that disables hatch length     ###
    ###         if the hatch filter option is not checked.            ###
    ###                                                               ###
    #####################################################################
    #####################################################################
    
    def _callback_hatch(self):
        
        try:
            _hatch = self.var07.get()
            if _hatch == 0:
                self.var08.set(10)
                self.entry08.config(state = 'disabled')
            if _hatch == 1:
                self.entry08.config(state = 'normal')
        except:
            pass
        finally:
            return None
    
    #########################################################################
    #########################################################################
    ###                                                                   ###
    ###          Method to load default values from config.txt            ###
    ###                                                                   ###
    #########################################################################
    #########################################################################
    
    def cfg_R(self):
        
        # First, ask the user if he/she wishes to proceed.
        cfg_R_msg = 'Load parameters from the "config.txt" file? \n'
        cfg_R_msg += '(This will overwrite existing inputs in the GUI!)'
        cfg_R_ask = tk.messagebox.askyesno('Load Config', cfg_R_msg)
        if cfg_R_ask == False:
            return None
        
        # Else, continue with loading the configuration file.
        cwd = dirname(dirname(abspath(__file__))) # Current working directory
        iwd = join(cwd, 'config', 'config.txt') # Inputs files
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
                        errmsg = 'Error, expected an integer when reading '
                        errmsg = errmsg + line_inp[0] + ' in config.txt! \n'
                        print(errmsg)
                        self.error_msgprint += errmsg
                        inps[ line_inp[0] ] = 0
                
                # then we parse parameters meant to be floats.
                elif line_inp[0] in floats:
                    
                    try:
                        inps[ line_inp[0] ] = float(line_inp[1])
                    except ValueError:
                        errmsg = 'Error, expected a float when reading '
                        errmsg = errmsg + line_inp[0] + ' in config.txt! \n'
                        print(errmsg)
                        self.error_msgprint += errmsg
                        inps[ line_inp[0] ] = 0.0
                        
                # For all other parameters, just log them down as they are.
                else:
                    inps[ line_inp[0] ] = line_inp[1]
        
        # Close the file when done
        inputfile.close()
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###  Parsing through the inputs dictionary to verify parameters   ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        # 1. Check for 4-letter ID of LEO-A
        
        self.var01.set(inps['name1'])
        if len(inps['name1']) != 4:
            errmsg = 'Invalid 4-letter ID for Spacecraft A! \n'
            self.errtx01.configure(text='!')
        else:
            errmsg = ''
            self.errtx01.configure(text='')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 2. Check for 4-letter ID of LEO-B
        
        self.var02.set(inps['name2'])
        if len(inps['name2']) != 4:
            errmsg = 'Invalid 4-letter ID for Spacecraft B! \n'
            self.errtx02.configure(text='!')
        else:
            if inps['name2'] == inps['name1']:
                errmsg = 'Error! Spacecraft A and B cannot have same ID! \n'
                self.errtx02.configure(text='!')
            else:
                errmsg = ''
                self.errtx02.configure(text='')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 3. Check for the starting epoch format
        
        self.var03.set(inps['dtstart'])
        if inps['dtstart'].count('-') == 5:
            
            inp_v03s = inps['dtstart'].split('-')
                
            # Check if it can be converted into a datetime object.
            try:
                inp_v03d = datetime.datetime(int(inp_v03s[0]),
                                             int(inp_v03s[1]),
                                             int(inp_v03s[2]),
                                             int(inp_v03s[3]),
                                             int(inp_v03s[4]),
                                             int(inp_v03s[5]))
                errmsg = ''
                self.errtx03.configure(text='')
            
            # If not, throw an exception and add it to the error log.
            except:
                errmsg = 'Error! Invalid date and time parameters! \n'
                self.errtx03.configure(text='!')
        
        # Else, throw a formatting error.
        else:
            errmsg = 'Error! Invalid date time format in config.txt! \n'
            self.errtx03.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 4. Check for the ending epoch format
        
        self.var04.set(inps['dtstop'])
        if inps['dtstop'].count('-') == 5:
            
            inp_v04s = inps['dtstop'].split('-')
            
            # Check if it can be converted into a datetime object.
            try:
                inp_v04d = datetime.datetime(int(inp_v04s[0]),
                                             int(inp_v04s[1]),
                                             int(inp_v04s[2]),
                                             int(inp_v04s[3]),
                                             int(inp_v04s[4]),
                                             int(inp_v04s[5]))
                
                # Check if the final epoch is after the initial epoch.
                if inp_v04d <= inp_v03d:
                    errmsg = 'Error! The epoch final is before start! \n'
                    self.errtx04.configure(text='!')
                else:
                    errmsg = ''
                    self.errtx04.configure(text='')
            
            # If not, throw an exception and add it to the error log.
            except:
                errmsg = 'Error! Invalid date and time parameters! \n'
                self.errtx04.configure(text='!')
        
        # Else, throw a formatting error.
        else:
            errmsg = 'Error! Invalid date time format! \n'
            self.errtx04.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 5. Check for the time step format
        
        self.var05.set(inps['timestep'])
        if type(inps['timestep']) == int:
            if inps['timestep'] >= 1:
                errmsg = ''
                self.errtx05.configure(text='')
            else:
                errmsg = 'Error! Timestep is zero or negative! \n'
                self.errtx05.configure(text='!')
        else:
            errmsg = 'Error! Timestep is a non-integer! \n'
            self.errtx05.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 6. Check for the frequency number.
        
        if inps['freq'] == 1:
            errmsg = ''
            self.entry06a.select()
            self.entry06b.deselect()
            self.errtx06.configure(text='')
        elif inps['freq'] == 2:
            errmsg = ''
            self.entry06a.deselect()
            self.entry06b.select()
            self.errtx06.configure(text='')
        else:
            errmsg = 'Invalid frequency! Check config.txt! \n'
            self.entry06a.deselect()
            self.entry06b.deselect()
            self.errtx06.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 7. Check box for enabling hatch filtering.
        
        if inps['hatchfilter'] == 'True':
            errmsg = ''
            self.entry07.select()
            self.errtx07.configure(text='')
        elif inps['hatchfilter'] == 'False':
            errmsg = ''
            self.entry07.deselect()
            self.errtx07.configure(text='')
        else:
            errmsg = 'Error! Input in config.txt must be True or False! \n'
            self.entry07.deselect()
            self.errtx07.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 8. Length of the hatch filter.
        
        self.var08.set(inps['hatchlength'])
        if type(inps['hatchlength']) == int:
            if inps['hatchlength'] > 1000:
                errmsg = 'Error! Hatch length exceeds 1000, too long! \n'
                self.errtx08.configure(text='!')
            elif inps['hatchlength'] <= 2:
                errmsg = 'Error! Hatch length too short! \n'
                self.errtx08.configure(text='!')
            else:
                errmsg = ''
                self.errtx08.configure(text='')
        else:
            errmsg = 'Error! Hatch length is non-integer! \n'
            self.errtx08.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 9. Number of standard deviations for cycle slip tolerance.
        
        self.var09.set(inps['cycsliptol'])
        if type(inps['cycsliptol']) == float and inps['cycsliptol'] > 0.0:
            errmsg = ''
            self.errtx09.configure(text='')
        else:
            errmsg = 'Error! Tolerance must be positive float! \n'
            self.errtx09.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 10. Length of cycle slip detection filter window.
        
        self.var10.set(inps['cycsliplen'])
        if type(inps['cycsliplen']) == int:
            if inps['cycsliplen'] > 100:
                errmsg = 'Error! Filter window too long! \n'
                self.errtx10.configure(text='!')
            elif inps['cycsliplen'] <= 2:
                errmsg = 'Error! Filter window is short! \n'
                self.errtx10.configure(text='!')
            else:
                errmsg = ''
                self.errtx10.configure(text='')
        else:
            errmsg = 'Error! Filter length invalid! \n'
            self.errtx10.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 11. Antenna offset in X direction.
        
        self.var11.set(inps['antoffsetX'])
        if type(inps['antoffsetX']) == float:
            errmsg = ''
            self.errtx11.configure(text='')
        else:
            errmsg = 'Error! Antenna offset must be float! \n'
            self.errtx11.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 12. Antenna offset in Y direction.
        
        self.var12.set(inps['antoffsetY'])
        if type(inps['antoffsetY']) == float:
            errmsg = ''
            self.errtx12.configure(text='')
        else:
            errmsg = 'Error! Antenna offset must be float! \n'
            self.errtx12.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 13. Antenna offset in Z direction.
        
        self.var13.set(inps['antoffsetZ'])
        if type(inps['antoffsetZ']) == float:
            errmsg = ''
            self.errtx13.configure(text='')
        else:
            errmsg = 'Error! Antenna offset must be float! \n'
            self.errtx13.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 14. Select the orbit ephemeris coordinate frame.
        
        if inps['frameOrb'] == 'ITRF':
            errmsg = ''
            self.entry14a.select()
            self.entry14b.deselect()
            self.errtx14.configure(text='')
        elif inps['frameOrb'] == 'ICRF':
            errmsg = ''
            self.entry14a.deselect()
            self.entry14b.select()
            self.errtx14.configure(text='')
        else:
            errmsg = 'Invalid frame! Check config.txt! \n'
            self.entry14a.deselect()
            self.entry14b.deselect()
            self.errtx14.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 15. Select the relative baseline coordinate frame.
        
        if inps['frameForm'] == 'ITRF':
            errmsg = ''
            self.entry15a.select()
            self.entry15b.deselect()
            self.entry15c.deselect()
            self.errtx15.configure(text='')
        elif inps['frameForm'] == 'ICRF':
            errmsg = ''
            self.entry15a.deselect()
            self.entry15b.select()
            self.entry15c.deselect()
            self.errtx15.configure(text='')
        elif inps['frameForm'] == 'Hill':
            errmsg = ''
            self.entry15a.deselect()
            self.entry15b.deselect()
            self.entry15c.select()
            self.errtx15.configure(text='')
        else:
            errmsg = 'Invalid frame! Check config.txt! \n'
            self.entry15a.deselect()
            self.entry15b.deselect()
            self.entry15c.deselect()
            self.errtx15.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # Finally, display an error textbox if there are any error messages.
        
        if len(self.error_msgprint) > 0:
            tk.messagebox.showerror("Error with Configuration File!",
                                    self.error_msgprint)
            self.error_msgprint = '' # Reset error message
            self.error_flag = True
        else:
            self.error_flag = False
        
        return None
    
    #########################################################################
    #########################################################################
    ###                                                                   ###
    ###     Method for writing the entries into the config.txt file.      ###
    ###                                                                   ###
    #########################################################################
    #########################################################################
    
    def cfg_W(self):
        
        # First, ask the user if he/she wishes to proceed.
        cfg_W_msg = 'Save GUI parameters into the "config.txt" file? \n'
        cfg_W_msg += '(This will overwrite parameters in "config.txt"!)'
        cfg_W_ask = tk.messagebox.askyesno('Save Config', cfg_W_msg)
        if cfg_W_ask == False:
            return None
        
        # Else, continue with saving the configurations.
        self.error_msgprint = ''
        
        # Get the directory paths.
        cwd = dirname(dirname(abspath(__file__))) # Current working directory
        iwd = join(cwd, 'config', 'config.txt') # Inputs files
        input_r = open(iwd,'r') # Open the config.txt file
        record = [] # Array to record the strings
        
        # Variables to be recorded based on tkinter entries.
        var_arr = [self.var01, self.var02, self.var03, self.var04, self.var05,
                   self.var06, self.var07, self.var08, self.var09, self.var10, 
                   self.var11, self.var12, self.var13, self.var14, self.var15]
        
        # Key values (to be referred).
        
        key_arr = ['name1', 'name2', 'dtstart', 'dtstop', 'timestep', 'freq',
                   'hatchfilter', 'hatchlength', 'cycsliptol', 'cycsliplen',
                   'antoffsetX', 'antoffsetY', 'antoffsetZ',
                   'frameOrb', 'frameForm']
        
        frame_arr = ['frameOrb', 'frameForm']                   
        
        t_f_arr = ['hatchfilter']
        
        #####################################################################
        #####################################################################
        
        # 1. Check for 4-letter ID of LEO-A
        
        try:
            _v01 = self.var01.get() # Exception raised if entry is erroneous
            if len(_v01) != 4:
                errmsg = 'Invalid 4-letter ID for Spacecraft A! \n'
                self.errtx01.configure(text='!')
            else:
                errmsg = ''
                self.errtx01.configure(text='')
        except:
            _v01 = ''
            errmsg = 'Error! Expected a string for the spacecraft A ID! \n'
            self.errtx01.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 2. Check for 4-letter ID of LEO-B
        
        try:
            _v02 = self.var02.get() # Exception raised if entry is erroneous
            if len(_v02) != 4:
                errmsg = 'Invalid 4-letter ID for Spacecraft B! \n'
                self.errtx02.configure(text='!')
            elif _v02 == _v01:
                errmsg = 'Error! Spacecraft A and B cannot have same ID! \n'
                self.errtx02.configure(text='!')
            else:
                errmsg = ''
                self.errtx02.configure(text='')
        except:
            _v02 = ''
            errmsg = 'Error! Expected a string for the spacecraft B ID! \n'
            self.errtx02.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 3. Check for the starting epoch format
        
        try:
            _v03  = self.var03.get() # Exception raised if entry is erroneous
            
            # Check if the date and time are separated by hyphens
            if _v03.count('-') == 5:
                _v03s = _v03.split('-')
                
                # Check if it can be converted into a datetime object.
                try:
                    _v03d = datetime.datetime(int(_v03s[0]),
                                              int(_v03s[1]),
                                              int(_v03s[2]),
                                              int(_v03s[3]),
                                              int(_v03s[4]),
                                              int(_v03s[5]))
                    errmsg = ''
                    self.errtx03.configure(text='')
                
                # If not, throw an exception and add it to the error log.
                except:
                    errmsg = 'Error! Invalid date and time parameters! \n'
                    self.errtx03.configure(text='!')
            
            # Else, throw a formatting error.
            else:
                errmsg = 'Error! Invalid date time format! \n'
                self.errtx03.configure(text='!')
        
        # If it can't be read as a string for some reason.
        except:
            errmsg = 'Error! Unable to read date and time as a string! \n'
            self.errtx03.configure(text='!')
        
        # Add the error log to the error message print.
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 4. Check for the ending epoch format
        
        try:
            _v04  = self.var04.get() # Exception raised if entry is erroneous
            
            # Check if the date and time are separated by hyphens
            if _v04.count('-') == 5:
                _v04s = _v04.split('-')
                
                # Check if it can be converted into a datetime object.
                try:
                    _v04d = datetime.datetime(int(_v04s[0]),
                                              int(_v04s[1]),
                                              int(_v04s[2]),
                                              int(_v04s[3]),
                                              int(_v04s[4]),
                                              int(_v04s[5]))
                    
                    # Check if the final epoch is after the initial epoch.
                    if _v04d <= _v03d:
                        errmsg = 'Error! The epoch final is before start! \n'
                        self.errtx04.configure(text='!')
                    else:
                        _tdelta = _v04d - _v03d
                        _duration  = (_tdelta.days*86400) + (_tdelta.seconds)
                        errmsg = ''
                        self.errtx04.configure(text='')
                
                # If not, throw an exception and add it to the error log.
                except:
                    errmsg = 'Error! Invalid date and time parameters! \n'
                    self.errtx04.configure(text='!')
            
            # Else, throw a formatting error.
            else:
                errmsg = 'Error! Invalid date time format! \n'
                self.errtx04.configure(text='!')
        
        # If it can't be read as a string for some reason.
        except:
            errmsg = 'Error! Unable to read date and time as a string! \n'
            self.errtx04.configure(text='!')
        
        # Add the error log to the error message print.
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 5. Check for the time step format
        
        try:
            _v05 = self.var05.get() # Exception raised if entry is erroneous
            self.var05.set(_v05) # Needed to force as integer if float.
            if _v05 <= 0:
                errmsg = 'Error! Time step cannot be zero or negative! \n'
                self.errtx05.configure(text='!')
            elif _v05 > int(_duration):
                errmsg = 'Error! Step size is greater than whole duration! \n'
                self.errtx05.configure(text='!')
            elif _v05 > 31536000:
                errmsg = 'Error! Step size is longer than a year! \n'
                self.errtx05.configure(text='!')
            else:
                errmsg = ''
                self.errtx05.configure(text='')
        except:
            _v05 = 1
            self.var05.set(_v05) # Reset to default 86400s (1 day)
            errmsg = 'Error! Expected integer for scenario step size (s)! \n'
            self.errtx05.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 6. Check for the frequency number.
        
        try:
            _v06 = self.var06.get()
            if _v06 == 1 or _v06 == 2:
                errmsg = ''
                self.errtx06.configure(text='')
            else:
                _v06 = 1 # Default
                errmsg = 'Invalid frequency selection! \n'
                self.errtx06.configure(text='!')
                self.var06.set(_v06)
                self.entry06a.select()
                self.entry06b.deselect()
        except:
            _v06 = 1 # Default
            errmsg = 'Invalid frequency selection! \n'
            self.errtx06.configure(text='!')
            self.var06.set(_v06)
            self.entry06a.select()
            self.entry06b.deselect()
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 7. Check box for enabling hatch filtering.
        
        try:
            _v07 = self.var07.get()
            if _v07 == 1 or _v07 == 0:
                errmsg = ''
                self.errtx07.configure(text='')
            else:
                _v07 = 0 # Default
                errmsg = 'Error! Hatch filter option not True or False! \n'
                self.errtx07.configure(text='!')
                self.var07.set(_v07)
                self.entry07.deselect()
        except:
            _v07 = 0 # Default
            errmsg = 'Error! Hatch filter option unreadable! \n'
            self.errtx07.configure(text='!')
            self.var07.set(_v07)
            self.entry07.deselect()
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 8. Entry for the user to define the length of the hatch filter.
        
        try:
            _v08 = self.var08.get()
            self.var08.set(_v08) # Needed to force as integer if float.
            if _v08 > 1000:
                errmsg = 'Error! Hatch length > 1000, too long! \n'
                self.errtx08.configure(text='!')
            elif _v08 <= 2:
                errmsg = 'Error! Hatch length < 3, too short! \n'
                self.errtx08.configure(text='!')
            else:
                errmsg = ''
                self.errtx08.configure(text='')
        except:
            _v08 = 10
            errmsg = 'Error! Hatch filter length not integer! Reset to 10. \n'
            self.errtx08.configure(text='!')
            self.var08.set(_v08)
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 9. Number of standard deviations for cycle slip tolerance.
        
        try:
            _v09 = float(self.var09.get())
            if _v09 > 6.0:
                errmsg = 'Error! Tolerance > 6-Sigma is pointless! \n'
                self.errtx09.configure(text='!')
            elif _v09 <= 0.0:
                errmsg = 'Error! Tolerance cannot be negative or zero! \n'
                self.errtx09.configure(text='!')
            else:
                errmsg = ''
                self.errtx09.configure(text='')
        except:
            _v09 = 3.0
            errmsg = 'Error! Tolerance not float! Reset to 3-Sigma. \n'
            self.errtx09.configure(text='!')
            self.var09.set(_v09)
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 10. Length of cycle slip detection filter window.
        
        try:
            _v10 = self.var10.get()
            self.var10.set(_v10) # Needed to force as integer if float.
            if _v10 > 100:
                errmsg = 'Error! Cycle slip filter window > 100, too long! \n'
                self.errtx10.configure(text='!')
            elif _v10 <= 2:
                errmsg = 'Error! Cycle slip filter window <= 2, too short! \n'
                self.errtx10.configure(text='!')
            else:
                errmsg = ''
                self.errtx10.configure(text='')
        except:
            _v10 = 10
            errmsg = 'Error! Cycle slip filter window length unreadable! \n'
            self.errtx10.configure(text='!')
            self.var10.set(_v10)
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 11. Antenna offset in X direction.
        
        try:
            _v11 = float(self.var11.get())
            errmsg = ''
            self.errtx11.configure(text='')
        except:
            _v11 = 0.0
            errmsg = 'Error! Antenna X offset not float! Reset to 0.0. \n'
            self.errtx11.configure(text='!')
            self.var11.set(_v11)
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 12. Antenna offset in Y direction.
        
        try:
            _v12 = float(self.var12.get())
            errmsg = ''
            self.errtx12.configure(text='')
        except:
            _v12 = 0.0
            errmsg = 'Error! Antenna Y offset not float! Reset to 0.0. \n'
            self.errtx12.configure(text='!')
            self.var12.set(_v12)
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 13. Antenna offset in Z direction.
        
        try:
            _v13 = float(self.var13.get())
            errmsg = ''
            self.errtx13.configure(text='')
        except:
            _v13 = 0.0
            errmsg = 'Error! Antenna Z offset not float! Reset to 0.0. \n'
            self.errtx13.configure(text='!')
            self.var13.set(_v13)
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 14. Select the output orbit ephemeris coordinate frame
        
        try:
            _v14 = self.var14.get()
            if _v14 == 1 or _v14 == 2:
                errmsg = ''
                self.errtx14.configure(text='')
            else:
                _v14 = 1 # Default
                errmsg = 'Invalid coordinate frame! \n'
                self.errtx14.configure(text='!')
                self.var14.set(_v14)
                self.entry14a.select()
                self.entry14b.deselect()
        except:
            _v14 = 1 # Default
            errmsg = 'Invalid coordinate frame! \n'
            self.errtx14.configure(text='!')
            self.var14.set(_v14)
            self.entry14a.select()
            self.entry14b.deselect()
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 15. Select the output relative baseline coordinate frame
        
        try:
            _v15 = self.var15.get()
            if _v15 == 1 or _v15 == 2 or _v15 == 3:
                errmsg = ''
                self.errtx15.configure(text='')
            else:
                _v15 = 1 # Default
                errmsg = 'Invalid coordinate frame! \n'
                self.errtx15.configure(text='!')
                self.var15.set(_v15)
                self.entry15a.select()
                self.entry15b.deselect()
                self.entry15c.deselect()
        except:
            _v15 = 1 # Default
            errmsg = 'Invalid coordinate frame! \n'
            self.errtx15.configure(text='!')
            self.var15.set(_v15)
            self.entry15a.select()
            self.entry15b.deselect()
            self.entry15c.deselect()
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # Finally, display an error textbox if there are any error messages.
        if len(self.error_msgprint) > 0:
            tk.messagebox.showerror("Error with Inputs!",
                                    self.error_msgprint)
            self.error_msgprint = '' # Reset error message
            self.error_flag = True
            return None
        
        else:
            
            self.error_flag = False
            
            # Now we parse through the config.txt file.
            for line in input_r:
                
                if line[0] == 'I':
                    words = line.split() # Split string into list of words.
                    key   = words[1] # Get the key from config.txt
                    value = words[2] # Get the value from config.txt
                    
                    # Get the updated value based on the index from key list.
                    value_new = str(var_arr[ key_arr.index(key) ].get())
                    
                    if key in t_f_arr:
                        if value_new == '1': # This is from Tkinter
                            value_new = 'True'
                        if value_new == '0': # This is from Tkinter
                            value_new = 'False'
                    
                    if key in frame_arr:
                        if value_new == '1': # This is from Tkinter
                            value_new = 'ITRF'
                        if value_new == '2': # This is from Tkinter
                            value_new = 'ICRF'
                        if value_new == '3': # This is from Tkinter
                            value_new = 'Hill'
                    
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
    
    #########################################################################
    #########################################################################
    ###                                                                   ###
    ###    Clears all existing relative orbit plots in the LEOGPS GUI.    ###
    ###                                                                   ###
    #########################################################################
    #########################################################################
    
    def clr(self):
        
        # First, ask the user if he/she wishes to proceed.
        clr_msg = 'Clear the plot in the GUI?'
        clr_ask = tk.messagebox.askyesno('Clear Plots', clr_msg)
        if clr_ask == False:
            return None
        
        # Else, continue with clearing the plots.
        self.orbAxis.clear()
        self.orbPlot.draw()
        
        return None
    
    #########################################################################
    #########################################################################
    ###                                                                   ###
    ###         Run the LEOGPS program using the leorun.py script,        ###
    ###                 and plots the relative trajectory.                ###
    ###                                                                   ###
    #########################################################################
    #########################################################################
    
    def run(self):
        
        # First, ask the user if he/she wishes to proceed.
        run_msg = 'Save Parameters and Run LEOGPS?'
        run_ask = tk.messagebox.askyesno('Run LEOGPS?', run_msg)
        if run_ask == False:
            return None
        
        try:
            
            # Write the GUI config and run.
            self.cfg_W()
            
            # `inps` is a dictionary of inputs created by `inpxtr.inpxtr()`
            if self.error_flag == False:
                results, inps = leorun.run()
            else:
                print('There is an error with the configurations in LEOGPS!')
                return None
            
            # Extract out the results and the inputs.
            # results[i] = [pos1, vel1, dop1, cb1, pos2, vel2, dop2, cb2, base]
            # dimensions :   1x3   1x3   1x3  1x1   1x3   1x3   1x3  1x1   1x3
            
            rpx, rpy, rpz = [], [], []
            for t in results:
                rpx.append( results[t][8][0] ) # Radial Values
                rpy.append( results[t][8][1] ) # Cross-Track Values
                rpz.append( results[t][8][2] ) # In-Track Values
            
            # Perform relative orbit plots in the Euler-Hill Frame.       
            self.orbAxis.scatter( rpy, # Cross-Track Axis
                                  rpz, # In-Track Axis
                                  rpx, # Radial Axis
                                  s=8, # Marker size
                                  alpha = 0.25, # Transparency
                                  label='Relative Orbit in Hill-Frame')
            
            # Update the plots
            self.orbPlot.draw()
            
            # Get the current axes limits on relative orbit plots.
            axOrbR_axes_limits = [abs(self.orbAxis.get_xlim()[0]),
                                  abs(self.orbAxis.get_xlim()[1]),
                                  abs(self.orbAxis.get_ylim()[0]),
                                  abs(self.orbAxis.get_ylim()[1]),
                                  abs(self.orbAxis.get_zlim()[0]),
                                  abs(self.orbAxis.get_zlim()[1])]
            
            # Using axOrbR_axes_limits above, we can find the minimum and
            # maximum axes limits and the span of values.
            axOrbR_axes_max  = max(axOrbR_axes_limits)
            axOrbR_axes_span = axOrbR_axes_max * 2

            # Scale all axes equally for relative orbit plots
            self.orbAxis.set_xlim( -1 * axOrbR_axes_max, axOrbR_axes_max )
            self.orbAxis.set_ylim( -1 * axOrbR_axes_max, axOrbR_axes_max )
            self.orbAxis.set_zlim( -1 * axOrbR_axes_max, axOrbR_axes_max )
            
            # It is important that the XYZ axes in the VVLH (relative 
            # orbit) frame is scaled the same, else it is difficult to 
            # interpret the relative separations on different scales.
            
            # Plot the chief satellite as a tri-axial quiver in VVLH frame.
            axOrbR0 = axOrbR_axes_span * 0.1
            self.orbAxis.quiver( 0,0,0,1,0,0, length = axOrbR0,
                                color = 'r', arrow_length_ratio=0.3 )
            self.orbAxis.quiver( 0,0,0,0,1,0, length = axOrbR0,
                                color = 'r', arrow_length_ratio=0.3 )
            self.orbAxis.quiver( 0,0,0,0,0,1, length = axOrbR0,
                                color = 'r', arrow_length_ratio=0.3 )
            
            # Update the plots
            self.orbPlot.draw()
            
        except Exception as excpt:
            print('Error in running!')
            print(excpt)
            pass
        
        return None
