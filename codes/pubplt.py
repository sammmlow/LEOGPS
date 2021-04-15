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
## FILE DESCRIPTION:                                                         ##
##                                                                           ##
## Publishing and plotting of interpolated GPS ephemeris and clock data.     ##
##                                                                           ##
## REMARKS:                                                                  ##
## This program is run as a sub-routine in gpsxtr.py.                        ##
## Saving of GPS PVT graphs and reports can be disabled in config.txt,       ##
## when *savefigs* and *savereport* is set as True/False                     ##
##                                                                           ##
## AUTHOR MODIFIED: 10-03-2021, by Samuel Y.W. Low                           ##
##                                                                           ##
###############################################################################
###############################################################################
'''

import datetime
import matplotlib.pyplot as plt
from decimal import Decimal

def gps_report(gpsdata, goodsats, inps):
    
    cwd = inps['cwd'] # Get current main working directory
    
    file_path = open(cwd+'\\output\\gps_report\\GPS_Report.txt', 'w')
    
    line = 'G         '
    line += 'Pos_X (km)       '
    line += 'Pos_Y (km)       '
    line += 'Pos_Z (km)     '
    line += 'Vel_X (km/s)     '
    line += 'Vel_Y (km/s)     '
    line += 'Vel_Z (km/s)     '
    line += '    Clk_Bias \n'
    file_path.write(line)
    
    # It's all string formatting from here... nothing scientific.
    for t in gpsdata:
        
        line = '\n'
        line += '*  '
        line += str(t.year) + ' '
        line += str(t.month) + ' '
        line += str(t.day) + ' '
        line += str(t.hour) + ' '
        line += str(t.minute) + ' '
        line += str(t.second) + '\n'
        file_path.write(line)

        for p in goodsats:
            
            if p < 10:
                pstr = '0' + str(p)
            else:
                pstr = str(p)
                
            # Write in position information
            line = 'G' + pstr + ' '
            
            for coord in ['x','y','z']:
                pos = str(gpsdata[t][p]['p' + coord])
                dot = pos.index('.')
                if len(pos[dot:]) > 7:
                    pos = pos[:dot+7]
                while len(pos[dot:]) < 7:
                    pos = pos + '0'
                while len(pos[:dot]) < 9:
                    pos = ' ' + pos
                    dot = pos.index('.')
                pos = pos + ' '
                line += pos
                            
            for coord in ['x','y','z']:
                vel = str(gpsdata[t][p]['v' + coord])
                dot = vel.index('.')
                if len(vel[dot:]) > 7:
                    vel = vel[:dot+7]
                while len(vel[dot:]) < 7:
                    vel = vel + '0'
                while len(vel[:dot]) < 9:
                    vel = ' ' + vel
                    dot = vel.index('.')
                vel = vel + ' '
                line += vel
            
            b = '%.9E' % Decimal(str(gpsdata[t][p]['clkb']))
            dot = b.index('.')
            while len(b[:dot]) < 2:
                b = ' ' + b
                dot = b.index('.')
                
            line += str(b)
            line += ' \n'
            file_path.write(line)
            
    file_path.close()
    return None

def gps_graphs(SV, t_usr_dt, t_usr_ss, gpsdata, inps):
    
    cwd = inps['cwd'] # Get current main working directory
    
    # Turn off interactive plotting
    plt.ioff()
    
    # Initialise the 1x3 subplot for PVT data.
    fig, (ax1, ax2, ax3) = plt.subplots(3,1,figsize=(12,8))
        
    # Get the positions, velocities, and clock biases.
    px = [gpsdata[t][SV]['px'] for t in t_usr_dt]
    py = [gpsdata[t][SV]['py'] for t in t_usr_dt]
    pz = [gpsdata[t][SV]['pz'] for t in t_usr_dt]
    vx = [gpsdata[t][SV]['vx'] for t in t_usr_dt]
    vy = [gpsdata[t][SV]['vy'] for t in t_usr_dt]
    vz = [gpsdata[t][SV]['vz'] for t in t_usr_dt]
    clkb = [gpsdata[t][SV]['clkb'] for t in t_usr_dt]
    
    # Position plots
    ax1.set_title('SV ' + str(SV) + ' Position (km)')
    ax1.plot(t_usr_ss, px, c = 'r', label='X')
    ax1.plot(t_usr_ss, py, c = 'g', label='Y')
    ax1.plot(t_usr_ss, pz, c = 'b', label='Z')
    ax1.legend(loc='lower right')
    
    # Velocity plots
    ax2.set_title('SV ' + str(SV) + ' Velocity (km/s)')
    ax2.plot(t_usr_ss, vx, c = 'r', label='X')
    ax2.plot(t_usr_ss, vy, c = 'g', label='Y')
    ax2.plot(t_usr_ss, vz, c = 'b', label='Z')
    ax2.legend(loc='lower right')
    
    # Clock bias plots
    ax3.set_title('SV ' + str(SV) + ' Clock bias (s)')
    ax3.plot(t_usr_ss, clkb, c = 'k', label='Bias')
    ax3.legend(loc="right")
    
    # Tight-spaced plot
    plt.tight_layout()
    plt.savefig(cwd + '\\output\\gps_plots\\GPS_SV' + str(SV) + '_PVT.png')
    
    # Close this figure
    plt.close(fig)
    
    return None

def leo_results(results, inps):
    
    print('Saving final report on both LEOs and their baselines \n')
    
    cwd = inps['cwd'] # Get current main working directory
    file_path = open(cwd+'\\output\\LEOGPS_Results.txt', 'w')
    
    line  = 'Date      '
    line += '     Time     '
    
    # Headers for LEO 1
    line += inps['name1'] + '_PosX     '
    line += inps['name1'] + '_PosY     '
    line += inps['name1'] + '_PosZ     '
    line += inps['name1'] + '_VelX     '
    line += inps['name1'] + '_VelY     '
    line += inps['name1'] + '_VelZ     '
    line += inps['name1'] + '_GDOP     '
    line += inps['name1'] + '_PDOP     '
    line += inps['name1'] + '_TDOP         '
    line += inps['name1'] + '_ClkB     '
    
    # Headers for LEO 2
    line += inps['name2'] + '_PosX     '
    line += inps['name2'] + '_PosY     '
    line += inps['name2'] + '_PosZ     '
    line += inps['name2'] + '_VelX     '
    line += inps['name2'] + '_VelY     '
    line += inps['name2'] + '_VelZ     '
    line += inps['name2'] + '_GDOP     '
    line += inps['name2'] + '_PDOP     '
    line += inps['name2'] + '_TDOP         '
    line += inps['name2'] + '_ClkB     '
    
    # Headers for baseline information
    line += 'RelativeX     '
    line += 'RelativeY     '
    line += 'RelativeZ     '
    line += '\n'
    file_path.write(line)
    
    # It's all string formatting from here... nothing scientific.
    for t in results:
        
        line  = str(t) # Date-time string (dictionary key)
        
        # Within each vector...
        for vector in results[t]:
            
            # Check if the vector is a 1x3 POS/VEL/DOP
            if len(vector) >= 3:
                for value in vector[:3]:
                    svalue = str(value)
                    dot = svalue.index('.')
                    if len(svalue[dot:]) > 4:
                        svalue = svalue[:dot+4]
                    while len(svalue[dot:]) < 4:
                        svalue = svalue + '0'
                    while len(svalue[:dot]) < 10:
                        svalue = ' ' + svalue
                        dot = svalue.index('.')
                    line += svalue
                    
            # Or the clock bias entry (1x1)
            else:
                for value in vector[:1]:
                    svalue = str(value/299792458.0)
                    dot = svalue.index('.')
                    
                    # Check if clock bias is in standard notation
                    if 'e' in svalue:
                        edot = svalue.index('e')
                        svalue1 = svalue[:dot]
                        svalue2 = svalue[dot:edot]
                        if len(svalue2) > 7:
                            svalue2 = svalue2[:7]
                        svalue3 = svalue[edot:]
                        svalue = svalue1 + svalue2 + svalue3
                    
                    # Else, if it is in decimal...
                    else:
                        if len(svalue[dot:]) > 11:
                            svalue = svalue[:dot+11]
                        while len(svalue[dot:]) < 11:
                            svalue = svalue + '0'
                    
                    while len(svalue[:dot]) < 7:
                        svalue = ' ' + svalue
                        dot = svalue.index('.')
                    line += svalue
                
        line += ' \n'
        file_path.write(line)
    
    file_path.close()
    
    print('Completed processing in LEOGPS! Output file stored:')
    print(cwd+'\\output\\LEOGPS_Results.txt \n')
    
    return None