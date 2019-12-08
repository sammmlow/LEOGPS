#!/usr/bin/env python3
'''
###############################################################################
###############################################################################
##                                                                           ##
##     _    ___  ___   ___ ___ ___                                           ##
##    | |  | __ /   \ / __| _ | __|                                          ##
##    | |__| __  ( ) | (_ |  _|__ \                                          ##
##    |____|___ \___/ \___|_| \___/                                          ##
##                                    v 0.1 (Alpha)                          ##
##                                                                           ##
## FILE DESCRIPTION:                                                         ##
##                                                                           ##
## Publishing and plotting of interpolated GPS ephemeris and clock data.     ##
##                                                                           ##
## REMARKS:                                                                  ##
## Run as a sub-routine in gpsxtr.py, and the main leogps.py.                ##
## Saving of GPS PVT graphs and reports can be disabled in config.txt,       ##
## when *savefigs* and *savereport* is set as True/False                     ##
##                                                                           ##
## AUTHOR MODIFIED: 30-11-2019, by Samuel Y.W. Low                           ##
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
    
    line = 'G   '
    line += '  Pos_X (km)  '
    line += '  Pos_Y (km)  '
    line += '  Pos_Z (km)  '
    line += ' Vel_X (km/s) '
    line += ' Vel_Y (km/s) '
    line += ' Vel_Z (km/s) '
    line += '  Clk_Bias \n'
    file_path.write(line)
    
    # It's all string formatting from here... nothing scientific.
    for t in range(0,len(gpsdata['t'])):
        
        line = '\n'
        line += '*  '
        line += str(gpsdata['t'][t].year) + ' '
        line += str(gpsdata['t'][t].month) + ' '
        line += str(gpsdata['t'][t].day) + ' '
        line += str(gpsdata['t'][t].hour) + ' '
        line += str(gpsdata['t'][t].minute) + ' '
        line += str(gpsdata['t'][t].second) + '\n'
        file_path.write(line)

        for p in goodsats:
            
            if p < 10:
                pstr = '0' + str(p)
            else:
                pstr = str(p)
                
            # Write in position information
            line = 'G' + pstr + ' '
            
            for coord in ['x','y','z']:
                pos = str(gpsdata[p]['p' + coord][t])
                dot = pos.index('.')
                if len(pos[dot:]) > 4:
                    pos = pos[:dot+4]
                while len(pos[dot:]) < 4:
                    pos = pos + '0'
                while len(pos[:dot]) < 9:
                    pos = ' ' + pos
                    dot = pos.index('.')
                pos = pos + ' '
                line += pos
                            
            for coord in ['x','y','z']:
                vel = str(gpsdata[p]['v' + coord][t])
                dot = vel.index('.')
                if len(vel[dot:]) > 4:
                    vel = vel[:dot+4]
                while len(vel[dot:]) < 4:
                    vel = vel + '0'
                while len(vel[:dot]) < 9:
                    vel = ' ' + vel
                    dot = vel.index('.')
                vel = vel + ' '
                line += vel
            
            b = '%.9E' % Decimal(str(gpsdata[p]['clkb'][t]))
            dot = b.index('.')
            while len(b[:dot]) < 2:
                b = ' ' + b
                dot = b.index('.')
                
            line += str(b)
            line += ' \n'
            file_path.write(line)
            
    file_path.close()
    return None

def gps_graphs(prn, t_usr_ls, gpsdata, inps):
    
    cwd = inps['cwd'] # Get current main working directory
    
    # Turn off interactive plotting
    plt.ioff()
    
    # Initialise the 1x3 subplot for PVT data.
    fig, (ax1, ax2, ax3) = plt.subplots(3,1,figsize=(12,8))
    
    # Position plots
    ax1.set_title('PRN ' + str(prn) + ' Position (km)')
    ax1.plot(t_usr_ls, gpsdata[prn]['px'], c = 'r', label='X')
    ax1.plot(t_usr_ls, gpsdata[prn]['py'], c = 'g', label='Y')
    ax1.plot(t_usr_ls, gpsdata[prn]['pz'], c = 'b', label='Z')
    ax1.legend(loc='lower right')
    
    # Velocity plots
    ax2.set_title('PRN ' + str(prn) + ' Velocity (km/s)')
    ax2.plot(t_usr_ls, gpsdata[prn]['vx'], c = 'r', label='X')
    ax2.plot(t_usr_ls, gpsdata[prn]['vy'], c = 'g', label='Y')
    ax2.plot(t_usr_ls, gpsdata[prn]['vz'], c = 'b', label='Z')
    ax2.legend(loc='lower right')
    
    # Clock bias plots
    ax3.set_title('PRN ' + str(prn) + ' Clock bias (s)')
    ax3.plot(t_usr_ls, gpsdata[prn]['clkb'], c = 'k', label='Bias')
    ax3.legend(loc="right")
    
    # Tight-spaced plot
    plt.tight_layout()
    plt.savefig(cwd + '\\output\\gps_plots\\GPS_PRN' + str(prn) + '_PVT.png')
    
    # Close this figure
    plt.close(fig)
    
    return None

def leo_results(results, inps):
    
    print('Saving final report on both LEOs and their baselines \n')
    
    cwd = inps['cwd'] # Get current main working directory
    file_path = open(cwd+'\\output\\LEOGPS_Results.txt', 'w')
    
    line  = 'Date      '
    line += '    Time     '
    
    # Headers for LEO 1
    line += inps['name1'] + '_PosX     '
    line += inps['name1'] + '_PosY     '
    line += inps['name1'] + '_PosZ     '
    line += inps['name1'] + '_VelX     '
    line += inps['name1'] + '_VelY     '
    line += inps['name1'] + '_VelZ     '
    line += inps['name1'] + '_GDOP     '
    line += inps['name1'] + '_PDOP     '
    line += inps['name1'] + '_TDOP     '
    
    # Headers for LEO 2
    line += inps['name2'] + '_PosX     '
    line += inps['name2'] + '_PosY     '
    line += inps['name2'] + '_PosZ     '
    line += inps['name2'] + '_VelX     '
    line += inps['name2'] + '_VelY     '
    line += inps['name2'] + '_VelZ     '
    line += inps['name2'] + '_GDOP     '
    line += inps['name2'] + '_PDOP     '
    line += inps['name2'] + '_TDOP     '
    
    # Headers for baseline information
    line += 'RelativeX     '
    line += 'RelativeY     '
    line += 'RelativeZ     '
    line += '\n'
    file_path.write(line)
    
    # It's all string formatting from here... nothing scientific.
    for t in results:
        line  = str(t)
        for vector in results[t]:
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
        line += ' \n'
        file_path.write(line)
    
    file_path.close()
    return None