..
   ###########################################################################
   ###########################################################################
   ##                                                                       ##
   ##     _    ___  ___   ___ ___ ___                                       ##
   ##    | |  | __ /   \ / __| _ | __|                                      ##
   ##    | |__| __  ( ) | (_ |  _|__ \                                      ##
   ##    |____|___ \___/ \___|_| \___/                                      ##
   ##                                    v 1.2 (Stable)                     ##
   ##                                                                       ##
   ###########################################################################
   ###########################################################################

.. image:: /_static/leogps_logo.png

|

Processing Flow 
===============

A summary of the processing flow in LEOGPS is given by the flow chart below. This is the main processing tree that runs internally behind the GUI, in the **run()** function of the module **leorun.py**.

.. figure:: /_figures/leogps_flowchart.png
   :align: center

If the user wishes to write their own processing workflow (by overwriting **leorun.py**), it is advisable to refer to the API reference on the next page.