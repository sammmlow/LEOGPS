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

Installation
============

First, find the LEOGPS GitHub repository in this `GitHub link <https://github.com/sammmlow/LEOGPS>`_, and download it. Alternatively, if you have Git installed, you can open Command Prompt or your Git Bash, enter the directory of your choice, and type:

.. code-block:: bash
   
   git clone https://github.com/sammmlow/LEOGPS.git

Second, you should do a pip install in your terminal of Martin Valgur's Pythonic translation of `Hatanaka (de)compression in Python <https://pypi.org/project/hatanaka/>`_:

.. code-block:: bash

   pip install hatanaka

The Hatanaka library in Python was kindly contributed by Martin Valgur in v1.1, and replaces the older "RNX2CRX" and "GZIP" Windows-only executables, making the (de)compression possible across all operating systems.

That's it! No further setup is needed, unless you need any of the other package dependencies.

.. note:: Package dependencies include: 
   copy, datetime, decimal, hatanaka, math, matplotlib, numpy, os, pathlib, PIL, tkinter, unlzw3, urllib, warnings

Next, we will run and explain the setup behind the default native scenario packaged in LEOGPS - a formation flying scenario at ~200km in-track baseline of the GRACE mission.
