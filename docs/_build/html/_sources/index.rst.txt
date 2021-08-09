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

.. |docs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat-square
   :target: https://leogps.readthedocs.io/en/latest/

.. |license| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square
   :target: https://github.com/sammmlow/LEOGPS/blob/master/LICENSE
   
.. |orcid| image:: https://img.shields.io/badge/ID-0000--0002--1911--701X-a6ce39.svg
   :target: https://orcid.org/0000-0002-1911-701X/
   
.. |linkedin| image:: https://img.shields.io/badge/LinkedIn-sammmlow-blue.svg
   :target: https://www.linkedin.com/in/sammmlow

.. image:: /_static/leogps_logo.png

|

:Github: https://github.com/sammmlow/LEOGPS
:Documents: https://leogps.readthedocs.io/en/latest/
:Version: 1.2 (Latest)
:Author: Samuel Y. W. Low

|docs| |license| |linkedin| |orcid|

LEOGPS
======

LEOGPS is an open-source Python software which performs relative satellite navigation between **two formation flying satellites**, with the objective of high accuracy relative positioning. Specifically, LEOGPS solves for the double-differenced baseline (using float ambiguity resolution) between satellites flying in formation in Low Earth Orbit (LEO). As such, the relative positioning accuracy diminishes with increasing formation baseline lengths.

.. note:: For formations with baselines kept below 200km, the relative positioning accuracy using the double-differenced kinematic technique can be kept under 1 meter for an automotive or space grade GPS receiver.

LEOGPS currently supports only observations from the GPS constellation (L1/L2 frequency), with observation files in RINEX v2.XX format. LEOGPS also uses the precise ephemeris (.EPH) and clock bias and drift files (.CLK) provided by the University of Bern, Center for Orbit Determination in Europe (​CODE). As such, the coordinate frame used in the relative positioning is the International Terrestrial Reference Frame (ITRS) which is an Earth-Centered Earth-Fixed (ECEF) frame. Since the ephemeris files and RINEX v2 observations default to GPS Time, it is very important to also note that the time scale used in LEOGPS output files is GPS Time (as opposed to UTC).

This project also gives sincere appreciation and credit to the University of Bern, for their provision of the CODE FTP.

Documentation tree for LEOGPS is listed below.

.. toctree::
   :maxdepth: 1
   :caption: Getting Started
   
   docs_install.rst
   docs_firststep.rst
   docs_firstrun.rst

.. toctree::
   :maxdepth: 1
   :caption: Differential Navigation
   
   docs_cdgps1.rst
   docs_cdgps2.rst

.. toctree::
   :maxdepth: 1
   :caption: Advanced References
   
   docs_process.rst
   docs_functions.rst

For bugs, raise the issues in the `GitHub repository <https://github.com/sammmlow/LEOGPS/issues>`_. For collaborations, reach out to me (sammmlow@gmail.com). The project is licensed under the MIT license.

*Written By: Samuel Y. W. Low*

|linkedin| |orcid|
