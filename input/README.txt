LEOGPS will look for two RINEX files with the name conventions:

For LEO 1 - "<AAAA><DDD>0.<YY><Z>"
For LEO 2 - "<BBBB><DDD>0.<YY><Z>"

AAAA is the 4-letter ID of the spacecraft LEO1;
BBBB is the 4-letter ID of the spacecraft LEO2;
DDD is the day-of-year based on 1st RINEX observation;
YY is the 2-digit year number (i.e. 2020 --> 20);
Z is the RINEX observation extension
Z = 'O' for typical observation files
Z = 'D' for Hatanaka-compressed observation files

The naming convention follows standard RINEX conventions,
As followed by IGS and University of Bern CODE.

In <.../input/> sub-directory contains 2 default example files:
- GRCA2080.10D (GRACE-A RINEX file v2.11, on Y2010, 208th Day)
- GRCB2080.10D (GRACE-B RINEX file v2.11, on Y2010, 208th Day)

To use the GRACE examples, ensure they are found in <.../input/>.
Run LEOGPS, and you may validate results with the ground truths.
You can paste any other LEO RINEX observation files to run LEOGPS.

By default install, <.../config/config.txt> will already include
all parameters (LEO file names, and correct dates / times),
for running LEOGPS with the example LEO files.

Ensure also that you have installed the Hatanaka Python library with:

pip install Hatanaka

So that you can perform the Hatanaka (de)compression.