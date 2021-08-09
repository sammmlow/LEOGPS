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

Carrier Double Differencing
===========================

We can actually go one step further beyond single differencing, and that is to difference across two reference GPS satellites. This step is called double differencing, and it is the backbone algorithm used in LEOGPS. This removes the relative receiver clock bias estimation errors, and this section will detail the algorithm and the results below.

.. figure:: /_figures/leogps_2diff_fig1.png
   :align: center
   
   Figure C2.1: Carrier phase double differencing equation setup.

The key idea behind double differencing, is to create a differential between two single difference equations. For example, between two LEO units A and B, if some GPS satellite P was taken as the reference emitter satellite in the single difference equation, then a second single difference equation would be taken with respect to a different GPS satellite Q as the reference emitter. In both single difference equations, it makes no difference to the relative receiver clock bias estimation which reference GPS satellite is used, as these are the errors specific only to the receiver-side in the relative sense.

Thus, taking two single difference measurements, one can observe that the relative receiver clock bias estimation error has now become a common error source in the double differenced equation, which can be cancelled out (both relative receiver clock bias terms in the red box in Figure C2.1 are actually the same and thus they cancel).

.. figure:: /_figures/leogps_2diff_fig2.png
   :align: center
   
   Figure C2.2: Extraction of the true baseline vector via double differencing.

In a similar fashion to the single differencing algorithm, the double differencing algorithm requires the estimation of the unit direction vector from the LEO A to reference GPS satellite P, and from the LEO B to the reference GPS satellite Q. Once again, with the double-differenced relative range, and on the assumption that the signal path pairs from P to A and P to B are parallel, and the path pairs from Q to A and Q to B are also parallel, one can de-construct the baseline vector AB since each single difference observation forms the base of a right-angled triangle.

Since two reference satellites are now used in the creation of a double differenced equation, for N number of common GPS satellites in view, one can expect N-1 number of double differenced equations. Solving the entire system of equations leads to the baseline vector solution, with the error norm plotted below (error taken against a ground truth for two LEO units at a 100km baseline separation). It is observed that the final elimination of receiver clock biases creates an almost ten fold increase in the navigation accuracy.

.. figure:: /_figures/leogps_2diff_plot.png
   :align: center
   
   Figure C2.3: Relative position error norm of a double differenced 100km LEO baseline.

In its essence, the single differencing paradigm changed the structure of the navigation problem, eliminating most nuisance errors, while the double differencing paradigm went further and essentially eliminated the receiver clock bias at the expense of one equation in the entire system of equations, while also redefining the ambiguity term into a double-differential ambiguity.

We can actually go one step further and perform triple-differencing, that is, performing a differencing measurement between two double-differences. This eliminates the ambiguity entirely, but it amplifies the random errors by another factor of root-2, and it further reduces the degrees of freedom by one, from N-1 to N-2, N being the number of GPS satellites.

At this point though, it becomes arguable that estimating the ambiguity with the C/A range (float ambiguity resolution) on a double-differencing measurement would give a better accuracy than performing triple differencing. This is the reason why LEOGPS does not use triple-differencing.
