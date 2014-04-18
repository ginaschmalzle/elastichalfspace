Elastic Half Space Modeling of a Vertical Strike Slip fault done in Fortran and Python
=======================================================================================


BACKGROUND
The programs apply a simple "elastic half space" strike-slip  model to high precision 
GPS velocity data from a transect across the Carrizo segment of the San Andreas Fault (SAF), 
probably the most famous strike slip fault in the world. Strike-slip faults form when two tectonic 
plates slide past eachother.  These plates do not move underneath eachother nor are they moving apart.
In between major strike-slip earthquakes, the ground accrues stress from tectonic plates that are 
trying to slide past each other. The walls of the fault are stuck together (known as  "locked"), which 
causes the ground to deform. Both the fault rate (R) and the locking depth (d) are  important 
parameters in seismic hazard since it is thought that areas that are stuck will slip in the next big 
earthquake. Larger rates and larger slip areas would produce larger earthquakes.  


The surface velocity for a profile across a vertical strike-slip fault can be described by the function::

      vel = (R/pi) * arctan ( x / d)

where vel = calculated profile velocity field at the surface, pi = 3.14159..., x is distance
from the fault and d is the distance from the surface to a depth at which the fault is 'stuck'. 
Beneath this depth is a shear zone that moving at the fault rate, R

The equation is from Savage and Burford, 1973, JGR.  

To better understand what is observed at the surface (and estimated by this equation), imagine a fence 
built perpendicularly across a strike-slip fault.  When the fault is first built it is nice and staight, 
but over time it starts to deform and look kind of like an "S".  When the earthquake occurs the ground 
(and the fence) will snap, and the two sides of the fence will become straight again at some time after 
the earthquake, although displaced.  ehalf.f describes the "S" phase.  

FORTRAN PROGRAM DESCRIPTION
ehalf.f will read in parameters from the input file param.  Those parameters are used in calculating the 
velocity profile, which is output to vel.txt.  The program also takes the GPS velocity, position and 
uncertainty estimates contained in the file data.txt, and estimates how well the model fits the GPS data
by calculating the chi2 statistic which is given by::

	chi2 = SUM ((dataR -vel)/(sig))**2

where chi2 = chi2 statistic, dataR = GPS estimated rate for a position aling the profile, vel = model
calculated velocity, and sig = GPS velocity uncertainty.  SUM is the sum of all chi2 for each GPS datum.  
The reduced chi2 is also calculated and is given by:: 

	reduced chi2 = chi2 / (N-v-1) 

where N = number of data, v = number of variable parameters (in this case 2, R and d).  An ideal reduced chi2 
should equal 1.


Data and model profiles are plotted in the shell script display_me.  You will need General Mapping Tools (GMT) 5 
to display results.  To easily make changes in the fortran code, compile and plot results you can use the shell
script ./runme.

You can estimate the low misfit model -- ie., find the model with the smallest reduced chi2 by using the shell script wrapper.


PYTHON PROGRAM DESCRIPTION

The Python progam (savburf.py) provides the same functionality as the fortran program, except I include an alternative method 
to solve for the low misfit model given a defined locking depth. Here I invert the data using singular value decomposition.

