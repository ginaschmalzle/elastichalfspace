	program ehalf
* Program ehalf was written by Gina Schmalzle on December 12, 2013.
* 
* Calculate model for a profile at specified spatial intervals
* Model variables:
* 	R = fault rate
* 	Vel = Calculated Surface velocity
* 	x = Distance from the fault
* 	xmin, xmax, dint = min and max distance of the profile and the 
* 	interval at which a velocity is calculated
*	 d = depth of the elastic layer
	
	double precision R, Vel, x, d
	double precision xmin, xmax, xint
	double precision sig, dataR, datax
	double precision chi
	integer N

* Input values from a parameters file param
	open(13, file='param')
	rewind (13)
	    READ (13,*,err=50)
	    READ (13,*,err=50)
	    READ (13,*,err=50) xmin, xmax, xint
	    READ (13,*,err=50)
	    READ (13,*,err=50) R
	    READ (13,*,err=50)
	    READ (13,*,err=50) d 
	close(13)
	goto 51
50	print *, "Error reading param"
	goto 102
51	CONTINUE
* Compute the Surface Velocity and write to file
	pi = 3.141592654
	open(14, file='vel.txt')
	rewind (14)
	x = xmin 
60	  if ( x .le. xmax) then 
	         Vel = 0
	         Vel = -((R/pi)*ATAN(x/d))
	         WRITE (14,*) x, Vel
	         x = x + xint
	         goto 60      
	  endif
	close(14)
	
* Calculate chi2
* The second part of this code generates some statistics on how well the model 
* GPS data.  You will need a file called data.txt that containes uncertainty in 
* mm/yr, Distance in km and rate in mm/yr.
* The model is re-run for the given distance (x), and the chi2 is calculated with 
* the parameters R and d given in the input file. 
* Open Data file and model output file
*	chi = chi2 statistic
*	dataR = GPS observed rate in mm/yr
*	datax = GPS perpendicular distance from fault in km
*
* The Reduced chi2 is also calculated and is chi2 / (N-rho-1), which N = number of data and
* rho = the number of variables (R and d)  which = 2  

	chi = 0
	N = 0
	open(15,file='data.txt')
	rewind (15)
70	  READ (15,*,err=100, end=101) sig, datax, dataR
	      Vel = 0
	         Vel = -((R/pi)*ATAN(datax/d))
	         chi = chi + ((dataR -Vel)/(sig))**2
	         N = N + 1 
	      goto 70
100	print *, "Error reading data.txt"
	goto 102
101	print *, "R(mm/yr), d(km), Red chi2"
	print *, R, d, chi /(N-3) 
102	CONTINUE
	close (15)
	END
