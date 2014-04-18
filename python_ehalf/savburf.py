#!/usr/bin/python
# This script is the python version of the Savage and Burford model
# First, open parameter file and read the contents.  NOTE -- I changed the contents of the parameter file around.  Also, the param file needs to be renamed to param.py for importing variables.  
#param.py includes:
#xmin=-150
#xmax=150 
#int=1 
#Vo=50
#d=30
plot="contour"  # contour or lines
# Import math modules
import numpy as np
import math 
# Import plotting modules
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
# Import variables in param.py


import param
print "xmin xmax int Vo d <-- Forward Values"
print param.xmin, param.xmax, param.int, param.Vo, param.d

f = open('vel.txt','w')
listx = []
listVel = []
x=param.xmin 
while (x <= param.xmax):
     Vel=-((param.Vo/np.pi)*math.atan(x/param.d))
     print >> f, x, Vel
     listx.append(x)
     listVel.append(Vel)
     x = x + param.int

#Open GPS file 
g=np.loadtxt('data.py')
gx = g[:,0]
gVel = g[:,1]
gsig = g[:,2]

#Calculate chi2 for Forward model
VelC = -((param.Vo / np.pi) * np.arctan ([ gx/param.d ]))
chi = ((gVel - VelC)/ (gsig))**2

chi2 = sum(chi.T)
redchi = chi2/(len(gVel)-3)
print "Forward chi2 = ", chi2
print "Forward Reduced chi2 = ", redchi



#Invert data
b=np.matrix( [gVel] )
bsig=np.matrix( [gsig] )
bT=b.T
A = np.matrix(-((1 / np.pi) * np.arctan ([ gx/param.d ]))).T
AT = A.T
ATA = AT*A

# det(ATA) is singular -- therefore we need to do singular value decomposition to get an inverse
U,s,V = np.linalg.svd(ATA)
ATAinv = np.dot(np.dot(V.T,np.linalg.inv(np.diag(s))),U.T)

# Finally, solve for the best value of R
xx =  ATAinv * AT * bT


# Run Forward model with the inverse value
# Extract matrix value for use below
xxaa = np.matrix.item(xx[0:])

listxinv = []
listVelinv = []
x=param.xmin
while (x <= param.xmax):
     Velinv=-((xxaa/np.pi)*math.atan(x/param.d))
     listxinv.append(x)
     listVelinv.append(Velinv)
     x = x + param.int

# Calculate the Chi2 and reduced Chi2 for the inverse model
VelCinv = -((xxaa / np.pi) * np.arctan ([ gx/param.d ]))
chiinv = ((gVel - VelCinv)/ (gsig))**2

chi2inv = sum(chiinv.T)
redchiinv = chi2inv/(len(gVel)-3)

if plot == "lines":
 #Plot results 
 plt.plot(listx,listVel)
 plt.plot(listx,listVelinv)
 plt.errorbar(gx,gVel,gsig,fmt='r^',ms=10)

 plt.annotate( xxaa , xy =(2,1), xytext=(90,15))
 plt.annotate( "Inv R (green) =" , xy =(2,1), xytext=(30,15))
 plt.annotate( "Inv R chi2 =" , xy =(2,1), xytext=(30,10))
 plt.annotate( redchiinv, xy =(2,1), xytext=(100,10))
 
 plt.annotate( param.Vo , xy =(2,1), xytext=(100,5))
 plt.annotate( "Forward R (blue) =" , xy =(2,1), xytext=(30,5))
 plt.annotate( "Forward R chi2 =" , xy =(2,1), xytext=(30,0))
 plt.annotate( redchi, xy =(2,1), xytext=(100,0))
 
 plt.annotate( param.d , xy =(2,1), xytext=(90,-5))
 plt.annotate( "D =" , xy =(2,1), xytext=(30,-5))
 
 plt.show()


if plot == "contour":
# Perform a gridsearch to make a contour plot
 dmin=param.dmin
 dmax=param.dmax
 Vmin=param.Vmin
 Vmax=param.Vmax

 d=param.dmin
 gridredchi = np.array([V, d, chi])
 grc = []
 c = open('chi.py','w')

 while (d <= dmax):
  V=param.Vmin
  while (V <=  Vmax):
    gridVelC = -((V / np.pi) * np.arctan ([ gx/d ]))
    gridchi = ((gVel - gridVelC)/ (gsig))**2
    gridchisum = np.matrix.item(sum(gridchi.T))
    gridrchi= gridchisum/(len(gVel)-3)
    newrow =  [ V, d, gridrchi ]
    gridredchi = np.vstack([gridredchi, newrow])
    print >> c, V, d, gridrchi
    plt.scatter(V,d, c=gridrchi, marker='s',lw=0,  s=40, vmin=0, vmax=10)
    V = V + param.int
  d = d + param.int
 print gridredchi.shape
 #print gridredchi[:,0]
 #xspace=np.linspace(dmin,dmax,1)
 #yspace=np.linspace(Vmin,Vmax,1)
 #zspace=ml.griddata(gridredchi[:,0],gridredchi[:,1],gridredchi[:,2],xspace,yspace)
 #plt.contour(xspace,yspace,zspace,5)
 plt.scatter(xxaa,param.d, marker='*',color='white', s=100)
 plt.xlim([param.Vmin,param.Vmax])
 plt.ylim([param.dmin,param.dmax])
 plt.show()
