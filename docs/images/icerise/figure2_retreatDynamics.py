import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.patches as pat
from matplotlib.colors import LogNorm
from matplotlib import gridspec

from numpy import ma
from amrfile import io as amrio

#using shell command inside python
import subprocess

#scipy interpolation
#from scipy.interpolate import griddata
from scipy import interpolate

#1-need to make a loop over time steps for the two simulations
#2-for every time step, record the position of grounding line on the sides (less influence of the pinning point)
#3-also record the contribution to SLR
#4-plot the two as a function of time

#noIceRiseExperiment: increase the sea level rise from 5000 years
#seaLevelRise Experiment: from 5000 years too, and during 15000 years
#there is one record every 10 years

#save or not to save
save = 1

#for flotation
rhow = 1000.
rhoi = 900.

#define epsilon
eps = 0.1
step = 10
timestep = 10

plotschoof = 0
grav=9.8
A=9.4608e-18
C=24131.3
n=3.
m=0.33333333

def readfile(fnam):

	amrID = amrio.load(fnam)
        ztopcomp = "Z_surface"
        zbottomcomp = "Z_bottom"
        bedcomp = "Z_base"
        thkcomp = "thickness"
	xvelcomp = "xVel"
	yvelcomp = "yVel"
	ccomp = "basal_friction"
	level = 2 
	#first, work out the do	main corners
	lo,hi = amrio.queryDomainCorners(amrID, level)
	order = 1 # interpolation order, 0 for piecewise constant, 1 for linear
        x,y,ztop = amrio.readBox2D(amrID, level, lo, hi, ztopcomp, order)
        x,y,zbottom = amrio.readBox2D(amrID, level, lo, hi, zbottomcomp, order)
        x,y,bed = amrio.readBox2D(amrID, level, lo, hi, bedcomp, order)
        x,y,thk = amrio.readBox2D(amrID, level, lo, hi, thkcomp, order)
	x,y,xvel = amrio.readBox2D(amrID, level, lo, hi, xvelcomp, order)
        x,y,yvel = amrio.readBox2D(amrID, level, lo, hi, yvelcomp, order)
	x,y,c = amrio.readBox2D(amrID, level, lo, hi, ccomp, order)
	amrio.free(amrID)

        x = x/1000.
        y = y/1000.

	return x,y,ztop,zbottom,bed,thk,xvel,yvel,c

basis = "/home/lfavier/Programs/bisicles/BISICLES/examples/IceRiseSmallStudy/"

## firstly for the grounding line dynamics ##

#first for the noIceRiseExperiment
simulation = "noIceRiseExperiment/"

#create the list of the plot files, using an external program
subprocess.call(["./makelist.sh",simulation])

#loops over the listplotfiles text file, and take every line as an argument
listfiles = np.genfromtxt(basis+simulation+"listplotfiles",dtype='str')

#prepare the grounding line position vector to fill afterwards
glx = np.arange(0,len(listfiles),step)
#also contruct the time vector
time = np.arange(0,len(listfiles)*10,10*step)/1000.
#also contruct the flux vector
fluxgl = np.arange(0,len(listfiles),step)
#for the schoof flux
fluxsc = np.arange(0,len(listfiles),step)

cpt = 0
for i in range(0,len(listfiles),step):

	print listfiles[i]

	x,y,ztop,zbottom,bed,thk,xvel,yvel,c = readfile(basis+simulation+listfiles[i])
  
	#position of the grounding line, relies on the difference between the bedrock and the bottom surface
	difftoflot = zbottom-bed
	difftoflot = difftoflot[-1,:]
	condflot = difftoflot>eps
	posx = x[condflot==True]
	glx[cpt] = posx[0]

        #thickness along the flowline away from the pinning point
        thkfl = thk[-1,:]

        #take only the far flowline out of in the influence of the ice rise
        #velocity compute
        normvel = np.power(xvel*xvel+yvel*yvel,0.5)
        normvelfl = normvel[-1,:]

        #flux at the grounding line
        fluxgl[cpt] = normvelfl[x==posx[0]]*thkfl[x==posx[0]]

	#schoof flux at the grounding line
        fluxsc[cpt] = ((A*(rhoi*grav)**(n+1)*(1-rhoi/rhow)**n)/(4**n*C))**(1/(m+1))*(thkfl[x==posx[0]])**((m+n+3)/(m+1))

	cpt = cpt+1
	

#to calculate the retreat speed of the last grounded point
retreatspeedNoRise = (glx[1:]-glx[0:-1])/(timestep*step)
timeNoRise = time[1:]

#take the ymax for plotting purpose later on
width1 = max(y)-min(y)

#start the plot
fig = plt.figure(1,figsize=(7,4))

gs = gridspec.GridSpec(2,1,height_ratios=[7,2.5])

ax1 = fig.add_subplot(gs[0])

plt11 = ax1.plot(time,glx,'b',linewidth=2.0,label="without ice rise")

#second for the basic experiment
simulation = "seaLevelRise/"

#create the list of the plot files, using an external program
subprocess.call(["./makelist.sh",simulation])

#loops over the listplotfiles text file, and take every line as an argument
listfiles = np.genfromtxt(basis+simulation+"listplotfiles",dtype='str')
#prepare the grounding line position vector to fill afterwards
glx = np.arange(0,len(listfiles),step)
#also contruct the time vector
time = np.arange(0,len(listfiles)*10,10*step)/1000.
#also contruct the flux vector
fluxgl = np.arange(0,len(listfiles),step)
#for the schoof flux
fluxsc = np.arange(0,len(listfiles),step)

cpt = 0
for i in range(0,len(listfiles),step):

	print listfiles[i]

	x,y,ztop,zbottom,bed,thk,xvel,yvel,c = readfile(basis+simulation+listfiles[i])

	#position of the grounding line, relies on the difference between the bedrock and the bottom surface
	difftoflot = zbottom-bed
	difftoflot = difftoflot[-1,:]
	condflot = difftoflot>eps
	posx = x[condflot==True]
	glx[cpt] = posx[0]

        #thickness along the flowline away from the pinning point
        thkfl = thk[-1,:]

        #take only the far flowline out of in the influence of the ice rise
        #velocity compute
        normvel = np.power(xvel*xvel+yvel*yvel,0.5)
        normvelfl = normvel[-1,:]

        #flux at the grounding line
        fluxgl[cpt] = normvelfl[x==posx[0]]*thkfl[x==posx[0]]

        #schoof flux at the grounding line
        fluxsc[cpt] = ((A*(rhoi*grav)**(n+1)*(1-rhoi/rhow)**n)/(4**n*C))**(1/(m+1))*(thkfl[x==posx[0]])**((m+n+3)/(m+1))

        cpt = cpt+1


#to calculate the retreat speed of the last grounded point
retreatspeedIceRise = (glx[1:]-glx[0:-1])/(timestep*step)
timeIceRise = time[1:]

#take the ymax for plotting purpose later on
width2 = max(y)-min(y)

plt12 = ax1.plot(time,glx,'k',linewidth=2.0,label="with ice rise")
#ax1.set_title('Position of the last grounded point')
ax1.set_xlabel('Time (ka)')
ax1.set_ylabel('Position (km)')
ax1.set_xlim(5,30)
ax1.set_ylim(300,750)
ax1.set_yticks([300,400,500,600,700])
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top')

plt.legend(loc="upper right")

ax1.grid()

ax2 = fig.add_subplot(gs[1])

#sea level rise plot
slr=10*(time-5)
slr[slr>150] = 150
plt21 = ax2.plot(time,slr,'k',linewidth=2.0)

#ax2.set_title('Contribution to SLR for a 100 km wide glacier')
ax2.set_xlabel('Time (ka)')
ax2.set_ylabel('SLR (m)')
ax2.set_xlim(5,30)
ax2.set_ylim(-10,160)
ax2.set_yticks([0,50,100,150])
ax2.grid()

plt.tight_layout()

#saving stuff
#plt.savefig("figure2.pdf")
plt.savefig("figure2.png")

plt.show()


