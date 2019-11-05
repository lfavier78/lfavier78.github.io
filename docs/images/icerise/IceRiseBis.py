#! /usr/bin/python

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as pat
from matplotlib import cm
from matplotlib.colors import LogNorm
from matplotlib import gridspec

from mpl_toolkits.mplot3d import Axes3D

from numpy import ma
from amrfile import io as amrio

#using shell command inside python
import subprocess

#scipy interpolation
#from scipy.interpolate import griddata
from scipy import interpolate

#define epsilon
eps = 0.1
step = 10 #means every step*10 years 

rhoi = 900
rhow = 1000

# range limits in km
xmin = 300 
xmax = 800
ymin = -50
ymax = 50
zmin = -500
zmax = 1500

#new grid
xnew = np.arange(xmin,xmax,1)
ynew = np.arange(ymin,ymax,1)
xnewm,ynewm = np.meshgrid(xnew,ynew)

#for velocity colorbar
vmin=1.
vmax=600.

def readfile(fnam,lev):

        amrID = amrio.load(fnam)
        ztopcomp = "Z_surface"
        zbottomcomp = "Z_bottom"
        bedcomp = "Z_base"
        thkcomp = "thickness"
        xvelcomp = "xVel"
        yvelcomp = "yVel"
        ccomp = "basal_friction"
        level = lev
        #first, work out the do main corners
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

        #reconstruction du domaine entier
        yneg = y-max(y)
        yall = np.concatenate((yneg,y),axis=0)

        ztopneg = np.flipud(ztop)
        ztopall = np.concatenate((ztopneg,ztop),axis=0)

        zbottomneg = np.flipud(zbottom)
        zbottomall = np.concatenate((zbottomneg,zbottom),axis=0)

        bedneg = np.flipud(bed)
        bedall = np.concatenate((bedneg,bed),axis=0)

        thkneg = np.flipud(thk)
        thkall = np.concatenate((thkneg,thk),axis=0)

        xvelneg = np.flipud(xvel)
        xvelall = np.concatenate((xvelneg,xvel),axis=0)

        yvelneg = np.flipud(-yvel)
        yvelall = np.concatenate((yvelneg,yvel),axis=0)

        cneg = np.flipud(c)
        call = np.concatenate((cneg,c),axis=0)

        return x,yall,ztopall,zbottomall,bedall,thkall,xvelall,yvelall,call

basis = "/home/lfavier/Programs/bisicles/BISICLES/examples/IceRiseSmallStudy/"
simulation = "seaLevelRise/"

# create the list of the plot files, using an external program
subprocess.call(["./makelist.sh",simulation])

#loops over the listplotfiles text file, and take every line as an argument
listfiles = np.genfromtxt(basis+simulation+"listplotfiles",dtype='str')

#colormap
cmap1 = plt.get_cmap('gray')
cmap2 = plt.get_cmap('jet')

cpt = 0
for i in range(0,len(listfiles),step):

	#add time in the plots
	tt = i*10
	if (tt<10):
		disp = '0000'+str(tt)
	elif (tt<100):
		disp = '000'+str(tt)
	elif (tt<1000):
		disp = '00'+str(tt)
	elif (tt<10000):
		disp = '0'+str(tt)
	else:
		disp = str(tt)

	if (tt<25000 or tt>30000):
	#if (tt!=10000):
		continue

	print listfiles[i]

	x,y,ztop,zbottom,bed,thk,xvel,yvel,c = readfile(basis+simulation+listfiles[i],2)

	xm,ym = np.meshgrid(x,y)

	#interpolate on final grid
	f = interpolate.interp2d(x,y,ztop,kind='linear')
	ztopn = f(xnew,ynew)
	f = interpolate.interp2d(x,y,thk,kind='linear')
	thkn = f(xnew,ynew)
	f = interpolate.interp2d(x,y,bed,kind='linear')
        bedn = f(xnew,ynew)

	#topographic high summit
	bedmoins = bed[len(y)/2-1,:]
	bedmoins[x<xmin] = -9999.
	bedmoins[x>xmax] = -9999.
	maxlocal = np.max(bedmoins)
	maxind = bedmoins==np.max(bedmoins)
	summitx = x[maxind]

	#velocity compute
	normvel = np.power(xvel*xvel+yvel*yvel,0.5)
	f = interpolate.interp2d(x,y,normvel,kind='linear')
	normveln = f(xnew,ynew)
	#logarithmic velocity
	normvelnlog = np.log10(normveln)
	vminlog = np.log10(vmin)
	vmaxlog = np.log10(vmax)

	#for contour plots
	bedc = [-25] #np.arange(-40,-20,10)
	ztopc = np.arange(75.,225.,25.)

	#colors for the surface
	a = 1/(vmaxlog-vminlog)
	b = -a*vminlog
	colors = a*normvelnlog+b
	#print colors	
	colors[colors>1] = 1
	colors[colors<0] = 0

	#color for the bottom surface
	condfloat = (bedn<-thkn*rhoi/rhow)
	colorsb = np.empty(condfloat.shape)
	colorsb[condfloat==True] = 0.9
	colorsb[condfloat==False] = 0.5

	fig = plt.figure(i,figsize=(10,5))

	gs = gridspec.GridSpec(1,2,width_ratios=[20,1])

	ax1 = fig.add_subplot(gs[0],projection='3d')
	ax2 = fig.add_subplot(gs[1])

	plt11 = ax1.plot_surface(xnewm,ynewm,ztopn,linewidth=0,rstride=1,cstride=1,facecolors=cm.jet(colors))
	ax1.set_xlim(xmin,xmax)
	ax1.set_ylim(ymin,ymax)
	ax1.set_zlim(zmin,zmax)
	ax1.set_xlabel('x (km)')
	ax1.set_ylabel('y (km)')
	ax1.set_zlabel('z (m)')
	#ax1.set_yticks([-20,-10,0,10,20])
	#ax1.xaxis.tick_top()
	#ax1.xaxis.set_label_position("top")
	ax1.set_title('               after '+np.str(tt)+' years')
	plt12 = ax1.plot_surface(xnewm,ynewm,ztopn-thkn,linewidth=0,rstride=1,cstride=1,facecolors=cm.gray(colorsb))

	#ax3 colorbar construction
	xc = np.arange(2)
	yc = np.arange(np.log10(vmin),np.log10(vmax),0.01)
	barc = np.zeros((len(yc),2))
	barc[:,0] = yc
	barc[:,1] = yc

	plt21 = ax2.pcolormesh(xc,yc,barc,cmap=cmap2)
	ax2.set_xlim(min(xc),max(xc))
	ax2.set_ylim(min(yc),max(yc))
	ax2.set_xticks([])
	ax2.set_yticks([0.,1.,2.,np.log10(500)])
	ax2.set_yticklabels(['1','10','100','500'])
	ax2.yaxis.tick_right()
	ax2.set_title('u (m a$^{-1}$)')

	plt.tight_layout()

	#saving plot
	plt.savefig("movieicerise"+disp+".png")
	#subprocess.call(["convert","figure3.png", "figure3.pdf"])

	#plt.show()

