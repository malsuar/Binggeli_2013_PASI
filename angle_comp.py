#!/usr/bin/env python
from scipy import *
from pylab import *
import numpy as np
from sys import stderr
from pylab import * 
from sys import exit
from sys import *
from matplotlib.pylab import  hist, show
from matplotlib import *
plt=matplotlib.pyplot

def DIS(pos1, pos2):

    b = pos1-pos2
    b=np.sqrt(np.dot(b,b))
    return b

n=3                                             #Number of slices 
pdot=range(0,1)                             #Histogram cutoff 
frecang=[]
frecdis=[]                                         #histogram counts
for z in range (0,n,1):

    param=np.loadtxt('./data/datos%d.txt'%(z), delimiter=',')
    if "#" in param:continue
    BDMID=param[:,0]     # Halo index
    AXIS1X=param[:,1]    # x-component of major axis
    AXIS1Y=param[:,2]    # y-component of major axis
    AXIS1Z=param[:,3]    # z-component of major axis      
    MTHALO=param[:,4]    # 1/h Msun - halo mass - mass of all particles within Rvir
    MVHALO=param[:,5]    # 1/h Msun - halo mass - mass of bound particles within Rvir
    XXHALO=param[:,6]    # 1/h Mpc  - (comoving) position, x-component
    YYHALO=param[:,7]    # 1/h Mpc  - (comoving) position, y-component
    ZZHALO=param[:,8]    # 1/h Mp   - (comoving) position, z-component
    NPHALO=param[:,9]    # number of particles in halo

    angle=open('./results/slice%d.dat'%(z),'w')
    for i in range (0,len(param),1):
        for j in range (i+1,len(param),1):
            hvect1=array([AXIS1X[i], AXIS1Y[i], AXIS1Z[i]])
            hvect2=array([AXIS1X[j], AXIS1Y[j], AXIS1Z[j]])
            sdot=np.dot(hvect1,hvect2)**2
            pos1=array([XXHALO[i], YYHALO[i], ZZHALO[i]])
            pos2=array([XXHALO[j], YYHALO[j], ZZHALO[j]])
            cart=DIS(pos1,pos2) # UNIDADES? -  Distancia entre cumulos 
            print  pos1, pos2, cart
            angle.write('%.1e\n'%(sdot))
            #print hvect1, hvect2, sdot, distance
            frecang.append(sdot)
            frecdis.append(cart)
            
#HISTOGRAMA

    hist(frecang,30,(0,1), color= "w",normed = 1, histtype='stepfilled')
    ax=plt.axes()
    ax.set_title(u"squared dot at z=%d"%(z))
    ax.set_xlabel(r"square dot angle")
    ax.set_ylabel(r"Frecuency")
#   plt.text(7, 70, r'$\mu_t=%.4s$'%(media))
    plt.grid(True)
    plt.savefig('./graph/distt%d.eps'%(z))

    angle.close()
