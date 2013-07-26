#!/usr/bin/env python
from scipy import *
from pylab import *
import numpy as np
from subprocess import call
from sys import stderr
from pylab import * 
from sys import *
from matplotlib.pylab import  hist, show
from matplotlib import *
plt=matplotlib.pyplot

def DIS(pos1, pos2):

    b = pos1-pos2
    b=np.sqrt(np.dot(b,b))
    return b

def NORM(v1,v2):
    return np.sqrt(v1**2+v2**2)

n=40                                             #Number of slices 
pdot=range(0,1)                                    #Histogram cutoff 
frecang=[]
frecdis=[]                                         #histogram counts
for z in range (0,n,40):

    param=np.loadtxt('./data/MDR1_BDMW_snap_85_z_%.1f_%.1f.csv'%(z,z+40), delimiter=',', comments="#", skiprows=37)
#    if "#" in param:continue
    BDMID=param[:,0]     # Halo index
    SNAPNUM=param[:,1]
    NINCAT=param[:,2]
    HOSTFLAG=param[:,3]
    XXHALO=param[:,4]     # 1/h Mpc  - (comoving) position, x-component
    YYHALO=param[:,5]     # 1/h Mpc  - (comoving) position, y-component
    ZZHALO=param[:,6]     # 1/h Mp   - (comoving) position, z-component
    VXHALO=param[:,7]     # 1/h Mpc  - (comoving) position, x-component
    VYHALO=param[:,8]     # 1/h Mpc  - (comoving) position, y-component
    VZHALO=param[:,9]     # 1/h Mp   - (comoving) position, z-component
    NPHALO=param[:,10]    # number of particles in halo
    MVHALO=param[:,11]    # 1/h Msun - halo mass - mass of bound particles within Rvir
    MTHALO=param[:,12]    # 1/h Msun - halo mass - mass of all particles within Rvir 
    AXIS1X=param[:,23]    # x-component of major axis
    AXIS1Y=param[:,24]    # y-component of major axis
    AXIS1Z=param[:,25]    # z-component of major axis      

    angle=open('./results/slice%d.dat'%(z),'w')
    for i in range (0,len(param),1):
        for j in range (i+1,len(param),1):
#3D
            hvect1=array([AXIS1X[i], AXIS1Y[i], AXIS1Z[i]])
            hvect2=array([AXIS1X[j], AXIS1Y[j], AXIS1Z[j]])
            sdot=np.dot(hvect1,hvect2)**2
            pos1=array([XXHALO[i], YYHALO[i], ZZHALO[i]])
            pos2=array([XXHALO[j], YYHALO[j], ZZHALO[j]])
            vec_cm=pos1-pos2
            sdot_cm_pos=array([np.dot(vec_cm,hvect1),np.dot(vec_cm,hvect2)])
            delta=array([500-XXHALO[i],500-YYHALO[i],500-ZZHALO[i]])
            pos1=array([500,500,500])
            pos2=(pos2+delta)%1000
            cart=DIS(pos1,pos2) # UNIDADES? -  Distancia entre cumulos 
#2D
            hvect1_2D=array([AXIS1X[i],AXIS1Y[i]])/NORM(AXIS1X[i],AXIS1Y[i])
            hvect2_2D=array([AXIS1X[j],AXIS1Y[j]])/NORM(AXIS1X[j],AXIS1Y[j])

            h1_dot_h2= (np.dot(hvect1_2D,hvect2_2D))**2 #(AXIS1X[i]*AXIS1X[j]+ AXIS1X[i]*AXIS1X[j])**2
            pos1_2D=array([XXHALO[i], YYHALO[i]])
            pos2_2D=array([XXHALO[j], YYHALO[j]])
            cart_2D=DIS(pos1_2D,pos2_2D)

            delta_2D=array([500-XXHALO[i],500-YYHALO[i]])
            pos1_2D=array([500,500])
            pos2_2D=(pos2_2D+delta_2D)%1000
            cart_2D=DIS(pos1_2D,pos2_2D)

            #print  pos1, pos2, cart
            angle.write('%0.2e %e %0.2e %e \n'%(cart,sdot,cart_2D,h1_dot_h2))
            #print hvect1, hvect2, sdot, distance
            frecang.append(sdot)
            frecdis.append(cart)
            
#HISTOGRAMAS
    #angles
    hist(frecang,30,(0,1), color= "w",normed = 1, histtype='stepfilled')
    ax=plt.axes()
    ax.set_title(u"squared dot at z=%d"%(z))
    ax.set_xlabel(r"square dot angle")
    ax.set_ylabel(r"Frecuency")
#   plt.text(7, 70, r'$\mu_t=%.4s$'%(media))
    plt.grid(True)
    plt.savefig('./graph/distt%d.eps'%(z))
    angle.close()

 
"""
#Distances
    hist(frecdis,30,(0,100), color= "w",normed = 1, histtype='stepfilled')
    ax=plt.axes()
    ax.set_title(u"squared dot at z=%d"%(z))
    ax.set_xlabel(r"square dot angle")
    ax.set_ylabel(r"Frecuency")
#   plt.text(7, 70, r'$\mu_t=%.4s$'%(media))
    plt.grid(True)
    plt.savefig('./graph/distt%d.eps'%(z))

"""


