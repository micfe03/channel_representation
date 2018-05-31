# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 14:58:27 2018

@author: micfe03
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec as gridspec 
import matplotlib
import chanpy as cp

if __name__ == '__main__':
    matplotlib.rcParams.update({'font.size': 20})
    matplotlib.rcParams.update({'font.family': 'serif'})
    matplotlib.rcParams.update({'text.usetex': True})

    nCCB = 5
    nBCB = 6
    stp = 16
    dN = 2
    rCMax = (nCCB-dN)
    rCMin = 0
    rBMax = (nBCB-dN)*stp
    rBMin = 0

    xi = np.arange(rBMin,rBMax+1)
    xo = np.arange(rBMin-stp,rBMax+stp+1)
    yi = np.arange(rCMin,rCMax+1./stp,1./stp)
    yo = np.arange(rCMin-1,rCMax+1+1./stp,1./stp)

    print("Cos2-channels with " + str(nCCB) + " coefficients from " + str(rCMin) + " to " + str(rCMax))

    CCB = cp.Cos2ChannelBasis()
    CCB.setParameters(nCCB, rCMin, rCMax)
    CCBo = cp.Cos2ChannelBasis()
    CCBo.setParameters(nCCB+dN, rCMin-1, rCMax+1)

    print("Cos2-channels with " + str(nBCB) + " coefficients from " + str(rBMin) + " to " + str(rBMax)) 

    BCB = cp.Cos2ChannelBasis()
    BCB.setParameters(nBCB, rBMin, rBMax)
    BCBo = cp.Cos2ChannelBasis()
    BCBo.setParameters(nBCB+dN, rBMin-stp, rBMax+stp)

    print("Cross product of channels")
    CBCB = cp.CombinedChannelBasis()
    chBV = [BCB, CCB]
    CBCB.setParameters(chBV)

    print("Encode and decode $y$ as cos2 channel")
    cCV = cp.ChannelVector(CCBo)
    cCV.addSample(yo)
#    print("Vector is " + str(cCV))
    cCVi = cp.ChannelVector(CCB)
    
    print("Encode and decode $x$ as cos2 channel")
    cBV = cp.ChannelVector(BCBo)
    cBV.addSample(xo)
#    print("Vector is " + str(cBV))
    cBVi = cp.ChannelVector(BCB)

    plt.close(1)
    fig = plt.figure(1) #, figsize=(15,8))
    gs = gridspec.GridSpec(2,2,height_ratios=[nCCB,1],width_ratios=[1,nBCB])
    gs.update(left=0, right=1, bottom=0, top=1, wspace=0, hspace=0)

    ax1 = plt.subplot(gs[0,0])
    ax1.plot(-cCV[:,1:-1],yo)
    ax1.set_ylim([yo[0],yo[-1]])
    ax1.axis('off')

    ax2 = plt.subplot(gs[1,1])
    ax2.plot(xo,-cBV[:,1:-1])
    ax2.set_xlim([xo[0],xo[-1]])
    ax2.axis('off')

    ax3 = plt.subplot(gs[0,1])
    f = (nCCB-dN)/(1.+np.exp(rBMax*2./stp-xi*4./stp))
    ax3.plot(xi,f)
#    cbCV = cp.ChannelVector(CBCB)
#    DM = np.array([xi.transpose(), f.transpose()]).transpose()
#    cbCV.addSample(DM)
    xv, yv = np.meshgrid(np.linspace(-stp/2.,rBMax+stp/2.,nBCB),np.linspace(-1./2.,rCMax+1./2.,nCCB))
#    ax3.scatter(xv.transpose().flatten(),yv.transpose().flatten(),cbCV)
    cBVi.addSample(xi)
    cCVi.addSample(f)
    cbCV=cCVi.transpose()@cBVi
    plt.scatter(xv,yv,4*cbCV**2)
    ax3.set_xlim([xo[0],xo[-1]])
    ax3.set_ylim([yo[0],yo[-1]])
    ax3.set_xticks(np.linspace(xi[0],xi[-1],nBCB-dN+1)) 
    ax3.set_xlabel('$x$',horizontalalignment='right',x=1.0,verticalalignment='bottom',y=1.0) # Force this empty !
    ax3.set_yticks(np.linspace(yi[0],yi[-1],nCCB-dN+1)) 
    ax3.set_ylabel('$f(x)$',horizontalalignment='right',y=1.0) # Force this empty !
    ax3.grid()
    fig.show()
    
#    plt.grid(clip_box=mtrans.Bbox([[rBMin,rCMin],[rBMax,rCMax]]))