# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 10:09:42 2018

@author: micfe03
"""

import numpy as np
import imageio
from matplotlib import pyplot as plt
import chanpy as cp

if __name__ == '__main__':
    ngc = 10
    nrc = 24
    ncc = 28
    image = imageio.imread("4.1.07.tiff")
    image = image[:,:,0]
    plt.imshow(image,'gray')
    plt.title("Original Image")
    plt.colorbar()
    plt.show()
    CCBB = cp.Cos2ChannelBasis()
    CCBB.setParameters(ngc, 0, 1)
    CCBC = cp.Cos2ChannelBasis()
    CCBC.setParameters(ncc, 0, image.shape[1]-1)
    CCBR = cp.Cos2ChannelBasis()
    CCBR.setParameters(nrc, 0, image.shape[0]-1)
    image = 1.0/255.0*image
    CCFM = cp.CombinedChannelBasis()
    chBCCFM = [CCBB, CCBR, CCBC]
    CCFM.setParameters(chBCCFM)
    cCCFM = cp.ChannelVector(CCFM)
    cCCFM.addSample(image)
    cCCBB = cp.ChannelVector(CCBB)
    cCCBB.addSample(image[:nrc,:ncc]) # just to set parameters correctly
    cCCBB[:] = cCCFM.reshape(cCCBB.shape,order='F')
    fImage3 = cCCBB.decode(1)
    print('img3',fImage3.shape)
    plt.figure()
    plt.imshow(fImage3[:,:,0].squeeze(),'gray')
    plt.colorbar()    
    plt.title("Mode1 Image")
    plt.show()