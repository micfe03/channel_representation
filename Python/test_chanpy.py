import numpy as np
import numpy.random as rnd
from scipy import misc, ndimage
from matplotlib import pyplot as plt
import chanpy as cp

if __name__ == '__main__':
    nMax = 30
    nCCB = int(round(rnd.rand() * (nMax-2)) + 3)
    nBCB = int(round(rnd.rand() * (nMax-2)) + 3)
    rCMax = rnd.rand()
    rCMin = -rnd.rand()
    rBMax = rnd.rand()
    rBMin = -rnd.rand()

    x = rCMin+(rCMax-rCMin)*rnd.rand(2,2)
    y = rBMin+(rBMax-rBMin)*rnd.rand(2,2)
    r = np.array([[1.0, 0.9],[0.8, 0.7]])
    xr = np.array([x.transpose(), r.transpose()]).transpose()
    yr = np.array([y.transpose(), r.transpose()]).transpose()
    xy = np.array([x.transpose(), y.transpose()]).transpose()
    xyr = np.array([x.transpose(), y.transpose(), r.transpose()]).transpose()

    print("chanpy test")

    print("Cos2-channels with " + str(nCCB) + " coefficients from " + str(rCMin) + " to " + str(rCMax))

    CCB = cp.Cos2ChannelBasis()
    CCB.setParameters(nCCB, rCMin, rCMax)

    print("Cos2-channels with " + str(nBCB) + " coefficients from " + str(rBMin) + " to " + str(rBMax)) 
    BCB = cp.Cos2ChannelBasis()
    BCB.setParameters(nBCB, rBMin, rBMax)

    print("Cross product of channels")
    CBCB = cp.CombinedChannelBasis()
    chBV = [CCB, BCB]
    CBCB.setParameters(chBV)

    print("Channel vectors (non-sparse)")

    print("Encode and decode " + str(x[0:1,0:1]) + " as cos2 channel")
    cCV = cp.ChannelVector(CCB)
    cCV.addSample(x[0:1,0:1])
    print("Vector is " + str(cCV))
    print("Decoding gives " + str(cCV.decode()))

    print("Add " + str(xr[0:1,1:2]) + " as cos2 channel and decode")
    cCV.addSample(xr[0:1,1:2])
    print("Vector is " + str(cCV))
    print("Decoding gives " + str(cCV.decode(2)))

    cCV.fill(0.)

    print("Encode and decode " + str(xr[0:1,0:2]) + " as cos2 channel")
    cCV.addSample(xr[0:1,0:2])
    print("Vector is " + str(cCV))
    print("Decoding gives " + str(cCV.decode()))

    cCVsum = cp.ChannelVector(CCB)
    cCVsum[:] = cCV[0:1,:]+cCV[1:2,:]
    print("Average the two channels gives " + str(cCVsum))
    print("Decoding gives " + str(cCVsum.decode(2)))

    print("Encode and decode " + str(y[0:1,0:1]) + " as cos2 channel")
    bCV = cp.ChannelVector(BCB)
    bCV.addSample(y[0:1,0:1])
    print("Vector is " + str(bCV))
    print("Decoding gives " + str(bCV.decode()))

    print("Add " + str(yr[0:1,1:2]) + " as cos2 channel and decode")
    bCV.addSample(yr[0:1,1:2])
    print("Vector is " + str(bCV))
    print("Decoding gives " + str(bCV.decode(2)))

    bCV.fill(0.)

    print("Encode and decode " + str(yr[0:1,0:2]) + " as cos2 channel")
    bCV.addSample(yr[0:1,0:2])
    print("Vector is " + str(bCV))
    print("Decoding gives " + str(bCV.decode()))

    bCVsum = cp.ChannelVector(BCB)
    bCVsum[:] = bCV[0:1,:]+bCV[1:2,:]
    print("Average the two channels gives " + str(bCVsum))
    print("Decoding gives " + str(bCVsum.decode(2)))

    print("Encode and decode " + str(xy[0:1,0:1]) + " as combined channels")
    cbCV = cp.ChannelVector(CBCB)
    cbCV.addSample(xy[0:1,0:1])
    print("Vector is " + str(cbCV))
    print("Decoding gives " + str(cbCV.decode()) + " compared to " + str(xyr[0:1,0:1]))

    print("Add " + str(xyr[0:1,1:2]) + " as combined channel and decode")
    cbCV.addSample(xyr[0:1,1:2])
    print("Vector is " + str(cbCV))
    print("Decoding gives " + str(cbCV.decode(2)) + " compared to " + str(xyr[0:1,0:2]))

    cbCV.fill(0.)
    
    print("Encode and decode " + str(xyr[0:1,0:2]) + " as combined channel")
    cbCV.addSample(xyr[0:1,0:2])
    print("Vector is " + str(cbCV))
    print("Decoding gives " + str(cbCV.decode()) + " compared to " + str(xyr[0:1,0:2]))

    cbCVsum = cp.ChannelVector(CBCB)
    cbCVsum[:] = cbCV[0:1,:]+cbCV[1:2,:]
    print("Average the two channels gives " + str(cbCVsum))
    print("Decoding gives " + str(cbCVsum.decode(2)) + " compared to " + str(xyr[0:1,0:2]))

    image = misc.imread("house_orig.png",flatten=True)
    image = image[0:248,:]
    plt.imshow(image,'gray')
    plt.title("Original Image")
    plt.colorbar()
    CCBB = cp.Cos2ChannelBasis()
    CCBB.setParameters(13, 0., 1.)
    cCVB = cp.ChannelVector(CCBB)
    cCVB.addSample(1./255.*image)
    tmpImg = ndimage.filters.gaussian_filter(cCVB.channelImage(), 1.5)
    cCVB.channelImage()[:] = tmpImg
    fImage2 = cCVB.decode(2)
#    plt.imshow(cv2.convertScaleAbs(fImage2[:,:,0].squeeze()))
    print(fImage2.shape)
    plt.figure()
    plt.imshow(fImage2[:,:,0].squeeze(),'gray')
    plt.colorbar()
    plt.title("Mode0 Image")
    plt.figure()
    plt.imshow(fImage2[:,:,2].squeeze(),'gray')
    plt.colorbar()    
    plt.title("Mode1 Image")

    CCBC = cp.Cos2ChannelBasis()
    CCBC.setParameters(14, 0, image.shape[1]-1)
    CCBR = cp.Cos2ChannelBasis()
    CCBR.setParameters(12, 0, image.shape[0]-1)
    CCBB.setParameters(10, 0, 1)
    image = 1.0/255.0*image
    print(image.shape, image[100,250])
    CCFM = cp.CombinedChannelBasis()
    chBCCFM = [CCBB, CCBR, CCBC]
    CCFM.setParameters(chBCCFM)
    cCCFM = cp.ChannelVector(CCFM)
    cCCFM.addSample(image)
    print("Decoding gives " + str(cCCFM.decode(1)))
    plt.show()

    print("Histogram with " + str(nCCB) + " coefficients from " + str(rCMin) + " to " + str(rCMax))

    CCB = cp.HistogramBasis()
    CCB.setParameters(nCCB, rCMin, rCMax)

    print("Histogram with " + str(nBCB) + " coefficients from " + str(rBMin) + " to " + str(rBMax)) 
    BCB = cp.HistogramBasis()
    BCB.setParameters(nBCB, rBMin, rBMax)

    print("Cross product of channels")
    CBCB = cp.CombinedChannelBasis()
    chBV = [CCB, BCB]
    CBCB.setParameters(chBV)

    print("Channel vectors (non-sparse)")

    print("Encode and decode " + str(x[0:1,0:1]) + " as histogram")
    cCV = cp.ChannelVector(CCB)
    cCV.addSample(x[0:1,0:1])
    print("Vector is " + str(cCV))
    print("Decoding gives " + str(cCV.decode()))

    print("Add " + str(xr[0:1,1:2]) + " as histogram and decode")
    cCV.addSample(xr[0:1,1:2])
    print("Vector is " + str(cCV))
    print("Decoding gives " + str(cCV.decode(2)))

    cCV.fill(0.)

    print("Encode and decode " + str(xr[0:1,0:2]) + " as histogram")
    cCV.addSample(xr[0:1,0:2])
    print("Vector is " + str(cCV))
    print("Decoding gives " + str(cCV.decode()))

    cCVsum = cp.ChannelVector(CCB)
    cCVsum[:] = cCV[0:1,:]+cCV[1:2,:]
    print("Average the two channels gives " + str(cCVsum))
    print("Decoding gives " + str(cCVsum.decode(2)))

    print("Encode and decode " + str(y[0:1,0:1]) + " as histogram")
    bCV = cp.ChannelVector(BCB)
    bCV.addSample(y[0:1,0:1])
    print("Vector is " + str(bCV))
    print("Decoding gives " + str(bCV.decode()))

    print("Add " + str(yr[0:1,1:2]) + " as histogram and decode")
    bCV.addSample(yr[0:1,1:2])
    print("Vector is " + str(bCV))
    print("Decoding gives " + str(bCV.decode(2)))

    bCV.fill(0.)

    print("Encode and decode " + str(yr[0:1,0:2]) + " as histogram")
    bCV.addSample(yr[0:1,0:2])
    print("Vector is " + str(bCV))
    print("Decoding gives " + str(bCV.decode()))

    bCVsum = cp.ChannelVector(BCB)
    bCVsum[:] = bCV[0:1,:]+bCV[1:2,:]
    print("Average the two channels gives " + str(bCVsum))
    print("Decoding gives " + str(bCVsum.decode(2)))

    print("Encode and decode " + str(xy[0:1,0:1]) + " as combined channels")
    cbCV = cp.ChannelVector(CBCB)
    cbCV.addSample(xy[0:1,0:1])
    print("Vector is " + str(cbCV))
    print("Decoding gives " + str(cbCV.decode()) + " compared to " + str(xyr[0:1,0:1]))

    print("Add " + str(xyr[0:1,1:2]) + " as combined channel and decode")
    cbCV.addSample(xyr[0:1,1:2])
    print("Vector is " + str(cbCV))
    print("Decoding gives " + str(cbCV.decode(2)) + " compared to " + str(xyr[0:1,0:2]))

    cbCV.fill(0.)
    
    print("Encode and decode " + str(xyr[0:1,0:2]) + " as combined channel")
    cbCV.addSample(xyr[0:1,0:2])
    print("Vector is " + str(cbCV))
    print("Decoding gives " + str(cbCV.decode()) + " compared to " + str(xyr[0:1,0:2]))

    cbCVsum = cp.ChannelVector(CBCB)
    cbCVsum[:] = cbCV[0:1,:]+cbCV[1:2,:]
    print("Average the two channels gives " + str(cbCVsum))
    print("Decoding gives " + str(cbCVsum.decode(2)) + " compared to " + str(xyr[0:1,0:2]))

    print("Phistogram with " + str(nCCB) + " coefficients from " + str(rCMin) + " to " + str(rCMax))

    CCB = cp.BilinearBasis()
    CCB.setParameters(nCCB, rCMin, rCMax)

    print("Phistogram with " + str(nBCB) + " coefficients from " + str(rBMin) + " to " + str(rBMax)) 
    BCB = cp.BilinearBasis()
    BCB.setParameters(nBCB, rBMin, rBMax)

    print("Cross product of channels")
    CBCB = cp.CombinedChannelBasis()
    chBV = [CCB, BCB]
    CBCB.setParameters(chBV)

    print("Channel vectors (non-sparse)")

    print("Encode and decode " + str(x[0:1,0:1]) + " as phistogram")
    cCV = cp.ChannelVector(CCB)
    cCV.addSample(x[0:1,0:1])
    print("Vector is " + str(cCV))
    print("Decoding gives " + str(cCV.decode()))

    print("Add " + str(xr[0:1,1:2]) + " as phistogram and decode")
    cCV.addSample(xr[0:1,1:2])
    print("Vector is " + str(cCV))
    print("Decoding gives " + str(cCV.decode(2)))

    cCV.fill(0.)

    print("Encode and decode " + str(xr[0:1,0:2]) + " as phistogram")
    cCV.addSample(xr[0:1,0:2])
    print("Vector is " + str(cCV))
    print("Decoding gives " + str(cCV.decode()))

    cCVsum = cp.ChannelVector(CCB)
    cCVsum[:] = cCV[0:1,:]+cCV[1:2,:]
    print("Average the two channels gives " + str(cCVsum))
    print("Decoding gives " + str(cCVsum.decode(2)))

    print("Encode and decode " + str(y[0:1,0:1]) + " as phistogram")
    bCV = cp.ChannelVector(BCB)
    bCV.addSample(y[0:1,0:1])
    print("Vector is " + str(bCV))
    print("Decoding gives " + str(bCV.decode()))

    print("Add " + str(yr[0:1,1:2]) + " as phistogram and decode")
    bCV.addSample(yr[0:1,1:2])
    print("Vector is " + str(bCV))
    print("Decoding gives " + str(bCV.decode(2)))

    bCV.fill(0.)

    print("Encode and decode " + str(yr[0:1,0:2]) + " as phistogram")
    bCV.addSample(yr[0:1,0:2])
    print("Vector is " + str(bCV))
    print("Decoding gives " + str(bCV.decode()))

    bCVsum = cp.ChannelVector(BCB)
    bCVsum[:] = bCV[0:1,:]+bCV[1:2,:]
    print("Average the two channels gives " + str(bCVsum))
    print("Decoding gives " + str(bCVsum.decode(2)))

    print("Encode and decode " + str(xy[0:1,0:1]) + " as combined channels")
    cbCV = cp.ChannelVector(CBCB)
    cbCV.addSample(xy[0:1,0:1])
    print("Vector is " + str(cbCV))
    print("Decoding gives " + str(cbCV.decode()) + " compared to " + str(xyr[0:1,0:1]))

    print("Add " + str(xyr[0:1,1:2]) + " as combined channel and decode")
    cbCV.addSample(xyr[0:1,1:2])
    print("Vector is " + str(cbCV))
    print("Decoding gives " + str(cbCV.decode(2)) + " compared to " + str(xyr[0:1,0:2]))

    cbCV.fill(0.)
    
    print("Encode and decode " + str(xyr[0:1,0:2]) + " as combined channel")
    cbCV.addSample(xyr[0:1,0:2])
    print("Vector is " + str(cbCV))
    print("Decoding gives " + str(cbCV.decode()) + " compared to " + str(xyr[0:1,0:2]))

    cbCVsum = cp.ChannelVector(CBCB)
    cbCVsum[:] = cbCV[0:1,:]+cbCV[1:2,:]
    print("Average the two channels gives " + str(cbCVsum))
    print("Decoding gives " + str(cbCVsum.decode(2)) + " compared to " + str(xyr[0:1,0:2]))
