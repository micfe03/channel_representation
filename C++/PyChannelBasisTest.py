import ChannelBasis as cb
import numpy as np
import numpy.random as rnd
import cv2

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

print "PyChannelBasisTest"

print "Cos2-channels with", nCCB, "coefficients from", rCMin, "to", rCMax

CCB = cb.PyCos2ChannelBasis()
CCB.setParameters(nCCB, rCMin, rCMax)

print "Bspline-channels with", nBCB, "coefficients from", rBMin, "to", rBMax 
BCB = cb.PyBsplineChannelBasis()
BCB.setParameters(nBCB, rBMin, rBMax)

print "Cross product of channels"
CBCB = cb.PyCombinedChannelBasis()
chBV = [CCB, BCB]
CBCB.setParameters(chBV)

print "Channel vectors (non-sparse)"

print "Encode and decode", x[0:1,0:1], "as cos2 channel"
cCV = cb.PyChannelVector(CCB)
cCV.addSample(x[0:1,0:1])
print "Vector is", cCV.asarray()
print "Decoding gives", cCV.decode()

print "Add", xr[0:1,1:2], "as cos2 channel and decode"
cCV.addSample(xr[0:1,1:2])
print "Vector is", cCV.asarray()
print "Decoding gives", cCV.decode(2)

cCV.asarray()[:] = 0

print "Encode and decode", xr[0:1,0:2], "as cos2 channel"
cCV.addSample(xr[0:1,0:2])
print "Vector is", cCV.asarray()
print "Decoding gives", cCV.decode()

cCVsum = cb.PyChannelVector(CCB)
cCVsum.asarray()[:] = cCV.asarray()[0:1,:]+cCV.asarray()[1:2,:]
print "Average the two channels gives", cCVsum.asarray()
print "Decoding gives",  cCVsum.decode(2)

print "Encode and decode", y[0:1,0:1], "as bspline channel"
bCV = cb.PyChannelVector(BCB)
bCV.addSample(y[0:1,0:1])
print "Vector is", bCV.asarray()
print "Decoding gives", bCV.decode()

print "Add", yr[0:1,1:2], "as bspline channel and decode"
bCV.addSample(yr[0:1,1:2])
print "Vector is", bCV.asarray()
print "Decoding gives",  bCV.decode(2)

bCV.asarray()[:] = 0

print "Encode and decode", yr[0:1,0:2], "as bspline channel"
bCV.addSample(yr[0:1,0:2])
print "Vector is", bCV.asarray()
print "Decoding gives", bCV.decode()

bCVsum = cb.PyChannelVector(BCB)
bCVsum.asarray()[:] = bCV.asarray()[0:1,:]+bCV.asarray()[1:2,:]
print "Average the two channels gives", bCVsum.asarray()
print "Decoding gives",  bCVsum.decode(2)

print "Encode and decode", xy[0:1,0:1], "as combined channels"
cbCV = cb.PyChannelVector(CBCB)
cbCV.addSample(xy[0:1,0:1])
print "Vector is", cbCV.asarray()
print "Decoding gives", cbCV.decode()

print "Add", xyr[0:1,1:2], "as combined channel and decode"
cbCV.addSample(xyr[0:1,1:2])
print "Vector is", cbCV.asarray()
print "Decoding gives", cbCV.decode(2)

cbCV.asarray()[:] = 0

print "Encode and decode", xyr[0:1,0:2], "as combined channel"
cbCV.addSample(xyr[0:1,0:2])
print "Vector is", cbCV.asarray()
print "Decoding gives", cbCV.decode()

cbCVsum = cb.PyChannelVector(CBCB)
cbCVsum.asarray()[:] = cbCV.asarray()[0:1,:]+cbCV.asarray()[1:2,:]
print "Average the two channels gives", cbCVsum.asarray()
print "Decoding gives", cbCVsum.decode(2)

image = cv2.imread("house_orig.png")
image = cv2.cvtColor(image,cv2.cv.CV_RGB2GRAY)
image = image[0:248,:]
cv2.imshow("Original Image", image)
CCBB = cb.PyCos2ChannelBasis()
CCBB.setParameters(10, 0, 255)
cCVB = cb.PyChannelVector(CCBB)
cCVB.addSample(image)
cCVB.channelImage()[:] = cv2.GaussianBlur(cCVB.channelImage(), ksize = (7, 7), sigmaX = 1.5, sigmaY = 1.5)
fImage2 = cCVB.decode(2)
cv2.imshow("Mode0 Image", cv2.convertScaleAbs(fImage2[:,:,0].squeeze()))
cv2.imshow("Mode1 Image", cv2.convertScaleAbs(fImage2[:,:,2].squeeze()))
cv2.waitKey(0)

CCBC = cb.PyCos2ChannelBasis()
CCBC.setParameters(14, 0, image.shape[1]-1)
CCBR = cb.PyCos2ChannelBasis()
CCBR.setParameters(12, 0, image.shape[0]-1)
CCBB.setParameters(10, 0, 1)
image = 1.0/255.0*image
print image.shape, image[100,250]
CCFM = cb.PyCombinedChannelBasis()
chBCCFM = [CCBB, CCBR, CCBC]
CCFM.setParameters(chBCCFM)
cCCFM = cb.PyChannelVector(CCFM)
cCCFM.addSample(image)
print "Decoding gives", cCCFM.decode(1)
