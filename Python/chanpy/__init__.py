import numpy as np
import cmath

class ChannelBasis:
    
    def __init__(self):
        self.setParameters(self)
    
    def setParameters(self, nChannels=11, minValue=0., maxValue=0., modularChannels=False, boundedChannels=True):
        self.__modular__ = modularChannels
        self.__bounded__ = boundedChannels
        self.__minV__ = minValue
        self.__maxV__ = maxValue
        self.__width__ = np.pi/3.0
        self.__nChans__ = nChannels
#        self.__norm__ = 0.
        self.__updateInputScaling__()
    
    def displayParameters(self):
        self.displayBasisType()
        print(' nchans       : ' + str(self.__nChans__))
        print(' [minv,maxv]  : [' + str(self.__minV__) + ',' + str(self.__maxV__) + ']')
        print(' modular      : ' + str(self.__modular__))
        print(' bounded      : ' + str(self.__bounded__))
    
    def displayBasisType(self):
        print(' Channel type : empty')
    
    def getNrChannels(self):
        return self.__nChans__
    
    def getNrChannelsVec(self):
        nChansVec = (self.getNrChannels(),)
        return nChansVec
    
    def getMaxV(self):
        return self.__maxV__
    
    def getMinV(self):
        return self.__minV__
    
    def getNorm(self):
        return self.__norm__
    
    def encode(self, values):
        pass
    
    def decode(self, chCoeff, nrModes=1):
        pass
    
    def __scaleInput__(self, value):
        if self.__bounded__:
            value = np.maximum(value,self.__minV__)
            value = np.minimum(value,self.__maxV__)
        value = self.__scaling__*(value-self.__minV__)+self.__offset__
        return value
    
def sortModes(resVec, confVec, res, certainty):
    nModes = res.size
    if nModes == 1:
        for i in range(resVec.size):
            if confVec[i] > certainty[0]:
                res[0] = resVec[i]
                certainty[0] = confVec[i]
        return
    for i in range(resVec.size):
        if confVec[i] > certainty[nModes-1]:
            certainty[nModes-1] = confVec[i]
            res[nModes-1] = resVec[i]
            for j in range(nModes-1, 0, -1):
                if certainty[j] < certainty[j-1]:
                    break
                tmpC = certainty[j-1]
                certainty[j-1] = certainty[j]
                certainty[j] = tmpC
                tmpC = res[j-1]
                res[j-1] = res[j]
                res[j] = tmpC

class HistogramBasis(ChannelBasis):
    
    def __init__(self):
        self.__norm__ = 1.0

    def __updateInputScaling__(self):
        self.__scaling__ = self.__nChans__/(self.__maxV__ - self.__minV__)
        self.__offset__ = -0.5
        
    def encode(self, cVal):
        val = cVal
        val = self.__scaleInput__(val)
        chVec = np.zeros((1,self.__nChans__))
        i_centerCoeff = np.round(val).astype(int)
        i_centerCoeff = np.minimum(i_centerCoeff,self.__nChans__ - 1)
        i_centerCoeff = np.maximum(i_centerCoeff,0)
        chVec[0,i_centerCoeff]   = 1
        return chVec
    
    def displayBasisType(self):
        print(' Channel type : rectangular')
    
    def decode(self, chCoeff, nModes=1):
        resVec    = np.zeros(self.__nChans__)
        result    = np.zeros(nModes)
        certainty = -np.ones(nModes)
        res       = np.zeros((1,2*nModes))
        for i in range(self.__nChans__ ):
            resVec[i] = (i - self.__offset__)/self.__scaling__ + self.__minV__
        sortModes(resVec,chCoeff,result,certainty)
        for i in range(nModes):
            res[0,2*i]   = result[i]
            res[0,2*i+1] = certainty[i]/self.__norm__
        return res

class BilinearBasis(ChannelBasis):
    
    def __init__(self):
        self.__norm__ = 1.0
        
    def __updateInputScaling__(self):
        if self.__modular__:
            self.__scaling__ = (self.__nChans__)/(self.__maxV__ - self.__minV__)
        else:
            self.__scaling__ = (self.__nChans__ - 1.0)/(self.__maxV__ - self.__minV__)
        self.__offset__ = 0.0
        
    def encode(self, cVal):
        val = cVal
        val = self.__scaleInput__(val)
        i_centerCoeff = np.floor(val).astype(int)
        i_centerCoeff = np.maximum(i_centerCoeff,0)
        if self.__modular__:
            i_centerCoeff = np.minimum(i_centerCoeff,self.__nChans__ - 1)
            chVec = np.zeros((1,self.__nChans__+1))
        else:
            i_centerCoeff = np.minimum(i_centerCoeff,self.__nChans__ - 2)
            chVec = np.zeros((1,self.__nChans__))
        d = val - i_centerCoeff
        chVec[0,i_centerCoeff]   = 1.0 - d
        chVec[0,i_centerCoeff+1]   = d
        if self.__modular__:
            chVec[0,0] += chVec[0,-1]
            chVec = chVec[:,:-1]
        return chVec
    
    def displayBasisType(self):
        print(' Channel type : linear B-spline')
    
    def decode(self, chCoeff, nModes=1):
        resVec    = np.zeros(self.__nChans__ - 1)
        confVec   = np.zeros(self.__nChans__ - 1)
        result    = np.zeros(nModes)
        certainty = -np.ones(nModes)
        res       = np.zeros((1,2*nModes))
        for i in range(self.__nChans__ - 1):
            confVec[i] = chCoeff[i] + chCoeff[i+1]
            resVec[i] = (chCoeff[i+1] - chCoeff[i])/(2*confVec[i]+0.000001) + 0.5 + i
            if (resVec[i]<i) or (resVec[i]>i+1.):
                resVec[i]  = 0
                confVec[i] = 0
            else:
                resVec[i] = (resVec[i] - self.__offset__)/self.__scaling__ + self.__minV__
        sortModes(resVec,confVec,result,certainty)
        for i in range(nModes):
            res[0,2*i]   = result[i]
            res[0,2*i+1] = certainty[i]/self.__norm__
        return res
        
class Cos2ChannelBasis(ChannelBasis):
    
    def __init__(self):
        self.__norm__ = 1.5

    def __updateInputScaling__(self):
        self.__scaling__ = (self.__nChans__ - 2.)/(self.__maxV__ - self.__minV__)
        self.__offset__ = 0.5
    
    def encode(self, cVal):
        val = cVal
        val = self.__scaleInput__(val)
        chVec = np.zeros((1,self.__nChans__))
        i_centerCoeff = np.round(val).astype(int)
        i_centerCoeff = np.minimum(i_centerCoeff,self.__nChans__ - 2)
        i_centerCoeff = np.maximum(i_centerCoeff,1)
        d = val - i_centerCoeff
        chVec[0,i_centerCoeff-1] = np.cos((1+d) * self.__width__)**2
        chVec[0,i_centerCoeff]   = np.cos(d * self.__width__)**2
        chVec[0,i_centerCoeff+1] = np.cos((1-d) * self.__width__)**2
        return chVec
    
    def displayBasisType(self):
        print(' Channel type : cos2')
    
    def decode(self, chCoeff, nModes=1):
        resVec    = np.zeros(self.__nChans__ - 2)
        confVec   = np.zeros(self.__nChans__ - 2)
        result    = np.zeros(nModes)
        certainty = -np.ones(nModes)
        res       = np.zeros((1,2*nModes))
        argument  = 0.
        for i in range(self.__nChans__ - 2):
            argument = cmath.phase(chCoeff[i] + chCoeff[i+1]*cmath.exp(2j*self.__width__) + chCoeff[i+2]*cmath.exp(4j*self.__width__))
            if argument<0:
                argument += 2*np.pi
            resVec[i]  = i + argument/(2*self.__width__)
            confVec[i] = chCoeff[i] + chCoeff[i+1] + chCoeff[i+2]
        for i in range(self.__nChans__ - 2):
            if (resVec[i]<i+0.5) or (resVec[i]>i+1.5):
                resVec[i]  = 0
                confVec[i] = 0
            else:
                resVec[i] = (resVec[i] - self.__offset__)/self.__scaling__ + self.__minV__
        sortModes(resVec,confVec,result,certainty)
        for i in range(nModes):
            res[0,2*i]   = result[i]
            res[0,2*i+1] = certainty[i]/self.__norm__
        return res
        
class CombinedChannelBasis(ChannelBasis):
    
    def __init__(self):
        pass
    
    def setParameters(self, chBasisVector):
        self.__chBasisVector__ = chBasisVector
        self.__nChans__        = 1
        self.__norm__          = 1.
        self.__nChansVec__     = np.zeros(len(chBasisVector))
        for k in range(len(chBasisVector)):
            nPartChannels = chBasisVector[k].getNrChannels()
            self.__nChans__ *= nPartChannels
            self.__norm__ *= chBasisVector[k].getNorm()
            self.__nChansVec__[k] = nPartChannels
    
    def getNrChannelsVec(self):
        nrChansVec = self.__nChansVec__
        return nrChansVec
    
    def encode(self, vals):
        tmpxChannel = ChannelVector(self.__chBasisVector__[0])
        tmpxChannel.addSample(np.atleast_2d(vals[0]))
        #tmpOP = ChannelVector(self.__chBasisVector__[0])
        for k in range(1,len(vals)):
            tmpyChannel = ChannelVector(self.__chBasisVector__[k])
            tmpyChannel.addSample(np.atleast_2d(vals[k]))
            #sizes = (1,tmpxChannel.shape[1]*tmpyChannel.shape[1])
            #tmpOP.resize(sizes,refcheck=0)
            #tmpOP[:] = np.kron(tmpyChannel,tmpxChannel)
            #tmpxChannel.resize(sizes,refcheck=0)
            #tmpxChannel[:] = tmpOP.copy()
            tmpxChannel=np.kron(tmpyChannel,tmpxChannel)
        tmpxChannel.setChannelBasis(self)
        tmpxChannel.normalize()
        return tmpxChannel
    
    def displayBasisType(self):
        print('Combined Channel Basis')
    
    def decode(self, chCoeff, nrModes=1):
        res = np.zeros((nrModes,len(self.__chBasisVector__)+1))
        for mode in range(nrModes):
            nrXChannels = chCoeff.shape[0]
            xChannel = chCoeff.copy()
            for k in range(len(self.__chBasisVector__)-1,-1,-1):
                nrYChannels = self.__chBasisVector__[k].getNrChannels()
                nrXChannels //= nrYChannels
                # xChannel = xChannel.reshape((nrYChannels,nrXChannels))
                # tmpCoeff = np.sum(xChannel,1)
                tmpCoeff = np.zeros((nrYChannels))
                for l in range(nrYChannels):
                    for xCi in range(l*nrXChannels, (l+1)*nrXChannels):
                        tmpCoeff[l] += xChannel[xCi]
                # 
                if k==len(self.__chBasisVector__)-1:
                    pRes = self.__chBasisVector__[k].decode(tmpCoeff,mode+1)
                    res[mode,k] = pRes[0,2*mode]
                    tmpCoeff = self.__chBasisVector__[k].encode(pRes[0,2*mode])[0,:]
                    tmpCoeff *= pRes[0,2*mode+1]
                else:
                    pRes = self.__chBasisVector__[k].decode(tmpCoeff,1)
                    res[mode,k] = pRes[0,0]
                    tmpCoeff = self.__chBasisVector__[k].encode(pRes[0,0])[0,:]
                    tmpCoeff *= pRes[0,1]
                elem_counter = 0
                for tCi in range(len(tmpCoeff)):
                    xChannel[nrXChannels*elem_counter:nrXChannels*(elem_counter+1)] *= tmpCoeff[tCi]
                    elem_counter += 1
                xChannel = np.sqrt(xChannel)
                tmpCoeff = np.zeros((nrXChannels))
                for r in range(nrYChannels):
                    tmpCoeff += xChannel[r*nrXChannels:(r+1)*nrXChannels]
                xChannel = (tmpCoeff**2)/self.__chBasisVector__[0].getNorm()
                #xChannel = tmpCoeff
            res[mode,len(self.__chBasisVector__)] = xChannel[0]**0.25/self.__chBasisVector__[0].getNorm()**1.25
        res = res.reshape((1,-1))
        return res
    
class ChannelVector(np.ndarray):

    def __new__(cls, chBasis=None, coeffs=None):
        if chBasis is None:
            obj = np.array([[]]).view(cls)
            obj.__chBasis__ = chBasis
            obj.__support__ = ()
        else:
            if coeffs is None:
                obj = np.zeros((1,chBasis.getNrChannels())).view(cls).copy() #penalty?
                obj.__support__ = (1,1)
            else:
                obj = coeffs.view(cls)
                obj.__support__ = (1,coeffs.shape[0])
            obj.__chBasis__ = chBasis
        return obj
        
    def addSample(self, vals):
        sizes = self.__chBasis__.getNrChannelsVec()
        vals = np.atleast_3d(vals)
        nrR, nrC, nrK = vals.shape
        sizeDiff = nrK-len(sizes)
        createCCFM = (sizeDiff<0)
        dedicatedCertainty = (abs(sizeDiff)==1)
        cVals = np.empty(len(sizes))
        if createCCFM:
            nrEls = 1
        else:
            nrEls = nrR*nrC
        if (self.shape[1] != self.__chBasis__.getNrChannels()) or (self.shape[0] != nrEls):
            #self = ChannelVector(self.__chBasis__,np.zeros((nrEls,self.__chBasis__.getNrChannels())))
            self.resize((nrEls,self.__chBasis__.getNrChannels()),refcheck=0)
            self.fill(0.)
            if createCCFM:
                self.__support__ = (1,1)
            else:
                self.__support__ = (nrR,nrC)
        for rdx in range(nrR):
            for cdx in range(nrC):
                if dedicatedCertainty:
                    cVals[:nrK-1] = vals[rdx,cdx,:-1]
                else:
                    cVals[:nrK] = vals[rdx,cdx,:]
                if createCCFM:
                    cVals[-2] = rdx
                    cVals[-1] = cdx
                tmp_chCoeff = self.__chBasis__.encode(cVals)
                if createCCFM:
                    corrW = (nrR*nrC)**(2./sizes.size) # maybe bug
                    if dedicatedCertainty:
                        self += tmp_chCoeff*vals[rdx,cdx,-1]/corrW
                    else:
                        self += tmp_chCoeff/corrW
                else:
                    if dedicatedCertainty:
                        self[rdx*nrC+cdx,:] += tmp_chCoeff[0,:]*vals[rdx,cdx,-1]
                    else:
                        self[rdx*nrC+cdx,:] += tmp_chCoeff[0,:]
                        
    def decode(self, nrModes=1):
        sizes = self.__chBasis__.getNrChannelsVec()
        nrChannels = nrModes*(len(sizes)+1)
        nrEls = self.__support__[0]*self.__support__[1]
        localRes = np.zeros((nrEls,nrChannels))
        for cdx in range(self.shape[0]):
            localRes[cdx,:] = self.__chBasis__.decode(self[cdx,:],nrModes)
        res = localRes.reshape(self.__support__ + (-1,))
        return res
    
    def normalize(self):
        nrm = np.sum(self)
        if nrm > 0:
            self *= self.__chBasis__.getNorm()/nrm
        else:
            self.fill(self.__chBasis__.getNorm()/self.shape[1])
    
    def setChannelBasis(self, chBasis):
        self.__chBasis__ = chBasis
        #self = self.base.resize((1,chBasis.getNrChannels())).view(self.__class__)
        #self.fill(0.)
        self.__support__ = (1,1)
    
    def channelImage(self):
        res = self.reshape(self.__support__ + (-1,))
        return res
        
    def histogramMatrix(self):
        pass
    # std::vector<int> sizes(1);
    # m_chBasis->getNrChannelsVec(sizes);
    # sizes.push_back(rows);
    # cv::Mat resTmp((int)sizes.size(),&sizes[0],CV_32F,data);
    # res = resTmp;
