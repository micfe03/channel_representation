# distutils: language = c++

from libcpp.vector cimport vector
cimport numpy as np
from cpython cimport PyObject
from cv2 import CV_32F
import numpy as np

np.import_array()

cdef extern from "opencv2/opencv.hpp" namespace "cv":
    cdef cppclass Mat:
        Mat()
        Mat(int,int,int,void*)
        int rows
        int cols
        float* data
        int channels()

cdef extern from "ChannelBasis.hpp" namespace "cvl":
    cdef cppclass ChannelBasis:
        void setParameters( int,  float,  float,  int,  int)
        ChannelBasis()
        
    cdef cppclass Cos2ChannelBasis:
        Cos2ChannelBasis()
        
    cdef cppclass BsplineChannelBasis:
        BsplineChannelBasis()
        
    cdef cppclass CombinedChannelBasis:
        CombinedChannelBasis()
        void setParameters(vector[ChannelBasis*])
        
    cdef cppclass ChannelVector:
        ChannelVector(ChannelBasis*)
#        ChannelVector(ChannelBasis*, Mat)
        int rows
        int cols
        float* data
#        ChannelBasis* m_chBasis
        void addSample(Mat)
        void decode(Mat, int)
        void channelImage(Mat)
     
cdef class PyChannelBasis:
    cdef ChannelBasis *thisptr
    def __cinit__(self):
        pass
#        self.thisptr = new ChannelBasis()
    def __dealloc__(self):
        pass
#        del self.thisptr
    def setParameters(self,  int nChannels,  float minValue,  float maxValue,  int modularChannels=0,  int boundedChannels=1):
        self.thisptr.setParameters(nChannels, minValue, maxValue, modularChannels, boundedChannels)

cdef class PyCos2ChannelBasis(PyChannelBasis):
    def __cinit__(self):
        self.thisptr = <ChannelBasis*> new Cos2ChannelBasis()
    def __dealloc__(self):
        del self.thisptr
    
cdef class PyBsplineChannelBasis(PyChannelBasis):
    def __cinit__(self):
        self.thisptr = <ChannelBasis*> new BsplineChannelBasis()
    def __dealloc__(self):
        del self.thisptr

cdef class PyCombinedChannelBasis(PyChannelBasis):
    def __cinit__(self):
        self.thisptr = <ChannelBasis*> new CombinedChannelBasis()
    def __dealloc__(self):
        del self.thisptr
    def setParameters(self, PychBasisVector):
        cdef vector[ChannelBasis*] chBasisVector
        cdef ChannelBasis* pdummy
        for x in PychBasisVector:
            pdummy = (<PyChannelBasis>x).thisptr
            chBasisVector.push_back(pdummy)
        (<CombinedChannelBasis*>self.thisptr).setParameters(chBasisVector)
    
cdef class PyChannelVector:
    cdef ChannelVector *thisptr
    def __cinit__(self, PyChannelBasis PychBasis):
        cdef ChannelBasis* chBasis = PychBasis.thisptr
        self.thisptr = new ChannelVector(chBasis)
    def __dealloc__(self):
        del self.thisptr
#    def fromArray(self, PyChannelBasis PychBasis, coeffs):
#        cdef ChannelBasis* chBasis = PychBasis.thisptr
#        del self.thisptr
#        cdef np.ndarray[np.float32_t, ndim = 2, mode = 'c'] toCV = np.ascontiguousarray(coeffs,dtype = np.float32)
#        cdef Mat cvCoeffs
#        cvCoeffs = Mat(np.size(coeffs,0),np.size(coeffs,1),CV_32F,toCV.data)        
#        self.thisptr = new ChannelVector(chBasis,cvCoeffs)
    def asarray(self):
        return np.PyArray_SimpleNewFromData(2,[self.thisptr.rows,self.thisptr.cols],np.NPY_FLOAT32,self.thisptr.data)
    def addSample(self, vals):
        if vals.ndim == 2:
            vals = np.array([vals.T]).T
        cdef Mat cvVals
        cdef np.ndarray[np.float32_t, ndim = 3, mode = 'c'] valsa = np.ascontiguousarray(vals,dtype = np.float32)
        cdef int CVtype
        CVtype = ((CV_32F&7) + ((np.size(vals,2)-1) << 3))
        cvVals = Mat(np.size(vals,0),np.size(vals,1),CVtype,valsa.data)
        self.thisptr.addSample(cvVals)
    def decode(self, nrModes=1):
        cdef Mat* cvRes
        cvRes = new Mat()
        self.thisptr.decode(cvRes[0], nrModes)
        return np.PyArray_SimpleNewFromData(3,[cvRes.rows,cvRes.cols,cvRes.channels()],np.NPY_FLOAT32,cvRes.data)
    def channelImage(self):
        cdef Mat cvRes
        self.thisptr.channelImage(cvRes)
        return np.PyArray_SimpleNewFromData(3,[cvRes.rows,cvRes.cols,cvRes.channels()],np.NPY_FLOAT32,cvRes.data)