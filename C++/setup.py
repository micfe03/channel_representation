from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
    
module1 = Extension("ChannelBasis",
                    ["PyChannelBasis.pyx"],
                    include_dirs=["/usr/local/include/",
                                  "/Library/Python/2.7/site-packages/numpy-1.8.0.dev_436a28f_20120710-py2.7-macosx-10.7-x86_64.egg/numpy/core/include/"],
                    libraries=["Channelbasis"],
                    library_dirs=["Release"],
                    language="c++")
                  
setup(cmdclass = {"build_ext": build_ext}, ext_modules = [module1])


