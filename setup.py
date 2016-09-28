from distutils.core import setup, Extension
import numpy.distutils.misc_util

ext = Extension("rotsampling", sources=["rotsamplingmodule.c"], include_dirs=[numpy.get_include()], extra_compile_args=["-std=c99"])
#setup(ext_modules=[ext], include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs())
setup(ext_modules=[ext], include_dirs=[numpy.get_include()])
