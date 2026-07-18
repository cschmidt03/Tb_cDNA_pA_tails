from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import os

CONDA_PREFIX = os.environ["CONDA_PREFIX"]

ext_modules = [
    Extension(
        name="cy_stringscan",
        sources=["./cython/cy_stringscan.pyx"],
        language="c++",
        extra_compile_args = [
            "-Wno-unused-result", "-Wno-sign-compare"
        ]
   ),
   Extension(
       name="cy_bamindexer",
       sources=["./cython/cy_bamindexer.pyx"],
        include_dirs=[os.path.join(CONDA_PREFIX, "include")],
        library_dirs=[os.path.join(CONDA_PREFIX, "lib")],
        libraries=["hts"],
        language="c++",
        extra_compile_args = [
            "-Wno-unused-result", "-Wno-sign-compare"
        ]
   )
]

setup(
    ext_modules = cythonize(ext_modules)
)
