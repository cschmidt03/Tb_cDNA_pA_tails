from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import os

include_dirs = []
library_dirs = []

if "CONDA_PREFIX" in os.environ:
    CONDA_PREFIX = os.environ["CONDA_PREFIX"]
    include_dirs.append(os.path.join(CONDA_PREFIX, "include"))
    library_dirs.append(os.path.join(CONDA_PREFIX, "lib"))

ext_modules = [
    Extension(
        name="cy_stringscan",
        sources=["./cy_stringscan.pyx"],
        language="c++",
        extra_compile_args = [
            "-Wno-unused-result", "-Wno-sign-compare"
        ]
   ),
   Extension(
       name="cy_bamindexer",
       sources=["./cy_bamindexer.pyx"],
        include_dirs=include_dirs,
        library_dirs=library_dirs,
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
