## Compile the Cython extensions

The extensions can be compiled by executing `analysis_setup/cython/compile_cython.sh`.
This script runs the Python setup script and moves the files to the `scripts/` directory.

### Requirements 
* `setuptools` python module (e.g. install via `pip`)
* `cython` python module (available via `pip install Cython`)
* a C++ compiler (on Linux usually present as `g++`, more information [here](https://cython.readthedocs.io/en/latest/src/quickstart/install.html))
* the `HTSlib` headers & libraries. These may be already present by a `samtools` installation. If the compilation script is run inside a Conda environment, it will check the environment's `include/` and `lib/` directories first. HTSlib can be installed for a conda environment with `conda install bioconda::htslib`