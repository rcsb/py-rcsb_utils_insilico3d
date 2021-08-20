# RCSB Python In Silico Model Access Utilities

## Introduction

This module contains utility methods for accessing in silico 3D models and metadata from external data resources, including AlphaFold, ModBase, SWISS-MODEL, and Model Archive.

### Installation

Download the library source software from the project repository:

```bash

git clone --recurse-submodules https://github.com/rcsb/py-rcsb_utils_insilico3d.git

```

**Important:** Setup will require an up-to-date version of [cmake](https://cmake.org/install/) to be installed on the machine and the executable to be in the system's PATH.

Optionally, run test suite (Python versions 3.9) using
[setuptools](https://setuptools.readthedocs.io/en/latest/) or
[tox](http://tox.readthedocs.io/en/latest/example/platform.html):

```bash
python setup.py test

or simply run

tox
```

Installation is via the program [pip](https://pypi.python.org/pypi/pip).

```bash
pip install rcsb.utils.insilico3d

or for the local repository:

pip install .
```

## References
* Jumper, J et al. Highly accurate protein structure prediction with AlphaFold. Nature (2021)
* Bienert, S., Waterhouse, A., de Beer, T.A.P., Tauriello, G., Studer, G., Bordoli, L., Schwede, T. The SWISS-MODEL Repository - new features and functionality. Nucleic Acids Res. 45, D313-D319 (2017). 
* Waterhouse, A., Bertoni, M., Bienert, S., Studer, G., Tauriello, G., Gumienny, R., Heer, F.T., de Beer, T.A.P., Rempfer, C., Bordoli, L., Lepore, R., Schwede, T. SWISS-MODEL: homology modelling of protein structures and complexes. Nucleic Acids Res. 46(W1), W296-W303 (2018). 

