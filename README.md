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
