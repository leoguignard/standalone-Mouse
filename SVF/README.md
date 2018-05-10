# SVF

This repository contains the Statistical Vector Flow and the tissue propoagation software proposed in our Mouse-Atlas article.

## Description of the repository
Folder:
  - csv-parameter-files: Example of parameterization csv files for each algorithms.
  
Python files:
  - SVF-prop.py: python script to build Statistical Vector Flow from a TGMM dataset.
  - tissue-bw-prop.py: python script to propagate tissue information from a manually annotated 3D image.

## Basic usage
Each of the python scripts proposed here can be run from a terminal the following way:

`python SVF-prop.py`

`python tissue-bw-prop.py`

It is then prompted to provide a parameter csv file (examples provided in the folder csv-parameter-files). The path to the parameter file should be then typed in the terminal.

**If the package have been installed using setup.py, the script can be call just by typing `SVF-prop.py` or `tissue-bw-prop.py` from anywhere in the terminal**


## Dependencies
Some dependecies are requiered:
  - general python dependecies:
    - numpy, scipy, pandas
  - SVF-prop.py:
     - TGMMlibraries has to be installed (https://github.com/leoguignard/TGMMlibraries)
  - tissue-bw-prop.py:
    - TGMMlibraries has to be installed (https://github.com/leoguignard/TGMMlibraries)
    - IO library has to be installed (https://github.com/leoguignard/IO)

## Quick install
To quickly install the script so it can be call from the terminal and install too the common dependecies one can run
```shell
python setup.py install [--user]
```
Still will be remaining to install IO and TGMMlibraries packages.
