# Time-registration

This repository contains the Time registration script proposed in our Mouse-Atlas article.

## Description of the repository
Folders:
  - csv-parameter-files: Example of parameterization csv file for standalone-registration.py.
  
Python files:
  - standalone-registration.py: python script to register in time a time-series of 3D intensity images

## Basic usage
The python scripts proposed here can be run from a terminal the following way:

`python standalone-registration.py`

It is then prompted to provide a parameter csv file (examples provided in the folder csv-parameter-files). The path to the parameter file should be then typed in the terminal.

**If the package have been installed using setup.py, the script can be call just by typing `TARDIS.py` from anywhere in the terminal**


## Dependencies
Some dependecies are requiered:
  - general python dependecies:
    - numpy, scipy
  - standalone-registration.py:
    - blockmatching installed such that it can be called in a terminal as > blockmatching -ref [...] (https://gitlab.inria.fr/greg/Klab-BlockMatching)
    - IO library has to be installed 
    - psutil python library
    - scikitimage
    
## Quick install
To quickly install the script so it can be call from the terminal and install too the common dependecies one can run
```shell
python setup.py install [--user]
```
Still will be remaining to install blockmathcing and IO packages.
