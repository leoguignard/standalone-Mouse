# I2AE

This repository contains the software to register an annotated image data onto our average embryo. This software is described in our Mouse-Atlas article: A dynamic atlas of mouse development: multi-scale analysis of gastrulation and early organogenesis

## Description of the repository
Folders:
  - Brachyury-embryo-data: folder place-holder meant to contain the test image data
  - Landmarks: folder place-holder meant to contain the manual landmarks
  - avg-embryo: folder place-holder meant to contain the Statistical Vector Flow of the average embryo
  - config-csv: folder containing an example of a configuration csv file

Files:
  - Readme.md: this file.
  - image-to-average-embryo.py: the script file.

## Basic usage
You can run the script the following way:
```shell
python image-to-average-embryo.py
```
It will then ask the path to the configuration file for the algorithm.

## Dependencies
The following basic dependencies are required:
  - numpy
  - scipy
  
To install these two dependencies, one can run:
```Shell
python setup.py install
```
  
Then, the following two "special" dependencies are required:
  - TGMMlibraries (from the other software available with that article)
  - IO (from the other software available with that article)
