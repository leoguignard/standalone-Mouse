# IO

This library allows to read and write 3D images from different sources. It also has the SpatialImage class that extend numpy arrays so they also contain the voxel size information.

This library is an extension of the one proposed by [VirtualPlants](https://team.inria.fr/virtualplants/) and that can be found there: https://github.com/openalea/openalea-components (and specifically [there](https://github.com/openalea/openalea-components/tree/master/image/src/openalea/image/serial))

The library can read and write 3D images from the following format: tiff, klb, hdf5 (only read, not write for the moment).
It can also read a stack of 2D images from a folder (useful when the stack is really big and hdf5/klb format are not available).

## Description of the repository
  - IO: folder containing the package
  - setup.py: Installation script
  - README.md: This file

## Basic usage
Once installed the library can be called the following way (as an example):
```python
from IO import imread, imsave, SpatialImage
```

## Dependencies
Some dependecies are requiered:
  - general python dependecies:
    - numpy, scipy
  - Image reading dependencies:
    - h5py (https://pypi.python.org/pypi/h5py)
    - pyklb (https://github.com/bhoeckendorf/pyklb)
    - pylibtiff (https://github.com/pearu/pylibtiff)

Note that all the dependecies but pyklb will be installed automatically with the setup.py (therefore you need to manually install pyklb correctly in order to read/write klb images).
    
## Quick install
To quickly install the script so it can be call from the terminal and install too the common dependecies one can run
```shell
python setup.py install [--user]
```
Still will be remaining to install pyklb package.
