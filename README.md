# Standalone version

This folder contains the scripts described in the article *In toto imaging and reconstruction of post-implantation mouse development at the single-cell level*. The algorithms present here are the main python algorithms:
  - Time-registration/standalone-registration.py that allows to stabilize in time 4D (3D + t) time series (probably the most complicated to install)
  - SVF/SVF-prop.py that build the Statistical Vector Flow from TGMM data
  - SVF/tissue-bw-prop.py that allows to propagate backward and forward masks from image data onto a corresponding SVF
  - I2AE/image-to-average-embryo.py that allows to register a masked image onto our average embryo.
  - svf2MaMuT/SVF2MaMuT.py that allows to rewrite an SVF binary together with its Database.csv (from the tissue-bw-prop.py algorithm) onto a MaMuT xml file.

The algorithms are presented here require the datasets provided [there]() to work with the example configuration files. They should be standalone once the necessary libraries have been installed (see README.md in all specific folders).

# Installation
To install these python scripts with ubuntu, one can simply run the script ubuntu-install.sh as follow:
```shell
source ubuntu-install.sh
```
Your admin password will be asked since some apt packages are required (see bellow).

# Pre installation
To run the python scripts described, the following libraries have to be installed:
  - git, optional, the codes could be directly download from the github website (```sudo apt install git```)
  - python-dev (```sudo apt install python-dev```)
  - pip (```sudo apt install python-pip```)
  - numpy (```pip install numpy [--user]```)
  - scipy (```pip install scipy [--user]```)
  - libhdf5 (```sudo apt install libhdf5-devpython ```, to allow read and write of hdf5 images)
  - unittest2 (```pip install unittest2 [--user]```, manual installation of this dependency for hdf5 read/write)
  - matplotlib (```pip install matplotlib [--user]```, specifically for TimeRegistration)
  - PyWavelets (```pip install PyWavelets [--user]```, specifically for TimeRegistration)
  - cmake-curses-gui (```sudo apt install cmake-curses-gui```, specifically for TimeRegistration)

To run the different scripts, it is necessary to pre-install the following libraries:
  - TGMMlibraries
  - IO
  - BlockMatching

## TGMMlibraries
TGMMlibraries is a class that allows to manipulate lineage trees un python.

To install TGMMlibraries, from a terminal go into the folder TGMMlibraries and run the following command:
```shell
python setup.py install [--user]
```

The option ```[--user]``` is not necessary but allows to install only for a particular user

## IO
IO is a library to read and write 3D images from different formats (tiff, inr, klb, h5).
To install IO, you need to first install pyklb the following way:
```shell
pip install cython [--user]
git clone https://github.com/bhoeckendorf/pyklb.git
cd pyklb
python setup.py bdist_wheel
pip install dist/<name_of_the_wheel>.whl [--user]
```

Then you need to manually link the pyklb libraries to your python environement:
```shell
path_lib='~/.local/lib' #This is an example that should work on linux
ln build/lib.<name-of-version>-2.7/pyklb.so $path_lib #links the two libraries to your libraries folder
ln build/lib/libklb.so $path_lib

echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'$path_lib >> ~/.bashrc #replace .bashrc by .profile for macOs for example
echo 'export PYTHONPATH=$PYTHONPATH:'$path_lib >> ~/.bashrc #replace .bashrc by .profile for macOs for example

source ~/.bashrc
```

Then you can go into the folder IO and run the following command:
```shell
python setup.py install [--user]
```

The option ```[--user]``` is not necessary but allows to install only for a particular user

## BlockMatching
BlockMatching is the C code that allows to register 3D images together. Only Time-registration is dependant on BlockMatching.

To install blockmatching a c/c++ compiler is necessary (gcc for example) together with Makefile, cMake and ccmake tools (see the pre-installation section).

First one have to install the external libraries. To do so, one can run the following commands in a terminal from the BlockMatching folder:
```shell
# To install klb read/write tools
cd external/KLBFILE/
mkdir build
cd build
cmake ../keller-lab-block-filetype
make -j<24> #specify the number of cores that you want to allow for the build

# To install tiff read/write tools, from the folder BlockMatching
cd external/TIFF
mkdir build
cd build
cmake ../tiff-4.0.6
make -j<24> #specify the number of cores that you want to allow for the build
```

Once these are installed one can run the following commands in a terminal from the BlockMatching folder:
```shell
mkdir build
cd build
ccmake ..
# press c then e
# Then enter the correct absolute paths for the tiff and klb builds (ie /path/to/BlockMatching/external/TIFF/build and /path/to/BlockMatching/external/KLBFILE/build)
# Then c then e then g
make -j<24> #specify the number of cores that you want to allow for the build
```

Then this newly built binaries (found in BlockMatching/build/bin) have to be accessible from the different scripts that will be ran. To do so one can add the following command to their \~/.bashrc (\~/.profile for mac users):
```shell
export PATH=$PATH:/path/to/BlockMatching/build/bin
```

One 'direct' way to do so is to run the following command:
```shell
echo 'export PATH=$PATH:/path/to/BlockMatching/build/bin' >> ~/.bashrc
```

or add a line to the csv configuration file for Time-registration/standalone-registration.py (the csv file is Time-registration/csv-parameter-files/standalone-registration.csv).

# Installation of the scripts
To install the scripts one can run the setup.py from each folder. For example, from this folder one can run this sequence of commands:
```shell
cd TGMMlibraries
python setup.py install --user

cd ../IO
python setup.py install --user

cd ../I2AE
python setup.py install --user

cd ../SVF
python setup.py install --user

cd ../Time-registration
python setup.py install --user

cd ../svf2MaMuT
python setup.py install --user
```
