#!/bin/sh
# This file is subject to the terms and conditions defined in
# file 'LICENCE', which is part of this source code package.
# Author: Leo Guignard (guignardl...@AT@...janelia.hhmi.org)
#
# Run this script in a shell in order to install all the required dependencies
# to be able to run the scripts

sudo apt install git python-dev python-pip libhdf5-dev cmake-curses-gui
pip install numpy scipy unittest2 matplotlib PyWavelets cython scikit-image --user

# git clone https://github.com/leoguignard/standalone-Mouse.git
git clone https://github.com/bhoeckendorf/pyklb.git

cd pyklb
python setup.py bdist_wheel
wheel_name=`ls dist/*.whl`
pip install $wheel_name --user
path_lib=~/.local/lib/ #This is an example that should work on linux
cd build/lib.*
cp pyklb.so $path_lib #links the two libraries to your libraries folder
cd ../..
cp build/lib/libklb.so $path_lib
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'$path_lib >> ~/.bashrc #replace .bashrc by .profile for macOs for example
echo 'export PYTHONPATH=$PYTHONPATH:'$path_lib >> ~/.bashrc #replace .bashrc by .profile for macOs for example

# To install klb read/write tools
cd ../BlockMatching/external/KLBFILE/
mkdir build
cd build
cmake ../keller-lab-block-filetype
make
klb_path=$(pwd | sed 's/\//\\\//g')

# To install tiff read/write tools, from the folder BlockMatching
cd ../../TIFF
mkdir build
cd build
cmake ../tiff-4.0.6
make
tiff_path=$(pwd | sed 's/\//\\\//g')

cd ../../..

mkdir build
cd build

sed_command=s/.path_to_tiff*/$tiff_path/
sed $sed_command < ../CMakeLists_to_fill.txt > tmp
sed_command=s/.path_to_klb*/$klb_path/
sed $sed_command < tmp > ../CMakeLists.txt

cmake ..
make
cd bin

echo 'export PATH=$PATH:'`pwd` >> ~/.bashrc

cd ../../../

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

cd ..
source ~/.bashrc
