#!/bin/bash

#
# Basic use of blockmatching
#

rm -rf EXAMPLE01
mkdir -p EXAMPLE01

#
# generates a transformed version of input image
#
applyTrsf data/irm2D-1.hdr  \
          EXAMPLE01/irm2D-trsf.hdr  \
          -trsf data/rigid2D.trsf

#
# registration
#
blockmatching -ref data/irm2D-1.hdr \
              -flo EXAMPLE01/irm2D-trsf.hdr  \
              -res EXAMPLE01/res.hdr  \
              -res-trsf EXAMPLE01/res.trsf \

#
# registration results
#
# EXAMPLE01/res.hdr is comparable to data/irm2D-1.hdr
# EXAMPLE01/res.trsf is comparable to the inverse of data/rigid2D.trsf
# 

#
# Once the transformation has been computed, the floating image
# can be resampled again with applyTrsf
#

applyTrsf EXAMPLE01/irm2D-trsf.hdr \
          EXAMPLE01/res-2.hdr\
          -trsf EXAMPLE01/res.trsf
