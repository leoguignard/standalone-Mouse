#!/bin/bash

#
# Use of -left-transformation option of blockmatching
#

rm -rf EXAMPLE02
mkdir -p EXAMPLE02

#
# registration
# it comes to register (data/irm2D-1.hdr o data/rigid2D.trsf)
# onto data/irm2D-1.hdr
# this is similar to example01.bash
#

blockmatching -ref data/irm2D-1.hdr \
              -flo data/irm2D-1.hdr  \
              -left-transformation data/rigid2D.trsf \
              -res EXAMPLE02/res.hdr  \
              -res-trsf EXAMPLE02/res.trsf

#
# registration results
#
# EXAMPLE02/res.hdr is comparable to data/irm2D-1.hdr
#   while EXAMPLE01/res.hdr was cropped wrt to data/irm2D-1.hdr
#   because of the first resampling by data/rigid2D.trsf
#   here EXAMPLE02/res.hdr isn't cropped.
# EXAMPLE02/res.trsf is comparable to the inverse of data/rigid2D.trsf
# 