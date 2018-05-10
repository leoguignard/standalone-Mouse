#!/bin/bash

#
# Use of -initial-result-transformation option of blockmatching
#

rm -rf EXAMPLE03
mkdir -p EXAMPLE03

#
# generates a transformed version of input image
#
applyTrsf data/irm2D-1.hdr  \
          EXAMPLE03/irm2D-trsf.hdr  \
          -trsf data/rigid2D.trsf

#
# registration
# only 2 iterations are performed at each level
# in order not to reach convergence
#
# -flo-frac 0.75 is required to ensure that the 
# same fraction of floating blocks are used 
# in both computation schemes
#

blockmatching -ref data/irm2D-1.hdr \
              -flo  EXAMPLE03/irm2D-trsf.hdr  \
              -res EXAMPLE03/res-1.hdr  \
              -res-trsf EXAMPLE03/res-1.trsf \
              -py-hl 3 -py-ll 2 \
              -flo-frac 0.75 \
              -max-iterations 2 

blockmatching -ref data/irm2D-1.hdr \
              -flo  EXAMPLE03/irm2D-trsf.hdr  \
              -res EXAMPLE03/res-2.hdr  \
              -init-res-trsf EXAMPLE03/res-1.trsf \
              -res-trsf EXAMPLE03/res-2.trsf \
              -py-hl 1 -py-ll 0 \
              -flo-frac 0.75 \
              -max-iterations 2 

blockmatching -ref data/irm2D-1.hdr \
              -flo  EXAMPLE03/irm2D-trsf.hdr  \
              -res EXAMPLE03/res-3.hdr  \
              -res-trsf EXAMPLE03/res-3.trsf \
              -py-hl 3 -py-ll 0 \
              -flo-frac 0.75 \
              -max-iterations 2 

#
# registration results
#
# EXAMPLE03/res-2.trsf and EXAMPLE03/res-3.trsf are equal
# (up to numerical "noise", e.g. writing/reading transformation
# matrices) as well as EXAMPLE03/res-2.hdr and EXAMPLE03/res-3.hdr
# 