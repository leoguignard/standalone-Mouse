#!/bin/bash

#
# Use of -initial-result-transformation option of blockmatching
# test with vector field transformations
#

rm -rf EXAMPLE04
mkdir -p EXAMPLE04

#
# generates a transformed version of input image
#
applyTrsf data/irm2D-1.hdr  \
          EXAMPLE04/irm2D-trsf.hdr  \
          -trsf data/vector2D.trsf

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
              -flo  EXAMPLE04/irm2D-trsf.hdr  \
              -res EXAMPLE04/res-1.hdr  \
              -res-trsf EXAMPLE04/res-1.trsf \
              -py-hl 3 -py-ll 2 \
              -flo-frac 0.75 \
              -trsf-type vectorfield2D \
              -max-iterations 10 

blockmatching -ref data/irm2D-1.hdr \
              -flo  EXAMPLE04/irm2D-trsf.hdr  \
              -res EXAMPLE04/res-2.hdr  \
              -init-res-trsf EXAMPLE04/res-1.trsf \
              -res-trsf EXAMPLE04/res-2.trsf \
              -py-hl 1 -py-ll 0 \
              -flo-frac 0.75 \
              -trsf-type vectorfield2D \
              -max-iterations 10 

blockmatching -ref data/irm2D-1.hdr \
              -flo  EXAMPLE04/irm2D-trsf.hdr  \
              -res EXAMPLE04/res-3.hdr  \
              -res-trsf EXAMPLE04/res-3.trsf \
              -py-hl 3 -py-ll 0 \
              -flo-frac 0.75 \
              -trsf-type vectorfield2D \
              -max-iterations 10 

#
# registration results
#
# EXAMPLE04/res-2.trsf and EXAMPLE04/res-3.trsf are equal
# (up to numerical "noise", e.g. writing/reading transformation
# matrices) as well as EXAMPLE04/res-2.hdr and EXAMPLE04/res-3.hdr
# 