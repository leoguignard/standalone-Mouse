

#
# creation of pairings
#

echo "35 35 0" > floating.pts
echo "55 60 0" >> floating.pts

echo "55 40 0" > reference.pts
echo "40 65 0" >> reference.pts

#
# creation of tests images
#

createGrid -dim 100 100 mosaic.mha -type mosaic -spacing 10 10
createGrid -dim 100 100 grid.mha -spacing 10 10
copy -norma mosaic.mha mosaic.mha

#
# vector field computation with pointmatching
#

pointmatching -flo floating.pts -ref reference.pts -trsf-type vectorfield \
  -template mosaic.mha \
  -propagation-distance 0 \
  -fading-distance 0 \
  -fluid-sigma 5 \
  -res-trsf vec-00-00-05.trsf

applyTrsf mosaic.mha mosaic-00-00-05.mha -trsf vec-00-00-05.trsf -interpolation nearest
applyTrsf grid.mha grid-00-00-05.mha -trsf vec-00-00-05.trsf



pointmatching -flo floating.pts -ref reference.pts -trsf-type vectorfield \
  -template mosaic.mha \
  -propagation-distance 10 \
  -fading-distance 0 \
  -fluid-sigma 0 \
  -res-trsf vec-10-00-00.trsf

applyTrsf mosaic.mha mosaic-10-00-00.mha -trsf vec-10-00-00.trsf -interpolation nearest
applyTrsf grid.mha grid-10-00-00.mha -trsf vec-10-00-00.trsf



pointmatching -flo floating.pts -ref reference.pts -trsf-type vectorfield \
  -template mosaic.mha \
  -propagation-distance 10 \
  -fading-distance 0 \
  -fluid-sigma 5 \
  -res-trsf vec-10-00-05.trsf

applyTrsf mosaic.mha mosaic-10-00-05.mha -trsf vec-10-00-05.trsf -interpolation nearest
applyTrsf grid.mha grid-10-00-05.mha -trsf vec-10-00-05.trsf



pointmatching -flo floating.pts -ref reference.pts -trsf-type vectorfield \
  -template mosaic.mha \
  -propagation-distance 00 \
  -fading-distance 10 \
  -fluid-sigma 0 \
  -res-trsf vec-00-10-00.trsf

applyTrsf mosaic.mha mosaic-00-10-00.mha -trsf vec-00-10-00.trsf -interpolation nearest
applyTrsf grid.mha grid-00-10-00.mha -trsf vec-00-10-00.trsf



pointmatching -flo floating.pts -ref reference.pts -trsf-type vectorfield \
  -template mosaic.mha \
  -propagation-distance 00 \
  -fading-distance 10 \
  -fluid-sigma 5 \
  -res-trsf vec-00-10-05.trsf

applyTrsf mosaic.mha mosaic-00-10-05.mha -trsf vec-00-10-05.trsf -interpolation nearest
applyTrsf grid.mha grid-00-10-05.mha -trsf vec-00-10-05.trsf



pointmatching -flo floating.pts -ref reference.pts -trsf-type vectorfield \
  -template mosaic.mha \
  -propagation-distance 10 \
  -fading-distance 10 \
  -fluid-sigma 0 \
  -res-trsf vec-10-10-00.trsf

applyTrsf mosaic.mha mosaic-10-10-00.mha -trsf vec-10-10-00.trsf -interpolation nearest
applyTrsf grid.mha grid-10-10-00.mha -trsf vec-10-10-00.trsf



pointmatching -flo floating.pts -ref reference.pts -trsf-type vectorfield \
  -template mosaic.mha \
  -propagation-distance 10 \
  -fading-distance 10 \
  -fluid-sigma 5 \
  -res-trsf vec-10-10-05.trsf

applyTrsf mosaic.mha mosaic-10-10-05.mha -trsf vec-10-10-05.trsf -interpolation nearest
applyTrsf grid.mha grid-10-10-05.mha -trsf vec-10-10-05.trsf






for p in mosaic*.mha grid*.mha
do
    copy $p `basename $p .mha`.pgm
    convert `basename $p .mha`.pgm `basename $p .mha`.png
    rm `basename $p .mha`.pgm
done
