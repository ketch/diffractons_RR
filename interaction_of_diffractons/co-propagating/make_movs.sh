#!/bin/bash

clear
mkdir animations
echo "*****************************"
echo "*** Script to create MOVs ***"
echo "*****************************"
echo " "

echo "sigma"
echo "  pcolor"
path="./_plots_to_animation/pcolor"
convert -delay 5 -loop 0 "$path/*.png" -delay 200 "$path/*500*" "./animations/co-propagation_pcolor.mov"
echo "  slices"
path="./_plots_to_animation/slices"
convert -delay 5 -loop 0 "$path/*.eps" -delay 200 "$path/*500*" "./animations/co-propagation_slices.mov"
