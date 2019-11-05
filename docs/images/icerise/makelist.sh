#!/bin/bash
#
# 1 argument, the name of the directory
#

namefile=listplotfiles

basis=/home/lfavier/Programs/bisicles/BISICLES/examples/IceRiseSmallStudy/$1

cd $basis
rm -f $namefile

for i in plot.iceRise.*
do
  echo $i >> $namefile
done

