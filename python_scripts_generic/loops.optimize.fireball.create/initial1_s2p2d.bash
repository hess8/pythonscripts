#!/bin/sh  

## Hydrogen
H="H"
n=$1
#cp paramatom parameters1  #for joinparam.x (not joining)
cd ../beginrun
cp initial.run.inputHDN initial.run.input
./initial.x <initial.run.input >initialH.out
mv initial.run.input ../output/initial.run.input$H${n}
cp begin.input begin.input1 # so it's available later to view
./begin.x >beginoutH
gnuplot statesHDN.plt
mv plot.eps plot$H${n}.eps
cp plot$H${n}.eps ../output
rm  plot$H${n}.eps


cd ../loops




