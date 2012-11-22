#!/bin/sh  
C="C"
n=$1

cd ../beginrun
## Carbon
#  parameters for V0 r0 are already written to initial.run.input by minfunc python code 

./initial.x <initial.run.input >initialC.out 
mv initial.run.input ../output/initial.run.input$C${n}
cp begin.input begin.input2
./begin.x >beginoutC   
gnuplot statesDNPol.plt
mv plot.eps plot$C${n}.eps
cp plot$C${n}.eps ../output
rm plot$C${n}.eps

## get files to create
cp 0* ../createrun/cin/3splt/  #files with name starting with zero including wavefunctions


## create initialization
cd ../createrun
cp ../input/createinput/* ../createrun  #theory.input, create.input, etc (not folders)

cd ../loops




