#!/bin/sh  

# Uses etot and xcomp (geometry) to minimize.  s2p2d   Fix d parameters for now
# Scale of xcomp variation is 0.01 ang.  Scale of etot variation is 0.1 eV. 
#  Scale each to about 1.
# So make the fit error:  diff2/1e-2 + (etot + 86.1)/0.1
#  Loop to optimize

###  All bash modules are started in the optimize/loops folder  Do everything relative to this

######################## Initialize


# Clean out files. 
curdir=$PWD  #loop dir
rm optim*.e*
rm optim*.o*
#pwd
rm ../output/*

## Make everything
cd $curdir #loops/
#cd ../../progs
#make
cd $curdir #loops/
cd ../../cprogs
#make
cd $curdir #loops/
cd ../../begin_2007  
make all

cd $curdir #loops/

##### Initialize Begin
cd ../beginrun
rm *


cp ../input/begininput/*.* ../beginrun
ln -s ../../begin_2007/initial.x initial.x
ln -s ../../begin_2007/begin.x begin.x


##### Initialize Create
cd ../createrun
pwd
rm *
rm coutput/*
rm cin/3splt/*
cp ../input/createinput/* ../createrun
ln -s ../../cprogs/create.x create.x
cd ..
pwd

##### Initialize Fireball
rm Fdata
ln -s createrun/coutput/ Fdata 
cd fireballrun
pwd
rm fireball.x
ln -s ../../progs/fireball.x fireball.x
cp ../input/fireballinput/* ../fireballrun
cp ../fireballrun/script.input ../fireballrun2/
cp  ../diamond1/script.input ../diamond2/ 
rm CHARGES
cp benzene.bas answer.bas


cd ../loops
pwd
