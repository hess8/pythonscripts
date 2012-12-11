#!/bin/sh  

#file='optimize1s/input/begininput/initial.run.input.sp.splitPol'
#file='optimize1s/input/begininput/initial.run.inputHDN'
#file='optimize1s/input/createinput/C.3splitinput'
#file='optimize1s/input/createinput/H.3splitinput'
#file='optimize/input/begininput/statesHDN.plt'
#file=optimize2/beginrun/statesDNPol.plt

#Loop

 
n=622890 #starts at this
nmax=622933
limit=`expr $nmax + 1`
while [ "$n" -lt "$limit" ]  ############## Loop
do
  echo "job $n"
  qdel $n
  n=`expr $n + 1`

done  ######################### end of loop


