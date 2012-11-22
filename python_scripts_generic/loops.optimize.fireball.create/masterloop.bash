#!/bin/sh  
#cp ../../optimize/loops/masterloop.py   #can't propagate this to other runs as parameters are set here.
cp ../../optimize/loops/minfunc_s2p2d.py .
python masterloop.py > out
