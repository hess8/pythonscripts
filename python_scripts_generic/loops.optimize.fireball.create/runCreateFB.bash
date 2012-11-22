
############################## Run programs #######################################################
n=$1  
echo "run $n" 
##  
read NP < NP.info
read PBS_NODEFILE < pbs.info
######## Create
cd ../createrun
rm coutput/* # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
date
#
#/opt/mpich/gnu/bin/mpirun \  #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
#        -machinefile $PBS_NODEFILE \
#        -np $NP \
#	create.x >output.log #<create.run.input
#
#sleep 100
#date

#./create.x >output.log <create.run.input
/opt/mpiexec/bin/mpiexec  \
	create.x < create.run.input > output.log

######## Fireball
cd ..
pwd
rm Fdata
ln -s createrun/coutput/ Fdata 
cd fireballrun
pwd
rm CHARGES
rm fb  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!
cp benzene.bas answer.bas
./fireball.x >fb   #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cp Eparts.dat "Eparts1Ben${n}.dat"
mv  "Eparts1Ben${n}.dat" ../output
cp answer.bas "answer1Ben${n}.xyz"
mv   "answer1Ben${n}.xyz" ../output
cp CHARGES  "../output/CHARGES1Ben${n}"
grep prev fb > ../output/grepT1Ben${n}
cd ../loops


######## Fireball2
cd ../fireballrun2
pwd
rm CHARGES
rm fb
cp benzene.bas answer.bas
./fireball.x >fb   #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cp Eparts.dat "Eparts2Ben${n}.dat"
mv  "Eparts2Ben${n}.dat" ../output
cp answer.bas "answer2Ben${n}.xyz"
mv   "answer2Ben${n}.xyz" ../output
cp CHARGES  "../output/CHARGES2Ben${n}"
grep prev fb > ../output/grepT2Ben${n}
cd ../loops
