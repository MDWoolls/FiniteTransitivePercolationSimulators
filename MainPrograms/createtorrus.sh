#!/bin/sh -l

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem-per-cpu=7G
#SBATCH --mail-user=mwool001@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -p batch
#SBATCH --time=30-00:00:00

module load itpp
module load eigen

nmin=20
nmax=220
ns=20;
#nmax=10
#ns=5

for j in 2
do
    for ((i=ns-2;i<=ns;i++))
    do
	L=$(( i*(nmax-nmin)/ns+nmin ))
	echo $L
	fX="testX.mtx"
	./torrus L=$L file=$fX dual=false
	
	fZ="testZ.mtx"
	./torrus L=$L file=$fZ dual=true
	
	ofX="./Data/${j}/4_4_${L}_X.mtx"
	ofZ="./Data/${j}/4_4_${L}_Z.mtx"
	echo $ofX
	echo $ofZ
	
	mpirun -n 64 percolate.out NA=1000000 graphfile=$fX dualfile=$fZ datfile=$ofX dualdatfile=$ofZ Np=1001 scut=100
    done
done

echo "-DONE-"

