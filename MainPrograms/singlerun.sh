#!/bin/sh -l

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem-per-cpu=7G
#SBATCH --mail-user=mwool001@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -p batch
#SBATCH --time=30-00:00:00

source ../LoadModules.sh

df="./SingleData/FILL/"


Dim=$1
i=$2

#np=20
#dL=5
Lmin=$3
Lmax=$4

for (( L=Lmin;L<=Lmax;L++ ))
do
    dat=${df/FILL/$i}    
    dat+="hypercube_${Dim}_${L}.tsv"
    graph="tmp_${Dim}_${L}_${i}.mtx"
    
    echo $graph
    echo $dat
    
    ./hypercube file=$graph L=$L D=$Dim
    mpirun -n 64 singleperc.out NA=1000000 datfile=$dat Np=1001 scut=100 graphfile=$graph

    rm $graph

done

echo "-DONE-"
