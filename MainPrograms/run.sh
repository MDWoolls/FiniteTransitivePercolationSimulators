#!/bin/sh -l

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem-per-cpu=7G
#SBATCH --mail-user=mwool001@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -p batch
#SBATCH --time=20-00:00:00

module load itpp
module load eigen

#gf="./Q/"
#df="./PentagonData/"

gf="./Graphs/"
gfd="./DualGraphs/"

df="./Data/FILL/"
dfd="./DualData/FILL/"

#p=4
#q=6

#fX=./Graphs/5_6_34560_X.mtx
#fZ=./DualGraphs/6_5_34560_Z.mtx

#dati=./Data/3/5_6_34560.tsv
#datduali=./DualData/3/6_5_34560.tsv

#echo $fX
#echo $fZ
#echo $dati
#echo $datduali

#mpirun -n 64 percolate.out NA=1000000 graphfile=$fX dualfile=$fZ datfile=$dati dualdatfile=$datduali Np=1001 scut=100

for fX in ${gf}*
do
    #fZ=${fX/${p}_${q}/${q}_${p}}
    fZ=${fX:0:9}${fX:11:1}_${fX:9:1}${fX:12}
    fZ=${fZ/X/Z}
    fZ=${fZ/$gf/$gfd}

    dat=${fX/$gf/$df}
    dat=${dat/.mtx/.tsv}

    datdual=${fZ/$gfd/$dfd}
    datdual=${datdual/.mtx/.tsv}
    
    echo $fX
    echo $fZ
    if [[ ( -f "$fX")&&( -f "$fZ") ]]; then
	echo ""
	for i in 1 2 3
	do
	    datduali=${datdual/FILL/$i}
	    dati=${dat/FILL/$i}
	    echo $dati
	    echo $datduali
	    seed=$(date +%s)
	    #echo $seed
	    
    	    mpirun -n 64 percolate.out NA=1000000 graphfile=$fX dualfile=$fZ datfile=$dati dualdatfile=$datduali Np=1001 scut=100 seed=$seed
	done
	rm $fX
	rm $fZ
    else
	echo "$fX or $fZ do not exist"
    fi
    echo ""
done

#gf="./5_5/"

#df="./5_5_data/"

#for fX in ${gf}*_X.mtx
#do
#    fZ=${fX/_X/_Z}
#    dat=${fX/$gf/$df}
#    dat=${dat/.mtx/.tsv}
#    datdual=${fZ/$gf/$df}
#    datdual=${datdual/.mtx/.tsv}
#    
#    echo $fX
#    echo $fZ
#    echo $dat
#    echo $datdual
#    
#    if [[ ( -f "$fX")&&( -f "$fZ") ]]; then
#	echo ""
#	mpirun -n 60 percolate.out NA=1000000 graphfile=$fX dualfile=$fZ datfile=$dat dualdatfile=$datdual Np=1000
#    else
#	echo "$fX or $fZ do not exist"
#    fi
#done

#gf="./SemiHyperGraphsAll/"
#df="./SemiHyperGraphsData/"

#for fX in ${gf}semi_hyper_dual*
#do
#    fZ=${fX/semi_hyper_dual/semi_hyper}
#    dat=${fX/$gf/$df}
#    dat=${dat/.mtx/.tsv}
#    datdual=${fZ/$gf/$df}
#    datdual=${datdual/.mtx/.tsv}
    
#    echo $fX
#    echo $fZ
#    echo $dat
#    echo $datdual
    
#    if [[ ( -f "$fX")&&( -f "$fZ") ]]; then
#	echo ""
#	mpirun -n 60 percolate.out NA=100000 graphfile=$fX dualfile=$fZ datfile=$dat dualdatfile=$datdual Np=1000
#    else
#	echo "$fX or $fZ do not exist"
#    fi
#done

#gf="./Q/"
#df="./PentagonHistograms/"
#gf="./4_5/"
#df="./4_5_Histograms/"

#gf="./SemiHyperGraphsAll/"
#df="./SemiHyperGraphsData/"

#for f in ${gf}*
#do
    #For {5,5} lattice
    #dat=${f/$gf/$df}
    #dat=${dat/.mtx/.tsv}
    #dat=${dat/Q/hist}
    
#    echo $f
#    echo $dat
    
#    if [[ ( -f "$f") ]]; then
#	mpirun -n 60 histperc.out NA=1000000 graphfile=$f datfile=$dat Np=1000
#    else
#	echo "$f do not exist"
#    fi
#done

echo "-DONE-"
