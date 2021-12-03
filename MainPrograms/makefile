opt = -O3
pkgs = -L/opt/linux/centos/7.x/x86_64/pkgs
itpp = ${pkgs}/itpp/4.3.1/lib/ -litpp
linker = -lm
debug = -O0
boost = -L/opt/linux/centos/7.x/x86_64/pkgs/boost/1.70.0/lib -I/opt/linux/centos/7.x/x86_64/pkgs/boost/1.70.0/include  -L/opt/linux/centos/7.x/x86_64/pkgs/boost/1.70.0/lib/libboost_mpi.so.1.70.0 -L/opt/linux/centos/7.x/x86_64/pkgs/boost/1.70.0/lib/libboost_math_c99.so.1.70.0 -lboost_mpi -lboost_serialization

CXX = g++ -Wall -g -m64 -std=c++0x ${opt} ${itpp} #${debug}

MPI = mpic++ -Wall -g -m64 -std=c++0x ${opt} ${itpp} #${debug}

all: torrus, hypercube, percolate, singleperc, randomperc

torrus: torrus.cpp
	${CXX} -o torrus $< -litpp -lm

hypercube: hypercube.cpp
	${CXX} -o hypercube $< -litpp -lm

percolate: percolate.cpp ../headers/cluster.h ../headers/functions.h makefile
	${MPI} -o percolate.out  $< -litpp -lm -lgmpxx -lgmp ${boost}

singleperc: singleperc.cpp ../headers/cluster.h ../headers/functions.h makefile
	${MPI} -o singleperc.out $< -litpp -lm -lgmpxx -lgmp ${boost}

randomperc: randomperc.cpp ../headers/functions.h ../headers/cluster.h ../headers/random.h makefile
	${MPI} -o randomperc.out $< -litpp -lm -lgmpxx -lgmp
