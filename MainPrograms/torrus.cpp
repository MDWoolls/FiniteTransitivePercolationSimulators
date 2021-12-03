#include <itpp/itcomm.h>
//#include <Eigen/Sparse>
//#include <Eigen/Core>

#include <fstream>
//#include <vector>
#include <list>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <climits>
//#include <algorithm>
//#include <string>
//#include <time.h>
//#include <chrono>
#include <stdio.h>

//#include <gmp.h>
//#include <gmpxx.h>
#include <assert.h>

//#include <mpi.h>

//#include "../headers/bigInt.h"
//#include "../headers/bigInt.cpp"
//#include "../headers/Lattice.h"
//#include "../headers/cluster.h"

//typedef Eigen::SparseMatrix<int> SpMati;
//typedef Eigen::Triplet<int> Trip;
//typedef SpMati::InnerIterator Spit;
//typedef std::vector< std::vector<int> > varray;
typedef unsigned long int ulint;

//max L = 1518500249
//max long int = 18446744073709551615

int main(int argc, char **argv){
    //Intializes things from inputs----------------------------------------------------
    itpp::Parser p; 
    p.init(argc,argv);
    if(argc>1){
	if((strcmp(argv[1],"--help")==0)||(strcmp(argv[1],"-h")==0)){
	    std::cout << argv[0] << ":Outputs an torrus as a incidence matrix:\n"
		      << "By default:\n"
		      << "   file=default.mtx\n"
		      << "   L=10\n"
		      << "   dual=false\n"
                      << std::endl;
      	    exit(-1);
	}
    }
    p.set_silentmode(true);

    std::string file="default.mtx";  
    p.get(file,"file");

    bool dual=false;
    p.get(dual,"dual");

    ulint L=10;
    p.get(L,"L");
 
    //std::cout << "size of long int is " << sizeof(unsigned long int) << std::endl;

    //SpMati A=Torrus(L,dual); 
    std::ofstream dat;
    dat.open(file.c_str());
    dat << "%%MatrixMarket matrix coordinate integer general\n"
        << "%%Torrus with L=" << L << " and dual is " << dual << "\n"
        << L*L << " " << 2*L*L << " " << 4*L*L << std::endl;
    
    if(!dual){
  	for(ulint i=0; i<L; i++){
	    for(ulint j=0; j<L; j++){
  		//Diagonal values of R      
      	        //r=i;
    	        //c=i;
    	        dat << L*i+j+1 << " " << L*i+j+1 << " " << 1 << "\n";
	        dat << L*j+i+1 << " " << L*j+i+L*L+1 << " " << 1 << "\n";

    	        //Off Diagonal values of R
    	        //r=(i+1)%L;
    	        //c=i;
	        dat << L*j+(i+1)%L+1 << " " << L*j+i+1 << " " << 1 << "\n";
     	        dat << L*i+j+1 << " " << L*((i+1)%L)+j+L*L+1 << " " << 1 << "\n";
	    }
	    dat.flush();
        }
    }
    else{
	for(ulint i=0; i<L; i++){
  	    //Kronecker Product
      	    //(AxB)[p*r+v][q*s+w]=a[r][s]b[v][w]
    	    //a is mxn b is pxq
    	    //A[L*r+j][L*c+j]=R[r][c]I[j][j]
    	    //A[L*j+r][L*j+c+L^2]=I[j][j]RT[r][c]
    	    for(ulint j=0; j<L; j++){
    	        //Diagonal values of R      
    	        //r=i;
    	        //c=i;
	        dat << L*i+j+1 << " " << L*i+j+L*L+1 << " " << 1 << "\n";
     	        dat << L*j+i+1 << " " << L*j+i+1 << " " << 1 << "\n";

    	        //Off Diagonal values of R
    	        //r=i;
    	        //c=(i+1)%L;
	        dat << L*j+i+1 << " " << L*j+(i+1)%L+L*L+1 << " " << 1 << "\n";
    	        dat << L*((i+1)%L)+j+1 << " " << L*i+j+1 << " " << 1 << "\n";
	    }
	    dat.flush();
        }
    }
    dat.close();
    return 0;
}
