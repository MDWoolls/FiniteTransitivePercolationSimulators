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
//#include <stdio.h>
#include <stdio.h>

#include <gmp.h>
#include <gmpxx.h>
#include <assert.h>

#include <mpi.h>

//#include "../headers/bigInt.h"
//#include "../headers/bigInt.cpp"
//#include "../headers/Lattice.h"
#include "../headers/cluster.h"
#include "../headers/functions.h"
//#include "../headers/random.h"

//typedef Eigen::SparseMatrix<int> SpMati;
//typedef Eigen::Triplet<int> Trip;
//typedef SpMati::InnerIterator Spit;
typedef std::vector< std::vector<int> > varray;
typedef unsigned long int ulint;
typedef long int lint;

//be is the sum of all the bond counts where a a homology first appears
void OutputData(std::vector< std::vector<mpf_class> > &sum,
	const int Ns,const int Nb,const int NA,const int n,
	const std::string datfile,
	const long double pmin,	const long double pmax,	const int Np,
	const int id,const int P,const int scut
	);

//takes in distance to find data for and number of times to average over it
//creates/appends to a file the new data values averaged over for each data type
int main(int argc, char **argv){
  /*
  std::ofstream timing;
  timing.open("before.tsv",std::ios::app);
  itpp::Real_Timer clock;
  //clock.tic();
  */

  //Intializes things from inputs----------------------------------------------------
  itpp::Parser p; 
  p.init(argc,argv);
  if(argc>1){
    if((strcmp(argv[1],"--help")==0)||(strcmp(argv[1],"-h")==0)){
      std::cout<< argv[0] << ": Performs bond percolation on a degree regular graph "
	       << "with out homology calculations.\n"
	       << "It outputs the average of the following information:\n"
	       << "  p-value\n"
	       << "  Number of clusters\n"
	       << "  Size of largest cluster\n"
	       << "  Size of 2nd largest cluster\n"
	       << "  Size of 3rd largest cluster\n"
	       << "  ...Size of nth largest cluster\n"
	       << "  Ratio of 2nd/1st largest clusters\n"
	       << "  Ratio of 3rd/1st largest clusters\n"
	       << "  Ratio of nth/1st largest clusters\n"
	       << "  Erasure probability\n"
	       << "  Number of homologies\n"
	       << "  Averge cluster size of of site 0\n"
	       << "  The square of all previous values except Erasure probability\n"
	       << "  Bound count value such that a homology first occurs\n"
	       << "With defaults:\n"
	       << "  nsave=2 (number of top largest cluster sizes to keep)\n"
	       << "  NA=10 (How many averages to do)\n"
	       << "  Np=101 (How many points to save)\n"
	       << "  datfile=default.tsv (File to save data to)\n"
	       << "  N=10 (Number of sites)\n"
	       << "  deg=3 (Degree of each site)\n"
	       << "  pmin=0 (smallest p-value to calculate)\n"
	       << "  pmax=1 (largest p-value to calculate)\n"
	       << "  scut=-1 (number of sigmas away x can be before it is set to zero. Negative means no cut of)\n"
	       << "  seed=-1 (if seed is less than 0 then time(NULL) is used)\n"
	       << std::endl;
      exit (-1);
    }
  }
  p.set_silentmode(true);
  std::string datfile;
  datfile="default.tsv";
  
  p.get(datfile,"datfile");

  int scut=-1,seed=-1;
  size_t NA=10,nsave=2,Np=101,N=10,deg=3;
  
  p.get(NA,"NA");
  p.get(nsave,"nsave");
  p.get(Np,"Np");
  p.get(scut,"scut");
  p.get(seed,"seed");
  p.get(N,"N");
  p.get(deg,"deg");

  long double pmin=0,pmax=1;

  p.get(pmin,"pmin");
  p.get(pmax,"pmax");
    
  if(pmin<0||pmin>1||pmax<0||pmax>1||pmax<pmin){
      std::cout << "pmin or pmax is out of range\n"
	        << "pmin = " << pmin << "\n"
		<< "pmax = " << pmax << "\n";
      throw;
  }
  //Setting up problem-----------------------------------------------------------------
  //std::cout << "Setting Up" << std::endl;
  int id=0,P=1;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&id);
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  
  if(seed<0) rng.seed(time(NULL));
  else rng.seed(seed);

  for(int i=0;i<id;i++){
      rng.jump();
  }
  
  
  //initializes bonds from files
  cluster c,c0;
  
  //std::cout << "initializing" << std::endl;
  std::vector<size_t> bonds;
  c0.randominitialize(N,deg,nsave);
  
  size_t Nb=c0.Nb;

  //Initializing data--------------------------------------------------------------------------   
  //std::cout << "Initializing data" << std::endl;
  std::vector<std::vector<mpf_class> > sum; //where data is stored
  //std::vector<ulint> k(Nb+1),kd(Nb+1);
  sum.reserve(Nb+1);
  std::vector<mpf_class> v0(c0.nsave*4+6);
   
  for(ulint i=0; i<=Nb; i++){
      v0[0]=ceil(i);
      sum.push_back(v0);
  }
 
  //Loops over averages----------------------------------------------------------
  for(int l=id;l<NA;l+=P){
    //std::cout << "l=" << l << std::endl;
    //creates a new random graph
    bonds=RandomDegreeRegularGraph(N,deg);
    //Percolation for graph------------------------------------------------------
    c=c0;
    c.UpdateData(sum[0]);
    
    //loops over each bond
    for(ulint i=0; i<Nb; i++){
      //Adds bond to clst
      c.addbond(bonds[i],bonds[i+Nb]);
      
      //adds data to sum
      c.UpdateData(sum[i+1]);
    }
  }
  
  //sum contains the following sums in its list (n is nsave):
  //(0) N_bs
  //(1) N_c
  //(2) S_1 
  //(3) S_2
  //(4) S_3
  //    ...
  //(n+1) S_n
  //(n+2) S_2/S_1
  //(n+3) S_3/S_1
  //...
  //(2n) S_n/S_1
  //(2n+1) Pr(N_H>0)
  //(2n+2) N_H
  //(2n+3) N_c^2
  //(2n+4) S_1^2
  //(2n+5) S_2^2
  //(2n+6) S_3^2
  //...
  //(3n+3) S_n^2
  //(3n+4) (S_2/S_1)^2
  //(3n+5) (S_3/S_1)^2
  //...
  //(4n+2) (S_n/S_1)^2
  //(4n+3) N_H^2
  //(4n+4) |C(0)|
  //(4n+5) |C(0)|^2
  /*
void OutputData(std::vector< std::vector<mpf_class> > &sum,
	const int Ns,const int Nb,const int NA,const int n,
	const std::string datfile,
	const long double pmin,	const long double pmax,	const int Np,
	int id=0,int P=1
	){
    */

  //std::cout << "id=" << id << " non-dual outputting data" << std::endl;
  //outputs graph to file after reaveraging it---------------------------------------------------  
  OutputData(sum,c.Ns,c.Nb,NA,c.nsave,datfile,pmin,pmax,Np,id,P,scut);
  //outputs dual graph to file after reaverging it---------------------------------------------------  
  
  MPI_Finalize();
  //std::cout << "Done" << std::endl;
  return 0;
}

//be is the sum of all the bond counts where a a homology first appears
void OutputData(std::vector< std::vector<mpf_class> > &sum,
	const int Ns,const int Nb,const int NA,const int n,
	const std::string datfile,
	const long double pmin,	const long double pmax,	const int Np,
	const int id=0,const int P=1,const int scut=-1
	){
    //outputs graph to file---------------------------------------------------  
    //std::cout << "outputting to file" << std::endl;
    std::ofstream dat;
    long double local=0,global=0;

    //Sums and share data across processors
    //std::cout << "id " << id << " sharing data" << std::endl;
    for(ulint i=0; i<sum.size(); i++){
	for(size_t k=1; k<sum[i].size(); k++){
	    //std::cout << "(" << id << "," << i << "," << k << ")" << std::endl;
    	    local=sum[i][k].get_d()/NA;
	    MPI_Allreduce(&local,&global,1,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	    sum[i][k]=mpf_class(std::to_string(global));
	}
    }

    //reaverages the data
    //std::cout << "id " << id << " doing reaverage" << std::endl;
    std::vector<std::vector<long double> > pdat;
    pdat=reaverage(sum,pmin,pmax,Np,id,P,scut);
    
    //opens and clears file and outputs header information
    if(id==0){
       	dat.open(datfile.c_str(),std::ofstream::trunc);
	if(dat.is_open()){
    	    dat << "#Ns\tNb\tNA\n"
    		<< Ns << "\t"	    
    		<< Nb << "\t"
    		<< NA << "\n"
    		<< "#"
    		<< "(1)p-value\t"
    		<< "(2)N_c\t";
    	    for(int i=0; i<n; i++){
    		dat << "(" << i+3 << ")S_" << i+1 << "\t";
    	    }
    	    for(int i=0;i<n-1;i++){
    		dat << "(" << n+3+i << ")s_" << i+2 << "/s_1" << "\t";
    	    }
    	    dat << "(" << 2*n+2 << ")Pr(N_H>0)" << "\t"
    		<< "(" << 2*n+3 << ")N_H" << "\t"
    		<< "(" << 2*n+4 << ")N_c^2" << "\t";
    	    for(int i=0; i<n; i++){
    		dat << "(" << 2*n+5+i << ")S_" << i+1 << "^2" << "\t";
    	    }
    	    for(int i=0;i<n-1;i++){
    		dat << "(" << 3*n+5+i << ")(s_" << i+2 << "/s_1)^2" << "\t";
    	    }
    	    dat << "(" << 4*n+4 << ")N_H^2" << "\t"
    		<< "(" << 4*n+5 << ")|C(0)|" << "\t"
    		<< "(" << 4*n+6 << ")|C(0)|^2" << std::endl;
	}
	else{
	    std::cout << "id 0 could not open file " << datfile << std::endl;
	    throw;
	}
	dat.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    dat.open(datfile.c_str(),std::ofstream::app);

    //std::cout << "id " << id << " outputing data" << std::endl;
    if(dat.is_open()){
    	//Outputting data to file
    	for(int i=0;i<Np;i++){
	    if(id==i%P){
    		for(size_t j=0;j<pdat[i/P].size(); j++){
		    dat << std::scientific << std::setprecision(8) << pdat[i/P][j] << "\t";
    		}
		dat << std::endl;
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
	}
    }
    else{
    	std::cout << "id " << id << " could not openfile " << datfile << std::endl;
    	throw;
    }
    dat.close();
    return;
}


