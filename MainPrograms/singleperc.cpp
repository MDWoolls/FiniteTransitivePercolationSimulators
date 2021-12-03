#include <itpp/itcomm.h>

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <stdio.h>

#include <gmp.h>
#include <gmpxx.h>

#include <mpi.h>
#include <boost/mpi.hpp>

#include "../headers/cluster.h"
#include "../headers/functions.h"

//takes in distance to find data for and number of times to average over it
//creates/appends to a file the new data values averaged over for each data type
int main(int argc, char **argv){
  //Intializes things from inputs----------------------------------------------------
  itpp::Parser p; 
  p.init(argc,argv);
  if(argc>1){
    if((strcmp(argv[1],"--help")==0)||(strcmp(argv[1],"-h")==0)){
      std::cout<< argv[0] << ": Performs bond percolation using "
	       << "incidence matrix with out homology calculations.\n"
	       << "It outputs the average of the following information (in each file):\n"
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
	       << "  graphfile=test.mtx (File containing graph)\n"
	       << "  pmin=0 (smallest p-value to calculate)\n"
	       << "  pmax=1 (largest p-value to calculate)\n"
	       << "  scut=-1 (number of sigmas away x can be before it is set to zero. Negative means no cut of\n"
	       << "  seed=-1 (if seed is less than 0 then time(NULL) is used\n"
	       << "  silentmode=true (sets if information is to be outputted)\n"
	       << std::endl;
      exit (-1);
    }
  }
  bool silentmode=true;
  p.set_silentmode(true);
  p.get(silentmode,"silentmode");
  p.set_silentmode(silentmode);

  std::string datfile,graphfile;
  datfile="default.tsv";
  graphfile="test.mtx";
  
  p.get(datfile,"datfile");
  p.get(graphfile,"graphfile");

  int NA=10,nsave=2,Np=101,scut=-1,seed=-1;
  
  p.get(NA,"NA");
  p.get(nsave,"nsave");
  p.get(Np,"Np");
  p.get(scut,"scut");
  p.get(seed,"seed");

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
  
  boost::mpi::environment env(argc,argv);
  boost::mpi::communicator world;

  if(seed<0) rng.seed(time(NULL));
  else rng.seed(seed);

  for(int i=0;i<world.rank();i++){
      rng.jump();
  }
  
  
  //initializes bonds from files
  cluster c;
  
  std::vector<lint> bonds;
  c.initialize(bonds,graphfile,nsave);
  
  ulint Nb=c.Nb;
  
  //Initializing data--------------------------------------------------------------------------   
  std::vector<std::vector<mpf_class> > sum; //where data is stored
  sum.reserve(Nb+1);
  std::vector<mpf_class> v0(c.nsave*4+6);
  
  for(ulint i=0; i<=Nb; i++){
    v0[0]=ceil(i);
    sum.push_back(v0);
  }
  
  //Loops over averages----------------------------------------------------------
  for(int l=world.rank();l<NA;l+=world.size()){
    //randomizes the order of bonds
    permutation(bonds);
    //Percolation for graph------------------------------------------------------
    c.reset();   
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
  if(!silentmode) std::cout << world.rank() << " finished calculating data" << std::endl;
  
  //outputs graph to file after reaveraging it---------------------------------------------------  
  OutputData(sum,c.Ns,c.Nb,NA,c.nsave,datfile,pmin,pmax,Np,world,scut,silentmode);
  
  return 0;
}


