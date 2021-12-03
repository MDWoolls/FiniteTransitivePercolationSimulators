#include <itpp/itcomm.h>

#include <iostream>
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
      std::cout<< argv[0] << ": Performs bond percolation on "
	       << "incidence matrix and its dual matrix.\n"
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
	       << "  nsave=3 (number of top largest cluster sizes to keep)\n"
	       << "  NA=10 (How many averages to do)\n"
	       << "  Np=1000 (How many points to save)\n"
	       << "  datfile=default.tsv (File to save data to)\n"
	       << "  dualdatfile=dualdefault.tsv (File to save dual data to)\n"
	       << "  graphfile=testX.mtx (File containing graph)\n"
	       << "  dualfile=testZ.mtx (File containing dual graph)\n"
	       << "  pmin=0 (smallest p-value to calculate)\n"
	       << "  pmax=1 (largest p-value to calculate)\n"
	       << "  scut=-1 (number of sigmas away x can be before it is set to zero. Negative means no cut of\n"
	       << "  seed=-1 (if seed is less than 0 then time(NULL) is used\n"
	       << "  silentmode=true (sets silent mode for parser)\n"
	       << "The graphfile and dualfile must be duals "
	       << "and incidence matricies.\n"
	       << std::endl;
      exit (-1);
    }
  }
  
  p.set_silentmode(true);
  bool silentmode=true;
  p.get(silentmode,"silentmode");
  
  if(!silentmode) p.set_silentmode(false);
  
  std::string datfile, dualdatfile, graphfile, dualfile;
  datfile="default.tsv";
  dualdatfile="dualdefault.tsv";
  graphfile="testX.mtx";
  dualfile="testZ.mtx";
  p.get(datfile,"datfile");
  p.get(graphfile,"graphfile");
  p.get(dualfile,"dualfile");
  p.get(dualdatfile,"dualdatfile");

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
  /*
  int id=0,P=1;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&id);
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  */
  boost::mpi::environment env(argc,argv);
  boost::mpi::communicator world;

  if(seed<0) rng.seed(time(NULL));
  else rng.seed(seed);

  for(int i=0;i<world.rank();i++){
      rng.jump();
  }
  
  //initializes bonds from files
  cluster c,cd;//c0,cd0;
  
  //bonds contains 2 lists the 2nd is dual
  std::vector< std::vector<lint> > bonds;
  bonds.resize(2);
  
  //c0.initialize(bonds[0],graphfile,nsave);
  //cd0.initialize(bonds[1],dualfile,nsave);
  c.initialize(bonds[0],graphfile,nsave);
  cd.initialize(bonds[1],dualfile,nsave);

  if(!silentmode){
    std::cout << world.rank() << ": Initialized cluster\n";
    c.output();
    cd.output();
    std::cout << "bonds=";
    for(const auto& b : bonds[0]){
      std::cout << b << " ";
    }
    std::cout << std::endl;

    cd.output();
    std::cout << "bonds=\n";
    for(const auto& b : bonds[1]){
      std::cout << b << " ";
    }
    std::cout << std::endl;
  }

  ulint Nb=c.Nb;

  //Initializing data--------------------------------------------------------------------------   
  //std::cout << "Initializing data" << std::endl;
  std::vector<std::vector<mpf_class> > sum, sumd; //where data is stored
  std::vector<ulint> k(Nb+1),kd(Nb+1);
  sum.reserve(Nb+1);
  sumd.reserve(Nb+1);
  std::vector<mpf_class> v0(c.nsave*4+6);
   
  for(ulint i=0; i<=Nb; i++){
      v0[0]=i;
      sum.push_back(v0);
      sumd.push_back(v0);
  }
 
  //Loops over averages----------------------------------------------------------
  for(long int l=world.rank();l<NA;l+=world.size()){
    //if(!silentmode) std::cout << world.rank() << ": l=" << l << "\n";
    //randomizes the order of bonds
    permutation(bonds);
    
    //Percolation for graph------------------------------------------------------
    //c=c0;
    //cd=cd0;
    c.reset();
    cd.reset();

    c.UpdateData(sum[0]);
    cd.UpdateData(sumd[0]);
        
    k[0]=c.nclst;
    kd[0]=cd.nclst;

    //loops over each bond
    for(ulint i=0; i<Nb; i++){
      //Adds bond to clst
      c.addbond(bonds[0][i],bonds[0][i+Nb]);
      cd.addbond(bonds[1][Nb-i-1],bonds[1][Nb-1-i+Nb]);
      
      //adds data to sum
      c.UpdateData(sum[i+1]);
      cd.UpdateData(sumd[i+1]);
      k[i+1]=c.nclst;
      kd[i+1]=cd.nclst;
    }

    c.UpdateHomologies(sum,k,kd);
    cd.UpdateHomologies(sumd,kd,k);
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

  //outputs graph to file after reaveraging it---------------------------------------------------  
  OutputData(sum,c.Ns,c.Nb,NA,c.nsave,datfile,pmin,pmax,Np,world,scut);
  //outputs dual graph to file after reaverging it---------------------------------------------------  
  //std::cout << "id=" << id << " dual outputting data" << std::endl;
  OutputData(sumd,cd.Ns,cd.Nb,NA,cd.nsave,dualdatfile,pmin,pmax,Np,world,scut);
  
  //MPI_Finalize();
  return 0;
}   


