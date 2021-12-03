#include <iostream>
#include <cstdlib>
#include <map>
#include <iomanip>
#include <math.h>

#include <gmp.h>
#include <gmpxx.h>
#include "random.h"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/mpi.hpp>
#include <mpi.h>

#ifndef myfunctions
#define myfunctions

//Generates a random graph with N sites and each site as deg bonds
//returns the lists of bonds
std::vector<size_t> RandomDegreeRegularGraph(const size_t N, const size_t deg){
  //test if graph with these parameters are possible
  if(N<=deg){
    std::cout << "deg is too small" << std::endl;
    throw;
  }
  else if((N*deg)%2!=0){
    std::cout << "The number of bonds is not even" << std::endl;
    throw;
  }

  //the site label and its corresponding degree count
  long int Nb=N*deg/2;
  //list of bonds, list of bonds sorted
  std::vector<size_t> bonds(2*Nb);
  std::vector<long int> vec(2*Nb);
  //fills bonds with deg compies of each site
  for(size_t i=0;i<N;i++){
    for(size_t j=0;j<deg;j++){
      //std::cout << "i=" << i << ", j=" << j << std::endl;
      bonds[j+i*deg]=i;
    }
  }

  //bond being tested,site 1, site 2 keep s1<s2
  long int s1,s2;
  //true if the right graph is found
  bool complete=false;
  //loops until broken out of
  while(!complete){
    //std::cout << "Staring loop" << std::endl;
    complete=true;
    //resets vec
    std::fill(vec.begin(),vec.end(),-1);
    //randomly permutes bonds
    for(long int i=0;i<2*Nb-1;i++){
      s1=rng.RandomInRange(i,2*Nb-1);
      if(s1!=i){
	std::swap(bonds[i],bonds[s1]);
      }
    }

    for(long int k=0;k<Nb;k++){
      if(!complete) break;
      s1=bonds[k];
      s2=bonds[k+Nb];
      if(s1==s2){
	complete=false;
	break;
      }
      else if(s1>s2) std::swap(s1,s2);

      //loops over each site connected to s1
      for(size_t j=0;j<deg;j++){
	if(vec[j+s1*deg]==-1){
	  vec[j+s1*deg]=s2;
	  break;
	}
	else if(vec[j+s1*deg]==s2){
	  complete=false;
	  break;
	}

	if(j==deg-1){
	  std::cout << "More than deg=" << deg << " bonds are being added to site s1=" << s1 << std::endl;
	  throw;
	}
      }
    }
  }
  return bonds;
}

/*
Defines two binomial functions
The direct version calculates the binomial directly
the other version uses recursion and saves previous answers

//stores list of calculated binomials such that Binomial(n,k)=BINOMIALS[n][k]
std::map< unsigned long long int, std::map<unsigned long long int, mpz_class> > BINOMIALS;
//saves list of factorials such that FACTORIALS[n]=n!
std::vector< mpz_class > FACTORIALS = {1,1,2};

//stores list of allready calculated values of n!/k! such that FACTORIALDIV[n][k]=n!/k! and n>=k
std::map< unsigned long long int, std::map< unsigned long long int, mpz_class > > FACTORIALDIV;

//dynamic function which finds the Binomial = n!/(k!(n-k)!)
mpz_class Binomial(size_t n,size_t k)
{
    //checks to make sure k isn't too large
    if(n<k)
    {
        std::cerr << "ERROR: Binomial - n<k" 
		  << "n=" << n << "\n"
		  << "k=" << k << std::endl;
        throw;
    }

    //checks to make sure BINOMIALS already has a map to n
    //creates it if not
    if(BINOMIALS[n].count(k)==0){
      if(n==k) BINOMIALS[n][k]=1;
      else if(k==0) BINOMIALS[n][k]=1;
      else BINOMIALS[n][k]=Binomial(n-1,k-1)+Binomial(n-1,k);
    }

    return BINOMIALS[n][k];
}

//dynamic function which finds n!
mpz_class Factorial(size_t n)
{
  const size_t NMAX=70000;
  if(n>NMAX){
    mpz_class tmp=Factorial(NMAX);
    for(size_t i=NMAX;i<n;i++){
      tmp*=(i+1);
    }
    return tmp;
  }
  //if n is larger than size of FACTORIALS then finds previous factors up to n
  //otherwise just returns FACTORIALS[n]
  for(size_t i=FACTORIALS.size(); i<=n; i++)
    {
      FACTORIALS.push_back(i*FACTORIALS[i-1]);
    }
  
  if(FACTORIALS[n]==0)
    {
      std::cout << "ERROR: Factorial failed and returns 0 at n=" << n << std::endl;
      throw;
    }
  return FACTORIALS[n];
}

//finds n!/k! for n>k
mpz_class FactorialDivision(const size_t n, const size_t k)
{
    if(n<k)
    {
        std::cerr << "Factorial Division - n<k can't return non integer\n" 
		  << "n=" << n << "\n"
		  << "k=" << k  << std::endl;
        throw;
    }
    
    if(FACTORIALDIV[n].count(k)==0){
      if(n==k) FACTORIALDIV[n][k]=1;
      else FACTORIALDIV[n][k]=n*FactorialDivision(n-1,k);
    }

    return FACTORIALDIV[n][k];
}

//function which finds the Binomial = n!/(k!(n-k)!) directly
mpz_class BinomialDirect(size_t n,size_t k)
{
    if(k>n)
    {
        std::cout << "ERROR: k>n = " << k << ">" << n << std::endl;
        throw;
    }
    if(n-k>k)
    {
        return FactorialDivision(n,n-k)/Factorial(k);
    }
    else
    {    
        return FactorialDivision(n,k)/Factorial(n-k);
    }
}

//find p^x where x is a positive integer
mpf_class intpow(const mpf_class p, const size_t x){
    if(x==0) return 1;
    else if(x==1) return p;
    else if(x==2) return p*p;
    mpf_class tmp;
    mpf_pow_ui(tmp.get_mpf_t(),p.get_mpf_t(),x);
    return tmp;
}

//initializes BINOMIALS using BinomialDirect function for a given n
//does the calculation in paralled
void InitializeBinomialList(const size_t nt, const boost::mpi::communicator& w)
{
  //finds the Binomial values
  const size_t NMAX=1000;
  size_t n=nt;
  if(n>NMAX) n=NMAX;

  for(size_t k=w.rank();k<=n;k+=w.size())
    {
      if(BINOMIALS[n].count(k)==0){
	BINOMIALS[n][k]=BinomialDirect(n,k);
      }
      //std::cout << w.rank() << ": BINOMIALS[" << n << "][" << k << "]=" << BINOMIALS[n][k] << std::endl;
    }
  
  std::cout << w.rank() << ": Preparing to share" << std::endl;

  //broadcasts the binomial values
  std::string cast;
  for(size_t i=0;i<=n;i++){
    cast=BINOMIALS[n][i].get_str();
    broadcast(w,cast,i%w.size());
    BINOMIALS[n][i]=mpz_class(cast,10);
  }
}

//returns the p^x*(1-p)^(n-x) part of the binomial distribution
//returns 0 if the number of sigmas away from the mean is larger than cut
mpf_class BinomialDistpPower(const size_t n, const double p, const size_t x, const int cut=100){
  if(p==0){
    if(x==0) return 1;
    else return 0;
  }
  else if(p==1){
    if(x==n) return 1;
    else return 0;
  }
  else if(x>n||x<0) return 0;
  else if(cut>0 && (n*p-x)*(n*p-x)>cut*cut*n*(1-p)){
    return 0;
  }

  return intpow(p,x)*intpow(1-p,n-x);

}

//takes in n (trials),p (success probability) and x (wanted sucess) and calculates the probability
//assumes if the x value if more than cut sigmas away it is just 0
long double BinomialDistribution(const size_t n, const double p, const size_t x,int cut=100){
    if(p==0){
	if(x==0) return 1;
	else return 0;
    }
    else if(p==1){
	if(x==n) return 1;
	else return 0;
    }
    else if(x>n||x<0) return 0;
    else if(cut>0 && (n*p-x)*(n*p-x)>cut*cut*n*(1-p)){
	return 0;
    }

    mpf_class tmp=Binomial(n,x)*intpow(p,x)*intpow(1-p,n-x);

    return tmp.get_d();
}
*/

//takes in input dat and outputs data for pmin to pmax using Np points
//each processor will calculate only every w.size() p-values starting at w.rank()
std::vector<std::vector<long double> > reaverage(
						 const std::vector<std::vector<mpf_class> > &dat,
						 const long double pmin,
						 const long double pmax,
						 const long int Np,
						 const boost::mpi::communicator& w,
						 const int scut=-1,
						 const bool silentmode=true){
  if(pmin<0||pmin>1||pmax<0||pmax>1){
    std::cout << "pmin or pmax are out of range:\n"
	      << "pmin = " << pmin << "\n"
	      << "pmax = " << pmax << std::endl;
    throw;
  }
  
  size_t Nb=dat.size()-1,Nc=dat[0].size();//bond count, column count
  long double dp=(pmax-pmin)/(Np-1);

  for(size_t i=1;i<=Nb;i++){
    if(dat[i].size()!=Nc){
      std::cout << "All rows do not have the same column count\n"
		<< "row " << i << " has " << dat[i].size() << "columns,\n"
		<< "but row 0 has " << Nc << " columns" << std::endl;
      throw;
    }
  }

  //each processor only calculates every w.rank() p-value
  std::vector< std::vector<long double> > pdat(Np);
  //initializes pdat
  for(long int i=w.rank();i<Np;i+=w.size()){
    pdat[i].resize(Nc,0);
    pdat[i][0]=pmin+dp*i;
  }

  //InitializeBinomialList(Nb,w);
  //if(!silentmode) std::cout << w.rank() << ": Initialized Binomial List" << std::endl;

  long double bn;
  long long int xL=0;
  size_t xU=Nb;
  boost::math::binomial dist;
  //loops over each p-values
  for(long int i=w.rank();i<Np;i+=w.size()){
    dist=boost::math::binomial(Nb,pdat[i][0]);
    //loops over each row in dat (bond counts)
    if(scut > 0){
      xL=floor(mean(dist)-scut*sqrt(variance(dist)));
      xU=ceil(mean(dist)+scut*sqrt(variance(dist)));
      if(xL<0) xL=0;
      if(xU>Nb) xU=Nb;
    }

    //std::cout << "(" << pdat[i][0] << "," << xL << "," << xU << "), ";

    for(size_t x=xL;x<=xU;x++){
      //loops over each column in dat[i]
      bn=boost::math::pdf(dist,x);
      for(size_t c=1;c<Nc;c++){
	pdat[i][c]+=bn*dat[x][c].get_d();
      }
    }
  }
  //std::cout << std::endl;

  /*
  //reuses Bnx=Binomial Coeff until finished with it
  mpz_class Bnx;

  //true until Bnxp is non-zero then turns false
  bool lowerQ;

  if(!silentmode) std::cout << w.rank() << ": averaging over data" << std::endl;

  //binomial dist value
  mpf_class Bnxp;
  //loops over each x
  for(size_t x=0;x<=Nb;x++){
    //if(!silentmode) std::cout << w.rank() << ": x=" << x << std::endl;
    
    try{
      Bnx=BinomialDirect(Nb,x);
    }
    catch(...){
      std::cerr << w.rank() << ": Failed to calculate Bnx\n"
		<< "n=" << Nb << "\n"
		<< "x=" << x << std::endl;
      throw;
    }
    
    lowerQ=true;
    //loops over each probability which is stored in pdat
    for(long int i=w.rank();i<Np;i+=w.size()){
      //if(!silentmode) std::cout << w.rank() << ": i=" << i << std::endl;
      
      Bnxp=Bnx*BinomialDistpPower(Nb,pdat[i][0],x,scut);
      if(Bnxp==0){
	if(lowerQ) continue; //if still under pL continues up
	else break; //otherwise their will be no more non-zeros
      }
      else if(lowerQ) lowerQ=false;
      //adds up all the dat values
      for(size_t c=1;c<Nc;c++){
	pdat[i][c]+=Bnxp.get_d()*dat[x][c].get_d();
      } 
    }
  }
  */

  return pdat;
}

//be is the sum of all the bond counts where a a homology first appears
void OutputData(std::vector< std::vector<mpf_class> > &sum,
		const int Ns,
		const int Nb,
		const int NA,
		const int n,
		const std::string datfile,
		const long double pmin,	
		const long double pmax,	
		const int Np,
		const boost::mpi::communicator& w,
		const int scut=-1,
		const bool silentmode=true
	){
    //outputs graph to file---------------------------------------------------  
    std::ofstream dat;
    long double local=0,global=0;

    //Sums and share data across processors
    for(size_t i=0; i<sum.size(); i++){
      for(size_t k=1; k<sum[i].size(); k++){
	local=sum[i][k].get_d()/NA;
	//MPI_Allreduce(&local,&global,1,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	all_reduce(w,local,global, std::plus<long double>() );
	sum[i][k]=mpf_class(std::to_string(global));
      }
    }

    if(!silentmode) std::cout << w.rank() << ": Finished sharing data" << std::endl;

    //reaverages the data
    std::vector<std::vector<long double> > pdat=reaverage(sum,pmin,pmax,Np,w,scut,silentmode);
    
    if(!silentmode) std::cout << w.rank() << ": created pdat" << std::endl;

    //opens and clears file and outputs header information
    if(w.rank()==0){
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
  
    w.barrier();
    dat.open(datfile.c_str(),std::ofstream::app);

    if(dat.is_open()){
      //Outputting data to file
      for(long int i=0;i<Np;i++){
	if(w.rank()==i%w.size()){
	  for(size_t j=0;j<pdat[i].size(); j++){
	    dat << std::scientific << std::setprecision(8) << pdat[i][j] << "\t";
	  }
	  dat << std::endl;
	}
	w.barrier();
      }
    }
    else{
      std::cout << "id " << w.rank() << " could not openfile " << datfile << std::endl;
      throw;
    }
    dat.close();
    return;
}

#endif
