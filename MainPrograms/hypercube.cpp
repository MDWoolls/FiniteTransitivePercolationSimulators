#include <itpp/itcomm.h>

#include <list>
#include <iostream>
#include <iomanip>
#include <cstdlib>
//#include <climits>
#include <stdio.h>

//#include <assert.h>

//max L = 1518500249
//max long int = 18446744073709551615

int main(int argc, char **argv){
    //Intializes things from inputs----------------------------------------------------
    itpp::Parser p; 
    p.init(argc,argv);
    if(argc>1){
	if((strcmp(argv[1],"--help")==0)||(strcmp(argv[1],"-h")==0)){
	    std::cout << argv[0] << ":Outputs a hypercube with periodic boundry conditions as an incidence matrix:\n"
		      << "By default:\n"
		      << "   file=default.mtx - file hypercube is saved in\n"
		      << "   L=10 - size of each size of the hypercube\n"
		      << "   D=2 - dimension of the hypercube\n"
                      << std::endl;
      	    exit(-1);
	}
    }
    p.set_silentmode(true);

    std::string file="default.mtx";  
    p.get(file,"file");

    
    size_t L=10,D=2;
    p.get(L,"L");
    p.get(D,"D");

    //list of powers of L
    std::vector< unsigned long long int > Lpow={1};
    for(size_t i=1;i<=D+1;i++)
      {
	Lpow.push_back(L*Lpow.back());
      }

    /*
    std::cout << "Powers of L = ";
    for(const auto& pow : Lpow) std::cout << pow << " ";
    std::cout << std::endl;
    */

    //number of sites
    //current bond location
    unsigned long long int N=Lpow[D],bond=1;
    //std::cout << "N=" << N << std::endl;

    //site the current site has a bond to connect too
    unsigned long long int s;

    std::ofstream dat;
    dat.open(file.c_str());
    dat << "%%MatrixMarket matrix coordinate integer general\n"
        << "%%Hypercube with L=" << L << " and D=" << D << "\n"
        << N << " " << D*N << " " << 2*D*N << std::endl;
    
    //loops over every site
    for(size_t i=0; i<N; i++)
      {
	//std::cout << "i=" << i << std::endl;
	//loops over each bond direction
	for(size_t j=0; j<D; j++)
	  {
	    //std::cout << "j=" << j << std::endl;
	    //calculates site position to connect to site i
	    s=(i+Lpow[j])%Lpow[j+1];
	    for(size_t k=j;k<D;k++) s+=(Lpow[k+1]*(i/Lpow[k+1]))%Lpow[k+2];
	    //std::cout << "s=" << std::endl;
	    
	    dat << i+1 << " " << bond << " " << 1 << "\n"
		<< s+1 << " " << bond << " " << 1 << "\n";
	    bond++;
	  }
	dat.flush();
      }
    dat.close();
    return 0;
}
