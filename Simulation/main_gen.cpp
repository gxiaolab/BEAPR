#include "Sim.h"

#include <string>
#include <stdlib.h>
#include <stdio.h>


using namespace std;

int main(int argc, char **argv) {
	
	Sim sm ; 
	stringstream ss; 
	double lib_size = atof(argv[1]);
	double allele_bias = atof(argv[2]); ;
	
	map<int, double> cv2map =  sm.learnCV2( "/home/yyang/data/eclip/K562/SRSF1_test/tmp/reg.rc");
    //map<int, double> distmap =  sm.learnDist( "/home/yyang/data/eclip/K562/SRSF1_test/tmp/pred.rc");
    map<int, double> distmap =  sm.learnDist_new( "/home/yyang/projects/eCLIP/motif/chronical/08192018_SRSF1_readcoverage/dep.txt");
    
	//sm.createTestMean("/home/yyang/data/eclip/K562/SRSF1_test/sim/mean", cv2map, distmap, lib_size , allele_bias );
	ss.str("");
	ss<< "/home/yyang/data/eclip/K562/SRSF1_test/sim/"; 
	
	sm.createTestMean( ss.str() + "/mean", cv2map, distmap, lib_size , allele_bias );
	
	//sm.createOvTestMean( ss.str() + "/mean", cv2map, distmap, lib_size , allele_bias );
	
	//sm.createReg( "/home/yyang/data/eclip/K562/SRSF1_test/sim/reg.rc", "/home/yyang/data/eclip/K562/SRSF1_test/sim/mean.out" , "/home/yyang/data/eclip/K562/SRSF1_test/sim/reg_filtered.rc" ); 
	
	return 0;
	
}
