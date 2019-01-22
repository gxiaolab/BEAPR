#ifndef SIM_H
#define SIM_H

#include<string>
#include<map>
#include<set>
#include<vector>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include <sstream>
#include <list>
#include <climits>
#include <algorithm>
#include <math.h>
#include <stdio.h>      
#include <stdlib.h>     
#include <time.h> 
#include "TOK.h"
#include "TextIO.h"

class Sim{
	public :
	
	void createTestMean(std::string outfname, std::map<int, double> cv2map,std::map<int, double> distmap, double libsize, double allele_bias);
	void createOvTestMean(std::string outfname, std::map<int, double> cv2map,std::map<int, double> distmap, double libsize, double allele_bias);
	void createReg(std::string regfname, std::string predfname, std::string outfname);	
  std::map<int,double> learnDist(std::string predfname);
		 std::map<int,double> learnDist_new(std::string predfname);
	std::map<int,double> learnCV2(std::string regfname);
		
		
	private:
	TextIO tio;
	TOK tok;
			
};


#endif