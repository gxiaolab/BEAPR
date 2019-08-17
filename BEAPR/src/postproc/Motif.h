#ifndef MOTIF_H
#define MOTIF_H

#include<string>
#include <algorithm>
#include<map>
#include<set>
#include<vector>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include <sstream>
#include <list>
#include <climits>
#include <math.h>
#include <time.h> 
#include "TOK.h"
#include "TextIO.h"


class Motif{
	public:
		
	std::vector< std::string>	_filter_DGE_SNP(  std::vector<std::string> PredVec ,  std::string ASEfname , std::string gene_tab );
	
		
	private:
		TextIO tio;
		TOK tok;
		
};
#endif