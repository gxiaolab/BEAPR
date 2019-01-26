#ifndef REG_H
#define REG_H

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
#include <dirent.h>
#include "TOK.h"
#include "TextIO.h"

class Reg{
	public :
	std::vector<std::string > mergeRCs( std::string Rep_str ,std::string tmp_dir, std::string RBP_id , double min_rc);
		
	std::vector<std::string > calcVar(std::string rcfname );
		
		
	void calcFDR(std::string dir_str) ; 
			
	private:
	TextIO tio;
	TOK tok;
			
};


#endif