#ifndef TOK_H
#define TOK_H

#include<string>
#include<map>
#include<set>
#include<vector>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include <sstream>
#include <list>

class TOK{
	public :
	TOK();
	std::vector<std::string> Tokenize( std::string str,  std::string delimiters );
		
};


#endif