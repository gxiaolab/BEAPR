#ifndef TEXTIO_H
#define TEXTIO_H

#include<string>
#include<map>
#include<set>
#include<vector>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include <sstream>
#include <list>
#include <assert.h>   
#include "TOK.h"

class TextIO{
	public :
	TextIO();
	std::vector< std::vector<std::string> > dlread(std::string fname,  std::string delimiters );
	std::vector< std::vector<double> > dlread_double(std::string fname,  std::string delimiters );
	std::vector< std::vector<double> > dlread_sparse(std::string fname);
	std::vector<std::string> dlread_vector(std::string fname);

	void dlwrite(std::string fname, std::vector< std::vector<std::string> > tokvec , std::string delimiters );
	void dlwrite_double(std::string fname, std::vector< std::vector<double> > tokvec , std::string delimiters );
	void dlwrite(std::string fname,  std::vector<std::string> tokvec );	

	void dlwrite_sparse(std::string fname, std::vector< std::vector<double> > tokvec )	;
	void dlwrite_withHead(std::string fname, std::vector<std::string> HeadVec, std::vector< std::vector<std::string> > tokvec , std::string delimiters );	
	void dlwrite_withHead(std::string fname, std::vector<std::string> HeadVec, std::vector< std::vector<double> > tokvec , std::string delimiters );	
};


#endif