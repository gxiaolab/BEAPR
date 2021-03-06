#ifndef NORMALIZER_H
#define NORMALIZER_H

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
#include "GenomeSegment.h"
#include "TOK.h"
#include "TextIO.h"

class Normalizer{
	public :
		void normalize(std::string pkbed_fname, std::string input_fname,std::string r1_fname,std::string tmp_dir);
		std::vector<std::string> normalize_chr(std::string pkbed_fname, std::string input_fname,std::string r1_fname,  double input_size,double r1_size,double cut);	
		std::vector<std::string> normalize_chr_new(std::string pkbed_fname, std::string r1_fname, std::string input_fname, double input_size,double r1_size,double cut);	
		void filterPKs(std::string r1_peak_fname, std::string r1_mp_fname, std::string input_mp_fname , std::string out_peak_fname,double input_size,double r1_size ,double cut);
		void filterPKs(std::string r1_peak_fname, std::string out_peak_fname , double cut);	
	  //2017.7.17
	  std::map< int, std::vector<double> > estimateSeqError( std::string mpfname, std::string SNPfname,double lib_size,int pos_flag);
	  
			
	private:
		TextIO tio;
		TOK tok;
		int _mid_point(std::string sam_str);
		std::vector<std::string> _sepBaseStr(std::string baseStr);	
		double	_depth(std::string baseStr);
		double	_match(std::string baseStr);	
};


#endif