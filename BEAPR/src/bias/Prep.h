#ifndef PREP_H
#define PREP_H

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
#include "TOK.h"
#include "TextIO.h"

class Prep{
	public :
	void sep_mp_byChrms(std::string  mp_fname , std::string tmpdir , std::string prefix );	
	void sep_reads_by_chrs(std::string fnamestr, std::string out_dir , std::string prefix);
	void preProc(std::string inBamfname ,std::string Ref_fname ,std::string tmpdir ,std::string prefix); 	
  void sep_beds_by_chrs(std::string inBedfname, std::string tmpdir,std::string prefix);
  //for Gio's pipeline
  void preProc_gio(std::string inBamfname ,std::string Ref_fname ,std::string tmpdir ,std::string prefix);	 		
	
	private:
	TextIO tio;
	TOK tok;		
};


#endif