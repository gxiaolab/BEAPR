#ifndef CROSSLINKER_H
#define CROSSLINKER_H

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

class CrossLinker{
	public :
		void callCLs(std::string pk_fname,  std::string mp_fname, std::string out_fname, double fdr_cut, char delimiter);
		std::vector<std::vector<double> >  calcPrior(std::vector< std::string>  CLfFnameVec, std::vector< std::string> SEQfnameVec, std::string out_dir, int half_size, int pos_flag);
		void chkMut(std::string mp_fname, std::string mate_fname);	
		std::map<int, double> null_pdf(int len, int site_num, int itr_num);	
		void predCLs(std::string bed_fname, std::string mp_fname , std::string out_fname, std::string strand_str , std::string FDR_str  );
		void prior(std::string cls_dir,std::string seq_dir , std::string out_dir, std::string chr_str,std::string size_str);	
					
	private:
		TextIO tio;
		TOK tok;
		std::vector<int> _mpstr_calc(std::string mpstr);
		
};


#endif