#ifndef SNP_H
#define SNP_H

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
#include "GenomeSegment.h"

class SNP{
	public :
			
	void listKnownSNPs( std::string snp_fname, std::string pk_fname , std::string mp_fname, std::string out_fname,int pos_pks ,double min_score );
	void listAllOvSNPs( std::string snp_fname, std::string pk_fname , std::string mp_fname, std::string out_fname,int pos_pks ,double min_score );
	
	void listAllOvSNPs_VCF( std::vector<std::string> ChrSNPVec , std::string pk_fname , std::string mp_fname, std::string out_fname,int pos_pks ,double min_score );
	
	//1.3.2018 simpy ov hetero. SNPs to pks 
	void listAllOvSNPs_VCF_lite( std::vector<std::string> ChrSNPVec , std::string pk_fname , std::string out_fname);
	void listAllOvSNPs_denovo_lite( std::vector<std::string> ChrSNPVec , std::string pk_fname , std::string out_fname);

	void readRawRC(std::vector<std::string> ChrSNPVec, std::string mp_fname, std::string out_fname );
	void readRawRC_denovo(std::vector<std::string> ChrSNPVec, std::string mp_fname, std::string out_fname );
	
	void listInputSNPs_depth( std::string snp_fname,  std::string mp_fname, std::string out_fname, int min_score  );
	
	
	std::vector<int> listCLS_pos(  std::vector<int> SNPPosVec,  std::vector<int> ClsPosVec ,int pos_flag);

	std::vector<std::string>  normalize( std::string snp_fname,  std::string cls_fname	,  std::string prior_dir  , int pos_flag,double lib_size);		
	std::vector<std::string>  normalize_VCF( std::string snp_fname,  std::string cls_fname	,  std::string prior_dir  , int pos_flag , double lib_size);		
	std::vector<std::string>  normalize_denovo( std::string snp_fname,  std::string cls_fname	,  std::string prior_dir  , int pos_flag , double lib_size);		

	std::vector<std::string> _sepBaseStr(std::string baseStr);
		
	std::vector<double> getNullScoreDist(std::string snp_fname, std::string tmp_dir);
	void SNPtoRC(	std::string SNP_dir, 	std::string tmp_dir, std::string Rep_prefix , std::string chr_str, double lib_size);
		
	void SNPtoRC_denovo(	std::string SNP_dir, 	std::string tmp_dir, std::vector< std::string > RepIDVec , std::string chr_str, std::map<std::string, double> SizeMap);	
	void SNPtoRC_train(	std::string SNP_dir, 	std::string tmp_dir, std::vector< std::string > RepIDVec , std::string chr_str, std::map<std::string, double> SizeMap);	
	void SNPtoRC_VCF(	std::string VCFfname, 	std::string tmp_dir, std::vector< std::string> RepIDVec , std::string chr_str, std::map<std::string, double> SizeMap);	
	void SNPtoRC_test(	std::string VCFfname, 	std::string tmp_dir, std::vector< std::string> RepIDVec , std::string chr_str, std::map<std::string, double> SizeMap);
	private:
	TextIO tio;
	TOK tok;
			
};


#endif