#ifndef GENOMESEGMENT_H
#define GENOMESEGMENT_H

#include<set>
#include<list>
#include<vector>
#include<iostream>
#include<algorithm>
#include<string>
#include <sstream>
#include <float.h>
#include <math.h>
#include <assert.h>   
#include <fstream>
#include <ctype.h>

#include "TOK.h"
#include "TextIO.h"

class GenomeSegment{
	public:
	
	std::vector<int> PosVec;
	std::map<std::string, std::map<std::string, std::string> > ExonMap;
	std::map<std::string,std::map<std::string, std::string> >	ReadMap;
	std::map<std::string, double > ExonCountMap;
		
	std::vector<double> readscov(std::string readsfname, int read_len);	
	std::string find_interval(int exon_pos) ;

	std::vector<std::string> get_intervals_by_pos(int s_pos, int e_pos) ;
	void	initGenomeSegment(std::string fname,int chr_len);
	void	initGenomeSegment_pksbed(std::string fname,int chr_len);
	void	initGenomeSegment(std::string fname,int chr_len,std::string gname_chk);
	void map_reads_bed(std::string gtfname,std::string bedfname);
	std::vector<std::string> segCorr();		
	void _debug();	
	void printReads(std::string dir_str);
	
	void sep_gtfs_by_genes(std::string fnamestr,std::string out_dir);	
		
	std::map<std::string, std::string> parse_gid(std::string str);
	
		 
};

#endif