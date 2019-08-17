#include "Motif.h"

using namespace std;

std::vector< std::string>	Motif::_filter_DGE_SNP(  std::vector<std::string> PredVec ,  std::string ASEfname , std::string gene_tab ){
	
	vector<string> outVec; 
	map<string, int > gene_s_map, gene_e_map;
	map<string, string > gene_chr_map; 
	vector<string> ASEVec = tio.dlread_vector(ASEfname); 
	vector<string> TabVec = tio.dlread_vector(gene_tab);
	
	//init gene map;
	for(int i =0 ; i<ASEVec.size();  i++){
		if(ASEVec[i].length()>0){
			string tmpstr = ASEVec[i].substr(0,3);
			if(tmpstr.compare("chr") == 0){
				vector<string> tokvec = tok.Tokenize( ASEVec[i], "\t");
				gene_s_map.insert(make_pair( tokvec[1], INT_MAX ));
				gene_e_map.insert(make_pair( tokvec[1], -1 ));
				gene_chr_map.insert(make_pair( tokvec[1], "" ));
			}
		}
	}
	
	
	for(int i=1 ; i<TabVec.size() ; i++){
		vector<string> tokvec =  tok.Tokenize( TabVec[i], "\t");
		string gname =tokvec[0];
		if(gene_s_map.find(gname) != gene_s_map.end()){
			gene_chr_map[gname] = tokvec[2] +"_"+tokvec[3] ;
			int s_pos = atoi(tokvec[4].c_str());
			int e_pos = atoi(tokvec[5].c_str());
			if(s_pos<gene_s_map[gname] ){
				gene_s_map[gname] = s_pos;
			}
			if(e_pos>gene_e_map[gname]){
				gene_e_map[gname] = e_pos;
			}
		}
	}
	
	
	for(int i=0 ; i<PredVec.size() ; i++){
		vector<string> tokvec =  tok.Tokenize( PredVec[i], "\t");
		string snp_chr = tokvec[1]; 
		int snp_pos = atoi(tokvec[2].c_str());
		string snp_strand = tokvec[3];
		
		int hit_flag =0;
		for(map<string,string>::iterator m_itr = gene_chr_map.begin() ; m_itr != gene_chr_map.end() ; m_itr++){
			 
			string gname = m_itr->first;
			string chr_str = snp_chr +"_"+snp_strand;
			if(chr_str.compare(m_itr->second)==0){
				if(snp_pos >= gene_s_map[gname] && snp_pos <= gene_e_map[gname]){
					hit_flag =1;
					cout<<PredVec[i] << endl;
					
					cout<< "is in " << gname << endl;
				}
				
			}
			
		}
		
		if(hit_flag==0){
			outVec.push_back(PredVec[i]);
		}
	}
	

	return outVec; 
}


