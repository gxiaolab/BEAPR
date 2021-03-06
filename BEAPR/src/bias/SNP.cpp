#include "SNP.h"

using namespace std;

std::vector<std::string> SNP::_sepBaseStr(std::string baseStr){
	vector<string> outVec;
	stringstream ss;
	for(int i =0; i<baseStr.length(); ){
		int k=1;
		
		if(baseStr[i]=='.' || baseStr[i]==',' || baseStr[i]=='>' || baseStr[i]=='<'|| baseStr[i]=='A' || baseStr[i]=='a'|| baseStr[i]=='T' || baseStr[i]=='t'|| baseStr[i]=='C' || baseStr[i]=='c'|| baseStr[i]=='G' || baseStr[i]=='g' || baseStr[i]=='*'){
			ss.str("");
			ss<< baseStr[i];
			outVec.push_back(ss.str());
		}
		if(baseStr[i]=='^'){
			ss.str("");
			ss<< baseStr[i+2];
			outVec.push_back(ss.str());
			k=3;
		}
		if(baseStr[i]=='+' || baseStr[i]=='-'){
			ss.str("");
			int d_num=0;
			int str_len ;
			ss << baseStr[i];
			
			while(baseStr[i+1+d_num]>=48 && baseStr[i+1+d_num]<=57 ){
				ss << baseStr[i+1+d_num];
				d_num++;
			}
			
			string tmp_str = ss.str();
			
			tmp_str = tmp_str.substr(1);
			str_len = atoi(tmp_str.c_str()); 
			
			for(int j =1 ; j<=str_len; j++){
				ss<< baseStr[i+d_num+j];
			}
	
			k = d_num+ str_len+1;
		}
		
		i= i+k;
	}
	
	return outVec;
}
//how 

std::vector<double> SNP::getNullScoreDist(std::string snp_fname, std::string mp_fname){
	map<int,string> SNPsMap;
	vector<double> outVec;
	vector<string> ChrSNPVec = tio.dlread_vector(snp_fname);
	for(int i=0 ; i<=70; i++){
		outVec.push_back(0);
	}
	
	for(int i =0 ; i<ChrSNPVec.size(); i++){
		vector<string> tokvec = tok.Tokenize(ChrSNPVec[i], "\t");
		int snp_pos = atoi(tokvec[1].c_str());
		snp_pos = snp_pos +1;
		
		SNPsMap.insert(make_pair(snp_pos,""));
		
	}
	
	vector<string> mpVec = tio.dlread_vector(mp_fname );
	for(int i =0 ; i< mpVec.size() ; i++){
		vector<string> tokvec = tok.Tokenize(mpVec[i],"\t");
		int pos = atoi(tokvec[1].c_str()) ;
		int cov_num = atoi(tokvec[3].c_str());
		if(cov_num > 0 ){ 
			string baseStr = tokvec[4];
			string scoreStr =  tokvec[5];
			//if not a SNP site
			if(SNPsMap.find(pos)==SNPsMap.end()){
				vector<string> baseVec =  _sepBaseStr(baseStr);
				if(baseVec.size()== scoreStr.length() ){
					
					for(int j =0 ; j<scoreStr.length(); j++){
						if(baseVec[j].compare("a")==0 || baseVec[j].compare("A")==0 || baseVec[j].compare("t")==0 || baseVec[j].compare("T")==0 ||baseVec[j].compare("c")==0 || baseVec[j].compare("C")==0 || baseVec[j].compare("g")==0 || baseVec[j].compare("G")==0 ){
							//cout <<baseStr << ":"<< scoreStr << endl;
							
							int score_idx = (int)(scoreStr[j]);
							outVec[score_idx] = outVec[score_idx]+1;
							
						}
					}				
				}
			}
			
		}



	}
	
	return outVec;
}


//find SNPs in PKs
void SNP::listKnownSNPs( std::string snp_fname, std::string pk_fname , std::string mp_fname, std::string out_fname, int pos_pks ,double min_score  ){
	
	stringstream ss;
	GenomeSegment gs ; 
	gs.initGenomeSegment_pksbed(  pk_fname ,INT_MAX	 );

	vector<string> ChrSNPVec = tio.dlread_vector(snp_fname);
	vector<int> PosVec;
	vector<string> ASVec;
	map<int,string> SNPMap ; 
	map<int,string> MPMap ;
	vector<string> outVec;
	for(int i =0 ; i<ChrSNPVec.size(); i++){
		
		vector<string> tokvec = tok.Tokenize(ChrSNPVec[i], "\t");
		int snp_pos = atoi(tokvec[1].c_str());
		snp_pos ++;
		
		string idx_str = gs.find_interval(snp_pos);
		vector<string> tmpvec = tok.Tokenize(idx_str, "_");
		int a_idx = atoi(tmpvec[0].c_str());
		int b_idx = atoi(tmpvec[1].c_str());
		ss.str("");
		ss<< gs.PosVec[a_idx]<< "_"<< gs.PosVec[b_idx];
		
		if( gs.ExonMap[ss.str()].size() >0 && tokvec[6].compare("single")==0){
			int SNP_strand =1;
			if(tokvec[4].compare("-")==0){
				SNP_strand =-1;
			} 
			string A_str = tokvec[7];
			if(SNP_strand * pos_pks == -1){
				//reverse the alleles
				ss.str("");
				
				for(int j =0 ; j< A_str.length(); j++){
					
					if(A_str[j] == 'A' || A_str[j] == 'a' ){
						ss<< "T";
					}
					if(A_str[j] == 'T' || A_str[j] == 't' ){
						ss<< "A";
					}
					if(A_str[j] == 'C' || A_str[j] == 'c' ){
						ss<< "G";
					}
					if(A_str[j] == 'G' || A_str[j] == 'g' ){
						ss<< "C";
					}
					if(A_str[j] == ',' ){
						ss<< ",";
					}
				}
				A_str = ss.str();
				
			}
			
		//	cout <<ChrSNPVec[i] << "\t" << A_str<< endl;
			
			SNPMap.insert(make_pair( snp_pos ,A_str ));
		} 
	}	
	
	vector<string> MPVec = tio.dlread_vector(mp_fname);
	for(int i =0 ; i<MPVec.size(); i++){
			vector<string> tokvec = tok.Tokenize(MPVec[i],"\t");
			int mp_pos = atoi(tokvec[1].c_str()); 
			if(atoi(tokvec[3].c_str())>0){
				MPMap.insert( make_pair(mp_pos, MPVec[i]));
			}
	}
	
	for(map<int,string>::iterator m_itr = SNPMap.begin() ; m_itr != SNPMap.end(); m_itr++){
		//cout << "look at "<< m_itr->first << "\t" << m_itr->second<< endl;
		int snp_pos = m_itr->first;
		vector<string> tokvec = tok.Tokenize( m_itr->second , ",");
		
		map<string, double> chkSNPMap ; 
		
		for(int i =0 ; i< tokvec.size(); i++){
			if(tokvec[i].length()>0){
				chkSNPMap.insert(make_pair(tokvec[i] , 0.0));
			}
		}
		 
		double match_num =0;
		
		if(MPMap.find(snp_pos) != MPMap.end()){
			//cout << MPMap[snp_pos] << endl;
		 	vector<string> tmpvec  = tok.Tokenize( MPMap[snp_pos] , "\t") ; 
		 	string  base_str = tmpvec[4];
		 	string  score_str = tmpvec[5];
		 	vector<string> baseVec = _sepBaseStr(base_str);
		 	for( int j =0 ; j<baseVec.size() ; j++){
		 		int score =(int) score_str[j]; 
		 		if(score> min_score){
		 			if(baseVec[j].compare(".")==0 || baseVec[j].compare(",")==0 ){
		 				match_num = match_num +1.0;
		 			}
		 			if(baseVec[j].compare("A")==0 || baseVec[j].compare("a")==0 ){
		 				if(chkSNPMap.find("A") != chkSNPMap.end()){
		 					chkSNPMap["A"] = chkSNPMap["A"] + 1.0;
		 				}
		 			}
		 			
		 			if(baseVec[j].compare("T")==0 || baseVec[j].compare("t")==0 ){
		 				if(chkSNPMap.find("T") != chkSNPMap.end()){
		 					chkSNPMap["T"] = chkSNPMap["T"] + 1.0;
		 				}
		 			}
		 			
		 			if(baseVec[j].compare("C")==0 || baseVec[j].compare("c")==0 ){
		 				if(chkSNPMap.find("C") != chkSNPMap.end()){
		 					chkSNPMap["C"] = chkSNPMap["C"] + 1.0;
		 				}
		 			}
		 			if(baseVec[j].compare("G")==0 || baseVec[j].compare("g")==0 ){
		 				if(chkSNPMap.find("G") != chkSNPMap.end()){
		 					chkSNPMap["G"] = chkSNPMap["G"] + 1.0;
		 				}
		 			}
		 				
		 		}
		 		
		 	}
		 	ss.str("");
		 	ss<< snp_pos ;
		 	ss<< "\t"<< match_num ;
		 	ss<< "\t";
		 	double mut_num =0;
		 	for(map<string,double>::iterator m_itr = chkSNPMap.begin(); m_itr != chkSNPMap.end() ; m_itr++){
		 		
		 		ss<<  m_itr->first<< ":" << m_itr->second <<"," ;
		 		mut_num = mut_num+m_itr->second;
		 	}
		 	if(mut_num>0 && match_num>0 ){
		 		outVec.push_back(ss.str());
			}
		}
		
	}
	
	tio.dlwrite(out_fname, outVec);
	
}

//
//void SNP::listAllOvSNPs( std::string snp_fname, std::string pk_fname , std::string mp_fname, std::string out_fname, int pos_pks ,double min_score  ){
//	
//	stringstream ss;
//	GenomeSegment gs ; 
//	gs.initGenomeSegment_pksbed(  pk_fname ,INT_MAX	 );
//
//	vector<string> ChrSNPVec = tio.dlread_vector(snp_fname);
//	vector<int> PosVec;
//	vector<string> ASVec;
//	map<int,string> SNPMap ; 
//	map<int,string> MPMap ;
//	vector<string> outVec;
//	for(int i =0 ; i<ChrSNPVec.size(); i++){
//		
//		vector<string> tokvec = tok.Tokenize(ChrSNPVec[i], "\t");
//		int snp_pos = atoi(tokvec[1].c_str());
//		snp_pos ++;
//		
//		string idx_str = gs.find_interval(snp_pos);
//		vector<string> tmpvec = tok.Tokenize(idx_str, "_");
//		int a_idx = atoi(tmpvec[0].c_str());
//		int b_idx = atoi(tmpvec[1].c_str());
//		ss.str("");
//		ss<< gs.PosVec[a_idx]<< "_"<< gs.PosVec[b_idx];
//		if( gs.ExonMap[ss.str()].size() >0 && tokvec[6].compare("single")==0){
//			int SNP_strand =1;
//			if(tokvec[4].compare("-")==0){
//				SNP_strand =-1;
//			} 
//			string A_str = tokvec[7];
//			if(SNP_strand * pos_pks == -1){
//				//reverse the alleles
//				ss.str("");
//				
//				for(int j =0 ; j< A_str.length(); j++){
//					
//					if(A_str[j] == 'A' || A_str[j] == 'a' ){
//						ss<< "T";
//					}
//					if(A_str[j] == 'T' || A_str[j] == 't' ){
//						ss<< "A";
//					}
//					if(A_str[j] == 'C' || A_str[j] == 'c' ){
//						ss<< "G";
//					}
//					if(A_str[j] == 'G' || A_str[j] == 'g' ){
//						ss<< "C";
//					}
//					if(A_str[j] == ',' ){
//						ss<< ",";
//					}
//				}
//				A_str = ss.str();
//				
//			}
//			
//		//	cout <<ChrSNPVec[i] << "\t" << A_str<< endl;
//			
//			SNPMap.insert(make_pair( snp_pos ,A_str ));
//		} 
//	}	
//	
//	vector<string> MPVec = tio.dlread_vector(mp_fname);
//	for(int i =0 ; i<MPVec.size(); i++){
//			vector<string> tokvec = tok.Tokenize(MPVec[i],"\t");
//			int mp_pos = atoi(tokvec[1].c_str()); 
//			if(atoi(tokvec[3].c_str())>0){
//				MPMap.insert( make_pair(mp_pos, MPVec[i]));
//			}
//	}
//	
//	for(map<int,string>::iterator m_itr = SNPMap.begin() ; m_itr != SNPMap.end(); m_itr++){
//		//cout << "look at "<< m_itr->first << "\t" << m_itr->second<< endl;
//		int snp_pos = m_itr->first;
//		
//		map<string, double> chkSNPMap ; 
//		
//		chkSNPMap.insert(make_pair("A" , 0.0));
//		chkSNPMap.insert(make_pair("T" , 0.0));
//		chkSNPMap.insert(make_pair("C" , 0.0));
//		chkSNPMap.insert(make_pair("G" , 0.0));
//		
//		double match_num =0;
//		
//		if(MPMap.find(snp_pos) != MPMap.end()){
//			//cout << MPMap[snp_pos] << endl;
//		 	vector<string> tmpvec  = tok.Tokenize( MPMap[snp_pos] , "\t") ; 
//		 	string ref_str = tmpvec[2]; 
//		 	string  base_str = tmpvec[4];
//		 	string  score_str = tmpvec[5];
//		 	vector<string> baseVec = _sepBaseStr(base_str);
//		 	
//		 	if(ref_str.compare("a")==0){
//		 		ref_str="A";
//		 	}
//		 	if(ref_str.compare("t")==0){
//		 		ref_str="T";
//		 	}
//		 	if(ref_str.compare("c")==0){
//		 		ref_str="C";
//		 	}
//		 	if(ref_str.compare("g")==0){
//		 		ref_str="G";
//		 	}
//		 	
//		 	for( int j =0 ; j<baseVec.size() ; j++){
//		 		int score =(int) score_str[j]; 
//		 		if(score> min_score){
//		 			if(baseVec[j].compare(".")==0 || baseVec[j].compare(",")==0 ){
//		 				match_num = match_num +1.0;
//		 				chkSNPMap[ref_str] = chkSNPMap[ref_str] + 1.0;
//		 			}
//		 			if(baseVec[j].compare("A")==0 || baseVec[j].compare("a")==0 ){
//		 				if(chkSNPMap.find("A") != chkSNPMap.end()){
//		 					chkSNPMap["A"] = chkSNPMap["A"] + 1.0;
//		 				}
//		 			}
//		 			
//		 			if(baseVec[j].compare("T")==0 || baseVec[j].compare("t")==0 ){
//		 				if(chkSNPMap.find("T") != chkSNPMap.end()){
//		 					chkSNPMap["T"] = chkSNPMap["T"] + 1.0;
//		 				}
//		 			}
//		 			
//		 			if(baseVec[j].compare("C")==0 || baseVec[j].compare("c")==0 ){
//		 				if(chkSNPMap.find("C") != chkSNPMap.end()){
//		 					chkSNPMap["C"] = chkSNPMap["C"] + 1.0;
//		 				}
//		 			}
//		 			if(baseVec[j].compare("G")==0 || baseVec[j].compare("g")==0 ){
//		 				if(chkSNPMap.find("G") != chkSNPMap.end()){
//		 					chkSNPMap["G"] = chkSNPMap["G"] + 1.0;
//		 				}
//		 			}
//		 				
//		 		}
//		 		
//		 	}
//		 	ss.str("");
//		 	ss<< snp_pos ;
//		 	ss<< "\t"<< match_num ;
//		 	ss<< "\t";
//		 	double mut_num =0;
//		 	for(map<string,double>::iterator m_itr = chkSNPMap.begin(); m_itr != chkSNPMap.end() ; m_itr++){
//		 		
//		 		ss<<  m_itr->first<< ":" << m_itr->second <<"," ;
//		 		mut_num = mut_num+m_itr->second;
//		 	}
//		 	
//		 	outVec.push_back(ss.str());
//			
//		}
//		
//	}
//	
//	tio.dlwrite(out_fname, outVec);
//	
//}

void SNP::listAllOvSNPs( std::string snp_fname, std::string pk_fname , std::string mp_fname, std::string out_fname, int pos_pks ,double min_score  ){
	
	stringstream ss;
	GenomeSegment gs ; 
	gs.initGenomeSegment_pksbed(  pk_fname ,INT_MAX	 );
	string str;
	vector<string> ChrSNPVec = tio.dlread_vector(snp_fname);
	vector<int> PosVec;
	vector<string> ASVec;
	map<int,string> SNPMap ; 
	map<int,string> MPMap ;
	vector<string> outVec;
	for(int i =0 ; i<ChrSNPVec.size(); i++){
		
		vector<string> tokvec = tok.Tokenize(ChrSNPVec[i], " ");
		int snp_pos = atoi(tokvec[1].c_str());
		//snp_pos ++;
		
		string idx_str = gs.find_interval(snp_pos);
		vector<string> tmpvec = tok.Tokenize(idx_str, "_");
		int a_idx = atoi(tmpvec[0].c_str());
		int b_idx = atoi(tmpvec[1].c_str());
		ss.str("");
		ss<< gs.PosVec[a_idx]<< "_"<< gs.PosVec[b_idx];
		if( gs.ExonMap[ss.str()].size() >0 ){
			SNPMap.insert(make_pair( snp_pos ,"" ));
		} 
	}	
	
	fstream filestr;
	filestr.open(mp_fname.c_str(), ios_base::in);
	
	while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 	
			vector<string> tokvec = tok.Tokenize( str,"\t");					
			int mp_pos = atoi(tokvec[1].c_str()); 
			if(atoi(tokvec[3].c_str())>0){
				MPMap.insert( make_pair(mp_pos, str));
			}
		}
	}
	filestr.close();
	
	for(map<int,string>::iterator m_itr = SNPMap.begin() ; m_itr != SNPMap.end(); m_itr++){
		//cout << "look at "<< m_itr->first << "\t" << m_itr->second<< endl;
		int snp_pos = m_itr->first;
		
		map<string, double> chkSNPMap ; 
		
		chkSNPMap.insert(make_pair("A" , 0.0));
		chkSNPMap.insert(make_pair("T" , 0.0));
		chkSNPMap.insert(make_pair("C" , 0.0));
		chkSNPMap.insert(make_pair("G" , 0.0));
		
		double match_num =0;
		
		if(MPMap.find(snp_pos) != MPMap.end()){
			//cout << MPMap[snp_pos] << endl;
		 	vector<string> tmpvec  = tok.Tokenize( MPMap[snp_pos] , "\t") ; 
		 	string ref_str = tmpvec[2]; 
		 	string  base_str = tmpvec[4];
		 	string  score_str = tmpvec[5];
		 	vector<string> baseVec = _sepBaseStr(base_str);
		 	
		 	if(ref_str.compare("a")==0){
		 		ref_str="A";
		 	}
		 	if(ref_str.compare("t")==0){
		 		ref_str="T";
		 	}
		 	if(ref_str.compare("c")==0){
		 		ref_str="C";
		 	}
		 	if(ref_str.compare("g")==0){
		 		ref_str="G";
		 	}
		 	
		 	for( int j =0 ; j<baseVec.size() ; j++){
		 		int score =(int) score_str[j]; 
		 		if(score> min_score){
		 			if(baseVec[j].compare(".")==0 || baseVec[j].compare(",")==0 ){
		 				match_num = match_num +1.0;
		 				chkSNPMap[ref_str] = chkSNPMap[ref_str] + 1.0;
		 			}
		 			if(baseVec[j].compare("A")==0 || baseVec[j].compare("a")==0 ){
		 				if(chkSNPMap.find("A") != chkSNPMap.end()){
		 					chkSNPMap["A"] = chkSNPMap["A"] + 1.0;
		 				}
		 			}
		 			
		 			if(baseVec[j].compare("T")==0 || baseVec[j].compare("t")==0 ){
		 				if(chkSNPMap.find("T") != chkSNPMap.end()){
		 					chkSNPMap["T"] = chkSNPMap["T"] + 1.0;
		 				}
		 			}
		 			
		 			if(baseVec[j].compare("C")==0 || baseVec[j].compare("c")==0 ){
		 				if(chkSNPMap.find("C") != chkSNPMap.end()){
		 					chkSNPMap["C"] = chkSNPMap["C"] + 1.0;
		 				}
		 			}
		 			if(baseVec[j].compare("G")==0 || baseVec[j].compare("g")==0 ){
		 				if(chkSNPMap.find("G") != chkSNPMap.end()){
		 					chkSNPMap["G"] = chkSNPMap["G"] + 1.0;
		 				}
		 			}
		 				
		 		}
		 		
		 	}
		 	ss.str("");
		 	ss<< snp_pos ;
		 	ss<< "\t"<< match_num ;
		 	ss<< "\t";
		 	double mut_num =0;
		 	for(map<string,double>::iterator m_itr = chkSNPMap.begin(); m_itr != chkSNPMap.end() ; m_itr++){
		 		
		 		ss<<  m_itr->first<< ":" << m_itr->second <<"," ;
		 		mut_num = mut_num+m_itr->second;
		 	}
		 	
		 	outVec.push_back(ss.str());
			
		}
		
	}
	
	tio.dlwrite(out_fname, outVec);
	
}

void SNP::listAllOvSNPs_VCF( vector<string> ChrSNPVec , std::string pk_fname , std::string mp_fname, std::string out_fname, int pos_pks ,double min_score  ){
	
	stringstream ss;
	GenomeSegment gs ; 
	gs.initGenomeSegment_pksbed(  pk_fname ,INT_MAX	 );
	string str;
	//vector<string> ChrSNPVec = tio.dlread_vector(snp_fname);
	vector<int> PosVec;
	vector<string> ASVec;
	map<int,string> SNPMap ; 
	map<int,string> MPMap ;
	vector<string> outVec;
	for(int i =0 ; i<ChrSNPVec.size(); i++){
		
		vector<string> tokvec = tok.Tokenize(ChrSNPVec[i], "\t");
		int snp_pos = atoi(tokvec[1].c_str());
		string alt_str = tokvec[4];
		
		if(alt_str.compare("a")==0){
			alt_str= "A" ;
		}
		if(alt_str.compare("t")==0){
			alt_str= "T" ;
		}
		if(alt_str.compare("c")==0){
			alt_str= "C" ;
		}
		if(alt_str.compare("g")==0){
			alt_str= "G" ;
		}
		
		//snp_pos ++;
		
		string idx_str = gs.find_interval(snp_pos);
		vector<string> tmpvec = tok.Tokenize(idx_str, "_");
		int a_idx = atoi(tmpvec[0].c_str());
		int b_idx = atoi(tmpvec[1].c_str());
		ss.str("");
		ss<< gs.PosVec[a_idx]<< "_"<< gs.PosVec[b_idx];
		if( gs.ExonMap[ss.str()].size() >0 ){
			SNPMap.insert(make_pair( snp_pos ,alt_str ));
		} 
	}	
	
	fstream filestr;
	filestr.open(mp_fname.c_str(), ios_base::in);
	
	while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 	
			vector<string> tokvec = tok.Tokenize( str,"\t");					
			int mp_pos = atoi(tokvec[1].c_str()); 
			if(atoi(tokvec[3].c_str())>0){
				MPMap.insert( make_pair(mp_pos, str));
			}
		}
	}
	filestr.close();
	
	for(map<int,string>::iterator m_itr = SNPMap.begin() ; m_itr != SNPMap.end(); m_itr++){
		//cout << "look at "<< m_itr->first << "\t" << m_itr->second<< endl;
		int snp_pos = m_itr->first;
		string alt_str = m_itr->second;
		map<string, double> chkSNPMap ; 
		
		chkSNPMap.insert(make_pair("A" , 0.0));
		chkSNPMap.insert(make_pair("T" , 0.0));
		chkSNPMap.insert(make_pair("C" , 0.0));
		chkSNPMap.insert(make_pair("G" , 0.0));
		
		double match_num =0;
		
		if(MPMap.find(snp_pos) != MPMap.end()){
			//cout << MPMap[snp_pos] << endl;
		 	vector<string> tmpvec  = tok.Tokenize( MPMap[snp_pos] , "\t") ; 
		 	string ref_str = tmpvec[2]; 
		 	string  base_str = tmpvec[4];
		 	string  score_str = tmpvec[5];
		 	vector<string> baseVec = _sepBaseStr(base_str);
		 	
		 	if(ref_str.compare("a")==0){
		 		ref_str="A";
		 	}
		 	if(ref_str.compare("t")==0){
		 		ref_str="T";
		 	}
		 	if(ref_str.compare("c")==0){
		 		ref_str="C";
		 	}
		 	if(ref_str.compare("g")==0){
		 		ref_str="G";
		 	}
		 	
		 	for( int j =0 ; j<baseVec.size() ; j++){
		 		int score =(int) score_str[j]; 
		 		if(score> min_score){
		 			if(baseVec[j].compare(".")==0 || baseVec[j].compare(",")==0 ){
		 				match_num = match_num +1.0;
		 				chkSNPMap[ref_str] = chkSNPMap[ref_str] + 1.0;
		 			}
		 			else{
		 				string mut_str ;
		 				if(baseVec[j].compare("A")==0 || baseVec[j].compare("a")==0 ){
		 					mut_str="A";
		 				}
		 			
		 				if(baseVec[j].compare("T")==0 || baseVec[j].compare("t")==0 ){
		 					mut_str="T";
		 				}
		 			
		 				if(baseVec[j].compare("C")==0 || baseVec[j].compare("c")==0 ){
		 					mut_str="C";
		 				}
		 				if(baseVec[j].compare("G")==0 || baseVec[j].compare("g")==0 ){
		 					mut_str="G";
		 				}
		 				
		 				if(mut_str.compare(alt_str) == 0){
		 					chkSNPMap[alt_str] = chkSNPMap[alt_str] + 1.0;;
		 				}
		 				
		 			}	
		 		}
		 		
		 	}
		 	ss.str("");
		 	ss<< snp_pos ;
		 	ss<< "\t"<< match_num ;
		 	ss<< "\t";
		 	double mut_num =0;
		 	for(map<string,double>::iterator m_itr = chkSNPMap.begin(); m_itr != chkSNPMap.end() ; m_itr++){
		 		
		 		ss<<  m_itr->first<< ":" << m_itr->second <<"," ;
		 		mut_num = mut_num+m_itr->second;
		 	}
		 	
		 	outVec.push_back(ss.str());
			
		}
		
	}
	
	tio.dlwrite(out_fname, outVec);
	
}

void SNP::listAllOvSNPs_denovo_lite( vector<string> ChrSNPVec , std::string pk_fname , std::string out_fname){
	
	stringstream ss;
	GenomeSegment gs ; 
	gs.initGenomeSegment_pksbed(  pk_fname ,INT_MAX	 );
	string str;
	//vector<string> ChrSNPVec = tio.dlread_vector(snp_fname);
	
	vector<string> outVec;
	for(int i =0 ; i<ChrSNPVec.size(); i++){
		
		vector<string> tokvec = tok.Tokenize(ChrSNPVec[i], " ");
		int snp_pos = atoi(tokvec[1].c_str());
		
		string idx_str = gs.find_interval(snp_pos);
		vector<string> tmpvec = tok.Tokenize(idx_str, "_");
		int a_idx = atoi(tmpvec[0].c_str());
		int b_idx = atoi(tmpvec[1].c_str());
		ss.str("");
		ss<< gs.PosVec[a_idx]<< "_"<< gs.PosVec[b_idx];
		if( gs.ExonMap[ss.str()].size() >0 ){
			outVec.push_back(ChrSNPVec[i]);
//			SNPMap.insert(make_pair( snp_pos ,alt_str ));
		} 
	}	
	
	tio.dlwrite(out_fname, outVec);
	
}

void SNP::listAllOvSNPs_VCF_lite( vector<string> ChrSNPVec , std::string pk_fname , std::string out_fname){
	
	stringstream ss;
	GenomeSegment gs ; 
	gs.initGenomeSegment_pksbed(  pk_fname ,INT_MAX	 );
	string str;
	//vector<string> ChrSNPVec = tio.dlread_vector(snp_fname);
	
	vector<string> outVec;
	for(int i =0 ; i<ChrSNPVec.size(); i++){
		
		vector<string> tokvec = tok.Tokenize(ChrSNPVec[i], "\t");
		int snp_pos = atoi(tokvec[1].c_str());
		
		string idx_str = gs.find_interval(snp_pos);
		vector<string> tmpvec = tok.Tokenize(idx_str, "_");
		int a_idx = atoi(tmpvec[0].c_str());
		int b_idx = atoi(tmpvec[1].c_str());
		ss.str("");
		ss<< gs.PosVec[a_idx]<< "_"<< gs.PosVec[b_idx];
		if( gs.ExonMap[ss.str()].size() >0 ){
			outVec.push_back(ChrSNPVec[i]);
//			SNPMap.insert(make_pair( snp_pos ,alt_str ));
		} 
	}	
	
	tio.dlwrite(out_fname, outVec);
	
}


std::vector<std::string>  SNP::normalize( std::string snp_fname,  std::string cls_fname	,  std::string prior_dir  , int pos_flag, double lib_size){
  
  vector<int> ClsPosVec, SNPPosVec;
	vector<string> outVec; 
	stringstream ss;
//read clsposvec 
	ClsPosVec.push_back(0);
	vector<string> CLSVec = tio.dlread_vector(cls_fname);
	for(int i =0 ; i<CLSVec.size() ; i++){
		vector<string> tokvec = tok.Tokenize( CLSVec[i], "\t");
		ClsPosVec.push_back(atoi(tokvec[0].c_str()));
	}
	ClsPosVec.push_back(INT_MAX);
	
//read snpposvec
	vector<string> SNPVec = tio.dlread_vector(snp_fname);
	for(int i =0 ; i<SNPVec.size() ; i++){
		vector<string> tokvec = tok.Tokenize( SNPVec[i], "\t");
		SNPPosVec.push_back(atoi(tokvec[0].c_str()));
	}
	
	vector<int> SNP_OffsetVec = listCLS_pos(  SNPPosVec,   ClsPosVec	, pos_flag );

	//read prior
	map<int,double> APriorMap,TPriorMap,CPriorMap,GPriorMap;
	vector<string> AVec,TVec,CVec,GVec  ;
	if(pos_flag == 1){
		AVec = tio.dlread_vector(prior_dir+"/a_pos.txt");
		TVec = tio.dlread_vector(prior_dir+"/t_pos.txt"); 
		CVec = tio.dlread_vector(prior_dir+"/c_pos.txt");
		GVec = tio.dlread_vector(prior_dir+"/g_pos.txt");
	}
	else{
		AVec = tio.dlread_vector(prior_dir+"/a_neg.txt");
		TVec = tio.dlread_vector(prior_dir+"/t_neg.txt"); 
		CVec = tio.dlread_vector(prior_dir+"/c_neg.txt");
		GVec = tio.dlread_vector(prior_dir+"/g_neg.txt");
	}
	int offset_max,offset_min; 
	
	for(int i =0; i< AVec.size() ; i ++){
		vector<string> tokvec = tok.Tokenize(AVec[i] ,",");
		APriorMap.insert(make_pair(atoi(tokvec[0].c_str()) , atof(tokvec[1].c_str())));
		
		tokvec = tok.Tokenize(TVec[i] ,",");
		TPriorMap.insert(make_pair(atoi(tokvec[0].c_str()) , atof(tokvec[1].c_str())));
		
		tokvec = tok.Tokenize(CVec[i] ,",");
		CPriorMap.insert(make_pair(atoi(tokvec[0].c_str()) , atof(tokvec[1].c_str())));
		
		tokvec = tok.Tokenize(GVec[i] ,",");
		GPriorMap.insert(make_pair(atoi(tokvec[0].c_str()) , atof(tokvec[1].c_str())));
		
		if(i==0){
			offset_min=  atoi(tokvec[0].c_str()); 
		}
		if(i==AVec.size()-1){
			offset_max=  atoi(tokvec[0].c_str()); 
		}
		
	}
	
	//output normalized allele-specific reads count
	for(int i =0 ; i<SNPVec.size(); i++){
		int M,m ;
		vector<double> RCVec ; 
		vector<string> tokvec = tok.Tokenize( SNPVec[i], "\t");
		vector<string> bpvec = tok.Tokenize( tokvec[2], ",");
		ss.str("");
		map<string,double> bpmap; 
		for(int j =0 ; j<bpvec.size(); j++){
			vector<string> tmpvec = tok.Tokenize(bpvec[j],":" );
			bpmap.insert(make_pair(tmpvec[0] , atof(tmpvec[1].c_str()))); 
		}
		int offset = SNP_OffsetVec[i];
		
			
		ss<< SNPVec[i];
		if(pos_flag ==1 ){
			ss<< "\t+";
		}
		else{
			ss<< "\t-";
		}
		
		ss<< "\t" <<  offset ; 
		
		if(offset <= offset_max && offset >= offset_min){
			
			bpmap["A"] = bpmap["A"]/(APriorMap[offset]*lib_size);
			ss<< "\t"<< APriorMap[offset]<< "\t"<< bpmap["A"];
			RCVec.push_back(bpmap["A"]);
			bpmap["T"] = bpmap["T"]/(TPriorMap[offset]*lib_size);
			ss<< "\t"<< TPriorMap[offset]<< "\t"<< bpmap["T"];
			RCVec.push_back(bpmap["T"]);
			bpmap["C"] = bpmap["C"]/(CPriorMap[offset]*lib_size);
			ss<< "\t"<< CPriorMap[offset]<< "\t"<< bpmap["C"];
			RCVec.push_back(bpmap["C"]);
			bpmap["G"] = bpmap["G"]/(GPriorMap[offset]*lib_size);
			ss<< "\t"<< GPriorMap[offset]<< "\t"<< bpmap["G"];
			RCVec.push_back(bpmap["G"]);
		}
		else{
			bpmap["A"] = bpmap["A"]/(0.25*lib_size);
			ss<< "\t"<< 0.25<< "\t"<< bpmap["A"];
			RCVec.push_back(bpmap["A"]);
			bpmap["T"] = bpmap["T"]/(0.25*lib_size);
			ss<< "\t"<< 0.25<< "\t"<< bpmap["T"];
			RCVec.push_back(bpmap["T"]);
			bpmap["C"] = bpmap["C"]/(0.25*lib_size);
			ss<< "\t"<< 0.25<< "\t"<< bpmap["C"];
			RCVec.push_back(bpmap["C"]);
			bpmap["G"] = bpmap["G"]/(0.25*lib_size);
			ss<< "\t"<< 0.25<< "\t"<< bpmap["G"];
			RCVec.push_back(bpmap["G"]);
		}
		
		sort (RCVec.begin(), RCVec.end());
		ss<< "\t"<<RCVec[3]<<"\t"<<RCVec[2];
		outVec.push_back(ss.str());
	}
	
	return outVec; 	
	
}

std::vector<std::string>  SNP::normalize_VCF( std::string snp_fname,  std::string cls_fname	,  std::string prior_dir , int pos_flag , double lib_size){
  
  vector<int> ClsPosVec, SNPPosVec;
	vector<string> outVec; 
	stringstream ss;
//read clsposvec 
	ClsPosVec.push_back(0);
	vector<string> CLSVec = tio.dlread_vector(cls_fname);
	for(int i =0 ; i<CLSVec.size() ; i++){
		vector<string> tokvec = tok.Tokenize( CLSVec[i], "\t");
		ClsPosVec.push_back(atoi(tokvec[0].c_str()));
	}
	
	ClsPosVec.push_back(INT_MAX);
	sort(ClsPosVec.begin(), ClsPosVec.end());
	
//read snpposvec
	vector<string> SNPVec = tio.dlread_vector(snp_fname);
	for(int i =0 ; i<SNPVec.size() ; i++){
		vector<string> tokvec = tok.Tokenize( SNPVec[i], "\t");
		SNPPosVec.push_back(atoi(tokvec[0].c_str()));
	}
	
	vector<int> SNP_OffsetVec = listCLS_pos(  SNPPosVec,   ClsPosVec	, pos_flag );

	//read prior
	map<int,double> APriorMap,TPriorMap,CPriorMap,GPriorMap;
	vector<string> AVec,TVec,CVec,GVec  ;
	AVec = tio.dlread_vector(prior_dir+"/a.txt");
	TVec = tio.dlread_vector(prior_dir+"/t.txt"); 
	CVec = tio.dlread_vector(prior_dir+"/c.txt");
	GVec = tio.dlread_vector(prior_dir+"/g.txt");
	
	
	int offset_max,offset_min; 
	
	for(int i =0; i< AVec.size() ; i ++){
		vector<string> tokvec = tok.Tokenize(AVec[i] ,",");
		APriorMap.insert(make_pair(atoi(tokvec[0].c_str()) , atof(tokvec[1].c_str())));
		
		tokvec = tok.Tokenize(TVec[i] ,",");
		TPriorMap.insert(make_pair(atoi(tokvec[0].c_str()) , atof(tokvec[1].c_str())));
		
		tokvec = tok.Tokenize(CVec[i] ,",");
		CPriorMap.insert(make_pair(atoi(tokvec[0].c_str()) , atof(tokvec[1].c_str())));
		
		tokvec = tok.Tokenize(GVec[i] ,",");
		GPriorMap.insert(make_pair(atoi(tokvec[0].c_str()) , atof(tokvec[1].c_str())));
		
		if(i==0){
			offset_min=  atoi(tokvec[0].c_str()); 
			//cout << "offset_min:" << tokvec[0] << endl;
		}
		if(i==AVec.size()-1){
			offset_max=  atoi(tokvec[0].c_str()); 
			//cout << "offset_max:" << tokvec[0] << endl;
		}
		
	}
	
	//output normalized allele-specific reads count
	for(int i =0 ; i<SNPVec.size(); i++){
		int M,m ;
		vector<double> RCVec ; 
		vector<string> tokvec = tok.Tokenize( SNPVec[i], "\t");
		string ref_str = tokvec[1];
		string alt_str = tokvec[2];
		vector<string> bpvec = tok.Tokenize( tokvec[3], ",");
		ss.str("");
		map<string,double> bpmap; 
		for(int j =0 ; j<bpvec.size(); j++){
			vector<string> tmpvec = tok.Tokenize(bpvec[j],":" );
			bpmap.insert(make_pair(tmpvec[0] , atof(tmpvec[1].c_str()))); 
		}
		int offset = SNP_OffsetVec[i];
		
			
		ss<< SNPVec[i];
		if(pos_flag ==1 ){
			ss<< "\t+";
		}
		else{
			ss<< "\t-";
		}
		
		ss<< "\t" <<  offset ; 
		
		if(offset <= offset_max && offset >= offset_min){
			//cout << "kk :"<<  SNPVec[i]<< "," << offset_max<< ","<<offset_min << endl;
			if(pos_flag ==1 ){
				bpmap["A"] = bpmap["A"]/(APriorMap[offset]*lib_size);
				ss<< "\t"<< APriorMap[offset]<< "\t"<< bpmap["A"];
				RCVec.push_back(bpmap["A"]);
				bpmap["T"] = bpmap["T"]/(TPriorMap[offset]*lib_size);
				ss<< "\t"<< TPriorMap[offset]<< "\t"<< bpmap["T"];
				RCVec.push_back(bpmap["T"]);
				bpmap["C"] = bpmap["C"]/(CPriorMap[offset]*lib_size);
				ss<< "\t"<< CPriorMap[offset]<< "\t"<< bpmap["C"];
				RCVec.push_back(bpmap["C"]);
				bpmap["G"] = bpmap["G"]/(GPriorMap[offset]*lib_size);
				ss<< "\t"<< GPriorMap[offset]<< "\t"<< bpmap["G"];
				RCVec.push_back(bpmap["G"]);
			}
			else{
				bpmap["A"] = bpmap["A"]/(TPriorMap[offset]*lib_size);
				ss<< "\t"<< TPriorMap[offset]<< "\t"<< bpmap["A"];
				RCVec.push_back(bpmap["A"]);
				bpmap["T"] = bpmap["T"]/(APriorMap[offset]*lib_size);
				ss<< "\t"<< APriorMap[offset]<< "\t"<< bpmap["T"];
				RCVec.push_back(bpmap["T"]);
				bpmap["C"] = bpmap["C"]/(GPriorMap[offset]*lib_size);
				ss<< "\t"<< GPriorMap[offset]<< "\t"<< bpmap["C"];
				RCVec.push_back(bpmap["C"]);
				bpmap["G"] = bpmap["G"]/(CPriorMap[offset]*lib_size);
				ss<< "\t"<< CPriorMap[offset]<< "\t"<< bpmap["G"];
				RCVec.push_back(bpmap["G"]);

			}
		}
		else{
			bpmap["A"] = bpmap["A"]/(0.25*lib_size);
			ss<< "\t"<< 0.25<< "\t"<< bpmap["A"];
			RCVec.push_back(bpmap["A"]);
			bpmap["T"] = bpmap["T"]/(0.25*lib_size);
			ss<< "\t"<< 0.25<< "\t"<< bpmap["T"];
			RCVec.push_back(bpmap["T"]);
			bpmap["C"] = bpmap["C"]/(0.25*lib_size);
			ss<< "\t"<< 0.25<< "\t"<< bpmap["C"];
			RCVec.push_back(bpmap["C"]);
			bpmap["G"] = bpmap["G"]/(0.25*lib_size);
			ss<< "\t"<< 0.25<< "\t"<< bpmap["G"];
			RCVec.push_back(bpmap["G"]);
		}
		
		//sort (RCVec.begin(), RCVec.end());
			
		ss<< "\t"<<bpmap[ref_str]<<"\t"<<bpmap[alt_str];
			
		
		outVec.push_back(ss.str());
	}
	
	return outVec; 	
	
}

std::vector<std::string>  SNP::normalize_denovo( std::string snp_fname,  std::string cls_fname	,  std::string prior_dir , int pos_flag , double lib_size){
  
  vector<int> ClsPosVec, SNPPosVec;
	vector<string> outVec; 
	stringstream ss;
//read clsposvec 
	ClsPosVec.push_back(0);
	vector<string> CLSVec = tio.dlread_vector(cls_fname);
	for(int i =0 ; i<CLSVec.size() ; i++){
		vector<string> tokvec = tok.Tokenize( CLSVec[i], "\t");
		ClsPosVec.push_back(atoi(tokvec[0].c_str()));
	}
	
	ClsPosVec.push_back(INT_MAX);
	sort(ClsPosVec.begin(), ClsPosVec.end());
	
//read snpposvec
	vector<string> SNPVec = tio.dlread_vector(snp_fname);
	for(int i =0 ; i<SNPVec.size() ; i++){
		vector<string> tokvec = tok.Tokenize( SNPVec[i], "\t");
		SNPPosVec.push_back(atoi(tokvec[0].c_str()));
	}
	
	vector<int> SNP_OffsetVec = listCLS_pos(  SNPPosVec,   ClsPosVec	, pos_flag );

	//read prior
	map<int,double> APriorMap,TPriorMap,CPriorMap,GPriorMap;
	vector<string> AVec,TVec,CVec,GVec  ;
	AVec = tio.dlread_vector(prior_dir+"/a.txt");
	TVec = tio.dlread_vector(prior_dir+"/t.txt"); 
	CVec = tio.dlread_vector(prior_dir+"/c.txt");
	GVec = tio.dlread_vector(prior_dir+"/g.txt");
	
	
	int offset_max,offset_min; 
	
	for(int i =0; i< AVec.size() ; i ++){
		vector<string> tokvec = tok.Tokenize(AVec[i] ,",");
		APriorMap.insert(make_pair(atoi(tokvec[0].c_str()) , atof(tokvec[1].c_str())));
		
		tokvec = tok.Tokenize(TVec[i] ,",");
		TPriorMap.insert(make_pair(atoi(tokvec[0].c_str()) , atof(tokvec[1].c_str())));
		
		tokvec = tok.Tokenize(CVec[i] ,",");
		CPriorMap.insert(make_pair(atoi(tokvec[0].c_str()) , atof(tokvec[1].c_str())));
		
		tokvec = tok.Tokenize(GVec[i] ,",");
		GPriorMap.insert(make_pair(atoi(tokvec[0].c_str()) , atof(tokvec[1].c_str())));
		
		if(i==0){
			offset_min=  atoi(tokvec[0].c_str()); 
			//cout << "offset_min:" << tokvec[0] << endl;
		}
		if(i==AVec.size()-1){
			offset_max=  atoi(tokvec[0].c_str()); 
			//cout << "offset_max:" << tokvec[0] << endl;
		}
		
	}
	
	//output normalized allele-specific reads count
	for(int i =0 ; i<SNPVec.size(); i++){
		int M,m ;
		vector<double> RCVec ; 
		vector<string> tokvec = tok.Tokenize( SNPVec[i], "\t");
		string ref_str = tokvec[1];
		string alt_str = tokvec[2];
		vector<string> bpvec = tok.Tokenize( tokvec[3], ",");
		ss.str("");
		map<string,double> bpmap; 
		for(int j =0 ; j<bpvec.size(); j++){
			vector<string> tmpvec = tok.Tokenize(bpvec[j],":" );
			bpmap.insert(make_pair(tmpvec[0] , atof(tmpvec[1].c_str()))); 
		}
		int offset = SNP_OffsetVec[i];
		
			
		ss<< SNPVec[i];
		if(pos_flag ==1 ){
			ss<< "\t+";
		}
		else{
			ss<< "\t-";
		}
		
		ss<< "\t" <<  offset ; 
		
		if(offset <= offset_max && offset >= offset_min){
			//cout << "kk :"<<  SNPVec[i]<< "," << offset_max<< ","<<offset_min << endl;
			if(pos_flag ==1 ){
				bpmap["A"] = bpmap["A"]/(APriorMap[offset]*lib_size);
				ss<< "\t"<< APriorMap[offset]<< "\t"<< bpmap["A"];
				RCVec.push_back(bpmap["A"]);
				bpmap["T"] = bpmap["T"]/(TPriorMap[offset]*lib_size);
				ss<< "\t"<< TPriorMap[offset]<< "\t"<< bpmap["T"];
				RCVec.push_back(bpmap["T"]);
				bpmap["C"] = bpmap["C"]/(CPriorMap[offset]*lib_size);
				ss<< "\t"<< CPriorMap[offset]<< "\t"<< bpmap["C"];
				RCVec.push_back(bpmap["C"]);
				bpmap["G"] = bpmap["G"]/(GPriorMap[offset]*lib_size);
				ss<< "\t"<< GPriorMap[offset]<< "\t"<< bpmap["G"];
				RCVec.push_back(bpmap["G"]);
			}
			else{
				bpmap["A"] = bpmap["A"]/(TPriorMap[offset]*lib_size);
				ss<< "\t"<< TPriorMap[offset]<< "\t"<< bpmap["A"];
				RCVec.push_back(bpmap["A"]);
				bpmap["T"] = bpmap["T"]/(APriorMap[offset]*lib_size);
				ss<< "\t"<< APriorMap[offset]<< "\t"<< bpmap["T"];
				RCVec.push_back(bpmap["T"]);
				bpmap["C"] = bpmap["C"]/(GPriorMap[offset]*lib_size);
				ss<< "\t"<< GPriorMap[offset]<< "\t"<< bpmap["C"];
				RCVec.push_back(bpmap["C"]);
				bpmap["G"] = bpmap["G"]/(CPriorMap[offset]*lib_size);
				ss<< "\t"<< CPriorMap[offset]<< "\t"<< bpmap["G"];
				RCVec.push_back(bpmap["G"]);

			}
		}
		else{
			bpmap["A"] = bpmap["A"]/(0.25*lib_size);
			ss<< "\t"<< 0.25<< "\t"<< bpmap["A"];
			RCVec.push_back(bpmap["A"]);
			bpmap["T"] = bpmap["T"]/(0.25*lib_size);
			ss<< "\t"<< 0.25<< "\t"<< bpmap["T"];
			RCVec.push_back(bpmap["T"]);
			bpmap["C"] = bpmap["C"]/(0.25*lib_size);
			ss<< "\t"<< 0.25<< "\t"<< bpmap["C"];
			RCVec.push_back(bpmap["C"]);
			bpmap["G"] = bpmap["G"]/(0.25*lib_size);
			ss<< "\t"<< 0.25<< "\t"<< bpmap["G"];
			RCVec.push_back(bpmap["G"]);
		}
		
		//sort (RCVec.begin(), RCVec.end());
		
		outVec.push_back(ss.str());
	}
	
	
	
	return outVec; 	
	
}

std::vector<int> SNP::listCLS_pos( std::vector<int> SNPPosVec,  std::vector<int> ClsPosVec	,int pos_flag ){

	vector<int> outVec;
	stringstream ss;
	stringstream int_ss;

	for(int i=0 ; i<SNPPosVec.size(); i++ ){
		int snp_pos = SNPPosVec[i];
		
		int a = 0;
		int b = ClsPosVec.size()-1;	
		
		while( (b-a)>1 ){
			int mid_idx =(b+a)/2; 
			if(ClsPosVec[mid_idx]<snp_pos){
				a = mid_idx;
			}
			else{
				b= mid_idx;
			}
		}
		int cls_pos_a, cls_pos_b;
		if(pos_flag == 1 ){
			cls_pos_a = ClsPosVec[a]-1;
			cls_pos_b = ClsPosVec[b]-1;
		}
		else{
			cls_pos_a = ClsPosVec[a]+1;
			cls_pos_b = ClsPosVec[b]+1;
		}
		int offset_a = snp_pos - cls_pos_a;
		int offset_b = snp_pos - cls_pos_b;
		int dist_a = abs(offset_a);	
		int dist_b = abs(offset_b);
		
		if(dist_a < dist_b){
			outVec.push_back(offset_a);
		}
		else{
			outVec.push_back(offset_b);
		}

	}
	
	return outVec;	
}

void SNP::listInputSNPs_depth( std::string snp_fname,  std::string mp_fname, std::string out_fname ,int min_score ){
	
	stringstream ss;
  vector<string> outVec;
	vector<string> ChrSNPVec = tio.dlread_vector(snp_fname);
	map<int,string> SNPMap ; 
	map<int,string> MPMap ;
	
	for(int i =0 ; i<ChrSNPVec.size(); i++){
		
		vector<string> tokvec = tok.Tokenize(ChrSNPVec[i], "\t");
		int snp_pos = atoi(tokvec[0].c_str());
		SNPMap.insert(make_pair(snp_pos,"")); 	
	}	
	
	vector<string> MPVec = tio.dlread_vector(mp_fname);
	for(int i =0 ; i<MPVec.size(); i++){
			vector<string> tokvec = tok.Tokenize(MPVec[i],"\t");
			int mp_pos = atoi(tokvec[1].c_str()); 
			if(atoi(tokvec[3].c_str())>0){
				MPMap.insert( make_pair(mp_pos, MPVec[i]));
			}
	}
	
	for(map<int,string>::iterator m_itr = SNPMap.begin() ; m_itr != SNPMap.end(); m_itr++){
		
		int snp_pos = m_itr->first;
		double depnum =0;
//		
		if(MPMap.find(snp_pos) != MPMap.end()){
			

		 	vector<string> tmpvec  = tok.Tokenize( MPMap[snp_pos] , "\t") ; 
		 	string  base_str = tmpvec[4];
		 	string  score_str = tmpvec[5];
		 	vector<string> baseVec = _sepBaseStr(base_str);
		 	
		 	for( int j =0 ; j<baseVec.size() ; j++){
		 		int score =(int) score_str[j]; 
		 		if(score> min_score){
		 			if(baseVec[j].compare(".")==0 || baseVec[j].compare(",")==0 ){
		 				depnum ++;
		 			}
		 			if(baseVec[j].compare("A")==0 || baseVec[j].compare("a")==0 ){
		 				depnum ++;
		 			}
		 			
		 			if(baseVec[j].compare("T")==0 || baseVec[j].compare("t")==0 ){
		 				depnum ++;
		 			}
		 			
		 			if(baseVec[j].compare("C")==0 || baseVec[j].compare("c")==0 ){
		 				depnum ++;
		 			}
		 			if(baseVec[j].compare("G")==0 || baseVec[j].compare("g")==0 ){
		 				depnum ++;
		 			}
		 				
		 		}
		 		
		 	}
		}
		
		ss.str("");
		ss<< snp_pos << "\t"<< depnum;
		outVec.push_back(ss.str());
	}
	
	tio.dlwrite(out_fname, outVec);
	
}

void SNP::SNPtoRC(	std::string SNP_dir, 	std::string tmp_dir, std::string Rep_prefix , std::string chr_str, double lib_size){
	stringstream ss;
	TextIO tio;
	TOK tok;
	vector<string> outVec; 
	
	vector<string> chrVec = tok.Tokenize( chr_str, ",") ; 
	for(int i =0 ; i<chrVec.size(); i++){
		string chr_name= chrVec[i];
		cout << "processing "<<Rep_prefix <<" SNPS in "<< chr_name<< endl;
	
	  listAllOvSNPs( SNP_dir+"/"+chr_name+".txt", tmp_dir+"/"+Rep_prefix+"_"+chr_name+"_pos.bed", tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+".mp", tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+".snp", 1,0 );
		listAllOvSNPs( SNP_dir+"/"+chr_name+".txt", tmp_dir+"/"+Rep_prefix+"_"+chr_name+"_neg.bed", tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+".mp", tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+".snp", -1,0 );
		
		//add SNP_id
//		
//		vector<string> snpvec = tio.dlread_vector(SNP_dir+"/"+chr_name+".txt");
//		map<string,string> snpmap;
//		for(int j=0 ; j<snpvec.size(); j++){
//			vector<string> tmpvec = tok.Tokenize(snpvec[j] ,"\t" );
//			snpmap.insert(make_pair(tmpvec[2] , tmpvec[3]) );
//		}
//			
		vector<string> SNP_RepVec= normalize( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+".snp",  tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls", tmp_dir , 1, lib_size);
		for(int j=0 ; j<SNP_RepVec.size(); j++){
			vector<string> tmpvec = tok.Tokenize( SNP_RepVec[j] , "\t");
			ss.str("");
			ss<< chr_name<<"_"<<tmpvec[0] << "\t"<< chr_name << "\t"<< SNP_RepVec[j];
			outVec.push_back(ss.str());
		}
//		
		SNP_RepVec= normalize( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+".snp",  tmp_dir+"/"+Rep_prefix+"_neg_2_"+chr_name+".cls", tmp_dir , 0 , lib_size);
		
	
		for(int j=0 ; j<SNP_RepVec.size(); j++){
			vector<string> tmpvec = tok.Tokenize( SNP_RepVec[j] , "\t");
  		ss.str("");
			ss<<chr_name<<"_"<<tmpvec[0] << "\t" <<  chr_name << "\t"<< SNP_RepVec[j];
			outVec.push_back(ss.str());
		}
	}
	tio.dlwrite( tmp_dir+"/"+Rep_prefix+".rc" , outVec);
	
}

void SNP::SNPtoRC_denovo(	std::string SNP_dir, 	std::string tmp_dir,  std::vector< std::string > RepIDVec, std::string chr_str, std::map<std::string, double> SizeMap){
	stringstream ss;
	TextIO tio;
	TOK tok;
	vector<string> outVec; 
	
	vector<string> chrVec = tok.Tokenize( chr_str, ",") ; 
	vector< vector<double> > cutoff_vec;
	//read seq cut-off 
	vector<string> ctvec = tio.dlread_vector(tmp_dir+"/seqerr.cut");
	
	//init cut-off vec 
	 
	for(int i =0 ; i<ctvec.size(); i++){
		vector<double> tmpvec;
		cutoff_vec.push_back(tmpvec);

	}
	
	for(int i =0 ; i<ctvec.size(); i++){
		vector<string> tokvec = tok.Tokenize(ctvec[i] , "\t");
		for(int j =0 ;j<tokvec.size() ; j++){
			double ct_value = atof(tokvec[j].c_str());
			cutoff_vec[i].push_back(ct_value); 
		} 
	}
	
	//read and normalize raw RC
	
	for(int i =0 ; i<chrVec.size(); i++){
		string chr_name= chrVec[i];
		//read SNP annotation 
		vector<string> SNPVec = tio.dlread_vector(SNP_dir+"/"+chr_name+".txt");
		// ov annotated all SNPs to peak  
	  cout << "read " <<SNP_dir+"/"+chr_name+".txt" <<endl;	  
	  for(int j=0 ; j<RepIDVec.size();j++){
  	
  		string Rep_prefix = RepIDVec[j];  
//			
			listAllOvSNPs_denovo_lite( SNPVec, tmp_dir+"/"+Rep_prefix+"_"+chr_name+"_pos.bed", tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.snp" );			
			cout << "print : "<< tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.snp" << endl;
			listAllOvSNPs_denovo_lite( SNPVec, tmp_dir+"/"+Rep_prefix+"_"+chr_name+"_neg.bed", tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.snp"  );
			cout << "print : "<< tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.snp" << endl;
		}
		
		
		map<string, string> possnpmap,negsnpmap ;
		
//		
	  for(int j=0 ; j<RepIDVec.size();j++){
  	
  		string Rep_prefix = RepIDVec[j];  
			
			vector<string> posvec = tio.dlread_vector(tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.snp");
			vector<string> negvec = tio.dlread_vector(tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.snp");
			
			for(int k =0 ; k<posvec.size() ; k++){
				
				if(possnpmap.find(posvec[k]) == possnpmap.end()){
					possnpmap.insert(make_pair( posvec[k], ""));
				}
			}
			
			for(int k =0 ; k<negvec.size() ; k++){
				
				if(negsnpmap.find(negvec[k]) == negsnpmap.end()){
					negsnpmap.insert(make_pair( negvec[k], ""));
				}
			}
			
		}
////		
		//print the union of SNP 
		vector<string> univec; 
		for(map<string, string>::iterator mitr = possnpmap.begin(); mitr != possnpmap.end(); mitr++){
			univec.push_back(mitr->first);
		}
		
		tio.dlwrite( tmp_dir+"/"+chr_name+"_pos_allovdenovo.snp" , univec);
		cout << "write: "<< tmp_dir+"/"+chr_name+"_pos_allovdenovo.snp"<< endl;
		univec.clear();
		
		for(map<string, string>::iterator mitr = negsnpmap.begin(); mitr != negsnpmap.end(); mitr++){
			univec.push_back(mitr->first);
		}
		
		tio.dlwrite( tmp_dir+"/"+chr_name+"_neg_allovdenovo.snp" , univec);
		cout << "write: "<< tmp_dir+"/"+chr_name+"_neg_allovdenovo.snp"<< endl;
		univec.clear();
		
		
		
		//read read cov
		for(int j=0 ; j<RepIDVec.size();j++){
			string Rep_prefix = RepIDVec[j];
			vector<string> SNPvec = tio.dlread_vector(tmp_dir+"/"+chr_name+"_pos_allovdenovo.snp");
			readRawRC_denovo(SNPvec, tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+".mp" , tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawdenovo.rc");
			cout << "calc raw RC :" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawdenovo.rc" << endl;
		
		//nomalize raw read counts			
			vector<string>  SNP_RepVec= normalize_denovo( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawdenovo.rc",  tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls", tmp_dir , 1 , SizeMap[Rep_prefix]);
			tio.dlwrite( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.rc" , SNP_RepVec);
			cout << "normalize raw RC :" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.rc\t"<< tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls" << "\t"<<SizeMap[Rep_prefix] << endl;
////			
			SNPvec = tio.dlread_vector(tmp_dir+"/"+chr_name+"_neg_allovdenovo.snp");
			readRawRC_denovo(SNPvec, tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+".mp" , tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawdenovo.rc");
			cout << "calc raw RC :" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawdenovo.rc" << endl;
////			
			SNP_RepVec= normalize_denovo( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawdenovo.rc",  tmp_dir+"/"+Rep_prefix+"_neg_2_"+chr_name+".cls", tmp_dir , 0 , SizeMap[Rep_prefix]);
			tio.dlwrite( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.rc" , SNP_RepVec);
			cout << "normalize raw RC :" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.rc\t"<< tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls" << "\t"<<SizeMap[Rep_prefix] << endl;
		}
		
	}
		
	//for regression training 
	
	
		vector<string> regvec ; 
		regvec.push_back("id\tchr\tpos\tstrand\tallele\txmean\tvar\tcv2");
		vector<string> testvec ; 
	  testvec.push_back("id\tchr\tpos\tstrand\traw_str\tcov\tM_str\tm_str\tM\tm");
	   
		for(int i =0 ; i<chrVec.size(); i++){
			
			string chr_name = chrVec[i];
			//merge RCs of pos SNPs
			
			vector<vector<double> > posNRC_Avec,posNRC_Tvec,posNRC_Cvec,posNRC_Gvec ;
			vector<vector<double> > posmeanvec;
			vector<vector<double > > posvarvec;
			vector<vector<string > > posRawRCvec;
			vector<string > posMstrvec,posmstrvec;
			
			vector<string > posposvec; 
			vector<int> posflagvec;
			
			for(int j=0 ; j<RepIDVec.size();j++){
				string Rep_prefix = RepIDVec[j];	
				vector<string>	rcvec = tio.dlread_vector( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.rc" )	;	
				cout << "merge:" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.rc" << endl;
				if(j==0){
//				//init the xmean vectors
					for(int q =0 ; q<rcvec.size(); q ++){
						vector<double> tmpAvec; 
						posNRC_Avec.push_back(tmpAvec);
						vector<double> tmpTvec; 
						posNRC_Tvec.push_back(tmpTvec);
						vector<double> tmpCvec; 
						posNRC_Cvec.push_back(tmpCvec);
						vector<double> tmpGvec; 
						posNRC_Gvec.push_back(tmpGvec);
						vector<double> tmpvec;
						posmeanvec.push_back(tmpvec);
						vector<double> tmpvarvec;
						posvarvec.push_back(tmpvarvec);
						
						posflagvec.push_back(0);
						
						vector<string> qqvec;
						posRawRCvec.push_back(qqvec);
						
					}	
				}
				for(int k =0 ; k<rcvec.size() ; k++){
					vector<string> tokvec = tok.Tokenize(rcvec[k], "\t" );
					posposvec.push_back(tokvec[0]);
					
					double d_a,d_t,d_c,d_g ;
					d_a = atof(tokvec[7].c_str());
					d_t = atof(tokvec[9].c_str());
					d_c = atof(tokvec[11].c_str());
					d_g = atof(tokvec[13].c_str());
					posNRC_Avec[k].push_back(d_a);
					posNRC_Tvec[k].push_back(d_t);
					posNRC_Cvec[k].push_back(d_c);
					posNRC_Gvec[k].push_back(d_g);
					
					//check if a hetero SNP in the rep
					
					string rawstr = tokvec[3] ;
					posRawRCvec[k].push_back(rawstr);
					vector<string> tmpvec = tok.Tokenize(rawstr, "," );
					vector<double> rawrcvec;
					double dep =0;
					for(int q=0 ; q<4; q++){
						vector<string> tmp1vec = tok.Tokenize(tmpvec[q], ":" );
						rawrcvec.push_back(atof(tmp1vec[1].c_str()));
						dep +=atof(tmp1vec[1].c_str() );
					}
					sort (rawrcvec.begin(), rawrcvec.end());
		    	int dep_idx = dep/10;
		    	if(dep_idx>=9){
		    		dep_idx =9;
		    	}
		    	double M_raw = rawrcvec[3];
		    	double m_raw = rawrcvec[2];
					if((m_raw/(m_raw+M_raw)) >= cutoff_vec[dep_idx][j] ){
						posflagvec[k] = posflagvec[k] +1;
						
					}
					
				}
				
				
				
				
					
			}
			
			//calc mean of ATCG
			
			for(int j =0 ; j< posNRC_Avec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Avec[j].size(); k++){
					sum = sum + posNRC_Avec[j][k] ;
				}
				posmeanvec[j].push_back(sum/(double)posNRC_Avec[j].size());
			}
			
			for(int j =0 ; j< posNRC_Tvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Tvec[j].size(); k++){
					sum = sum + posNRC_Tvec[j][k] ;
				}
				posmeanvec[j].push_back(sum/(double)posNRC_Tvec[j].size());
			}
			
			for(int j =0 ; j< posNRC_Cvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Cvec[j].size(); k++){
					sum = sum + posNRC_Cvec[j][k] ;
				}
				posmeanvec[j].push_back(sum/(double)posNRC_Cvec[j].size());
			}
			
			for(int j =0 ; j< posNRC_Gvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Gvec[j].size(); k++){
					sum = sum + posNRC_Gvec[j][k] ;
				}
				posmeanvec[j].push_back(sum/(double)posNRC_Gvec[j].size());
			}
			
			//calc var of ATCG
			
			
			for(int j =0 ; j< posNRC_Avec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Avec[j].size(); k++){
					sum = sum + (posNRC_Avec[j][k]-posmeanvec[j][0])*(posNRC_Avec[j][k]-posmeanvec[j][0]) ;
				}
				posvarvec[j].push_back(sum/(double)(posNRC_Avec[j].size()-1) );
			}
			
			for(int j =0 ; j< posNRC_Tvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Tvec[j].size(); k++){
					sum = sum + (posNRC_Tvec[j][k]-posmeanvec[j][1])*(posNRC_Tvec[j][k]-posmeanvec[j][1]) ;
				}
				posvarvec[j].push_back(sum/(double)(posNRC_Tvec[j].size()-1) );
			}
			
			for(int j =0 ; j< posNRC_Cvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Cvec[j].size(); k++){
					sum = sum + (posNRC_Cvec[j][k]-posmeanvec[j][2])*(posNRC_Cvec[j][k]-posmeanvec[j][2]) ;
				}
				posvarvec[j].push_back(sum/(double)(posNRC_Cvec[j].size()-1) );
			}
			
			for(int j =0 ; j< posNRC_Gvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Gvec[j].size(); k++){
					sum = sum + (posNRC_Gvec[j][k]-posmeanvec[j][3])*(posNRC_Gvec[j][k]-posmeanvec[j][3]) ;
				}
				posvarvec[j].push_back(sum/(double)(posNRC_Gvec[j].size()-1) );
			}
			
			//determine major and minor allele
			
			for(int j =0 ; j< posmeanvec.size() ; j++){
						
				vector<double> tmpmeanvec = posmeanvec[j];
				sort (tmpmeanvec.begin(), tmpmeanvec.end());
		    double M_rc = tmpmeanvec[3];
		    double m_rc = tmpmeanvec[2];
				
				//search M idx ;
				int Midx =-1;
				int midx =-1;
				for(int k =0 ; k< posmeanvec[j].size() ; k++){
					if(posmeanvec[j][k] == M_rc){
						Midx = k;
					}
					if(posmeanvec[j][k] == m_rc){
						midx = k;
					}
				}
				
				if(Midx ==0 ){
					posMstrvec.push_back("A");
				}else if(Midx ==1 ){
					posMstrvec.push_back("T");
				}else if(Midx ==2 ){
					posMstrvec.push_back("C");
				}else if(Midx ==3 ){
					posMstrvec.push_back("G");
				}else{
					posMstrvec.push_back("-");
				}
				
				
				if(midx ==0 ){
					posmstrvec.push_back("A");
				}else if(midx ==1 ){
					posmstrvec.push_back("T");
				}else if(midx ==2 ){
					posmstrvec.push_back("C");
				}else if(midx ==3 ){
					posmstrvec.push_back("G");
				}else{
					posmstrvec.push_back("-");
				}
				
				//select for reg
				if(Midx > -1){
					if(posmeanvec[j][Midx]>2.0){
						ss.str("");
						ss<< chr_name << "_"<< posposvec[j]<< "_+" << "\t"<< chr_name<< "\t" << posposvec[j] << "\t+\t" << posMstrvec[j] <<  "\t"<< posmeanvec[j][Midx] << "\t" <<posvarvec[j][Midx] <<"\t"<<posvarvec[j][Midx]/(posmeanvec[j][Midx]*posmeanvec[j][Midx])  ;
						regvec.push_back(ss.str());
					}
				}
				if(midx > -1){
					if(posmeanvec[j][midx]>2.0){
						ss.str("");
						ss<< chr_name << "_"<< posposvec[j]<< "_+" << "\t"<< chr_name<< "\t" << posposvec[j] << "\t+\t" << posmstrvec[j] <<  "\t"<< posmeanvec[j][midx] << "\t" <<posvarvec[j][midx] <<"\t"<<posvarvec[j][midx]/(posmeanvec[j][midx]*posmeanvec[j][midx])  ;
						regvec.push_back(ss.str());
					}
				}
				
				//select hetero SNPs
				if(posflagvec[j] ==RepIDVec.size()){
					
					int cov =0 ;
					ss.str("");
					for(int g =0 ; g< posRawRCvec[j].size(); g++){
						ss<< posRawRCvec[j][g];
				
						if(g< posRawRCvec[j].size() -1){
								ss<<"_";
						} 
				
						vector<string> tmpvec = tok.Tokenize( posRawRCvec[j][g] , ",") ;
						for(int s = 0 ; s< tmpvec.size() ; s++){
							vector<string> tmp1vec = tok.Tokenize( tmpvec[s] , ":") ;
							if(tmp1vec.size()>1){
								cov += atoi(tmp1vec[1].c_str());
							}
						}
					}
					string raw_str = ss.str();
					
					
					ss.str("");
					if(cov >=20){
						ss.str("");
			
						ss<< chr_name << "_"<< posposvec[j]<< "_+" << "\t"<< chr_name<< "\t" << posposvec[j] << "\t+\t"  << raw_str << "\t" << cov << "\t" << posMstrvec[j]<< "\t"<< posmstrvec[j]<< "\t"<< posmeanvec[j][Midx] << "\t" <<posmeanvec[j][midx] ;
						testvec.push_back(ss.str());
					}
						
						
				}
				
			}
			
			
			//merge RCs of neg SNPs
			vector<vector<double> > negNRC_Avec,negNRC_Tvec,negNRC_Cvec,negNRC_Gvec ;
			vector<vector<double> > negmeanvec;
			vector<vector<double > > negvarvec;
			vector<vector<string> >negRawRCvec; 
			vector<string > negMstrvec,negmstrvec;
			vector<string > negposvec; 
			vector<int> negflagvec;
		
			for(int j=0 ; j<RepIDVec.size();j++){
				string Rep_prefix = RepIDVec[j];	
				vector<string>	rcvec = tio.dlread_vector( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.rc" )	;	
				cout << "merge:" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.rc" << endl;
				if(j==0){
//				//init the xmean vectors
					for(int q =0 ; q<rcvec.size(); q ++){
						vector<double> tmpAvec; 
						negNRC_Avec.push_back(tmpAvec);
						vector<double> tmpTvec; 
						negNRC_Tvec.push_back(tmpTvec);
						vector<double> tmpCvec; 
						negNRC_Cvec.push_back(tmpCvec);
						vector<double> tmpGvec; 
						negNRC_Gvec.push_back(tmpGvec);
						vector<double> tmpvec;
						negmeanvec.push_back(tmpvec);
						vector<double> tmpvarvec;
						negvarvec.push_back(tmpvarvec);
						
						negflagvec.push_back(0);
						vector<string> qqvec;
						negRawRCvec.push_back(qqvec);
					}	
				}
				for(int k =0 ; k<rcvec.size() ; k++){
					vector<string> tokvec = tok.Tokenize(rcvec[k], "\t" );
					negposvec.push_back(tokvec[0]);
					
					double d_a,d_t,d_c,d_g ;
					d_a = atof(tokvec[7].c_str());
					d_t = atof(tokvec[9].c_str());
					d_c = atof(tokvec[11].c_str());
					d_g = atof(tokvec[13].c_str());
					negNRC_Avec[k].push_back(d_a);
					negNRC_Tvec[k].push_back(d_t);
					negNRC_Cvec[k].push_back(d_c);
					negNRC_Gvec[k].push_back(d_g);
					
					
					//check if a hetero SNP in the rep
					
					string rawstr = tokvec[3] ;
					negRawRCvec[k].push_back(rawstr);
					vector<string> tmpvec = tok.Tokenize(rawstr, "," );
					vector<double> rawrcvec;
					double dep =0;
					for(int q=0 ; q<4; q++){
						vector<string> tmp1vec = tok.Tokenize(tmpvec[q], ":" );
						rawrcvec.push_back(atof(tmp1vec[1].c_str()));
						dep +=atof(tmp1vec[1].c_str() );
					}
					sort (rawrcvec.begin(), rawrcvec.end());
		    	int dep_idx = dep/10;
		    	if(dep_idx>=9){
		    		dep_idx =9;
		    	}
		    	double M_raw = rawrcvec[3];
		    	double m_raw = rawrcvec[2];
					if((m_raw/(m_raw+M_raw)) >= cutoff_vec[dep_idx][j] ){
						negflagvec[k] = negflagvec[k] +1;
						
					}
					
				}
					
			}
			
			//calc mean of ATCG
			
			for(int j =0 ; j< negNRC_Avec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Avec[j].size(); k++){
					sum = sum + negNRC_Avec[j][k] ;
				}
				negmeanvec[j].push_back(sum/(double)negNRC_Avec[j].size());
			}
			
			for(int j =0 ; j< negNRC_Tvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Tvec[j].size(); k++){
					sum = sum + negNRC_Tvec[j][k] ;
				}
				negmeanvec[j].push_back(sum/(double)negNRC_Tvec[j].size());
			}
			
			for(int j =0 ; j< negNRC_Cvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Cvec[j].size(); k++){
					sum = sum + negNRC_Cvec[j][k] ;
				}
				negmeanvec[j].push_back(sum/(double)negNRC_Cvec[j].size());
			}
			
			for(int j =0 ; j< negNRC_Gvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Gvec[j].size(); k++){
					sum = sum + negNRC_Gvec[j][k] ;
				}
				negmeanvec[j].push_back(sum/(double)negNRC_Gvec[j].size());
			}
			
			//calc var of ATCG
			
			
			for(int j =0 ; j< negNRC_Avec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Avec[j].size(); k++){
					sum = sum + (negNRC_Avec[j][k]-negmeanvec[j][0])*(negNRC_Avec[j][k]-negmeanvec[j][0]) ;
				}
				negvarvec[j].push_back(sum/(double)(negNRC_Avec[j].size()-1) );
			}
			
			for(int j =0 ; j< negNRC_Tvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Tvec[j].size(); k++){
					sum = sum + (negNRC_Tvec[j][k]-negmeanvec[j][1])*(negNRC_Tvec[j][k]-negmeanvec[j][1]) ;
				}
				negvarvec[j].push_back(sum/(double)(negNRC_Tvec[j].size()-1) );
			}
			
			for(int j =0 ; j< negNRC_Cvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Cvec[j].size(); k++){
					sum = sum + (negNRC_Cvec[j][k]-negmeanvec[j][2])*(negNRC_Cvec[j][k]-negmeanvec[j][2]) ;
				}
				negvarvec[j].push_back(sum/(double)(negNRC_Cvec[j].size()-1) );
			}
			
			for(int j =0 ; j< negNRC_Gvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Gvec[j].size(); k++){
					sum = sum + (negNRC_Gvec[j][k]-negmeanvec[j][3])*(negNRC_Gvec[j][k]-negmeanvec[j][3]) ;
				}
				negvarvec[j].push_back(sum/(double)(negNRC_Gvec[j].size()-1) );
			}
			
			//determine major and minor allele
			
			for(int j =0 ; j< negmeanvec.size() ; j++){
						
				vector<double> tmpmeanvec = negmeanvec[j];
				sort (tmpmeanvec.begin(), tmpmeanvec.end());
		    double M_rc = tmpmeanvec[3];
		    double m_rc = tmpmeanvec[2];
				
				//search M idx ;
				int Midx =-1;
				int midx =-1;
				for(int k =0 ; k< negmeanvec[j].size() ; k++){
					if(negmeanvec[j][k] == M_rc){
						Midx = k;
					}
					if(negmeanvec[j][k] == m_rc){
						midx = k;
					}
				}
				
				if(Midx ==0 ){
					negMstrvec.push_back("A");
				}else if(Midx ==1 ){
					negMstrvec.push_back("T");
				}else if(Midx ==2 ){
					negMstrvec.push_back("C");
				}else if(Midx ==3 ){
					negMstrvec.push_back("G");
				}else{
					negMstrvec.push_back("-");
				}
				
				
				if(midx ==0 ){
					negmstrvec.push_back("A");
				}else if(midx ==1 ){
					negmstrvec.push_back("T");
				}else if(midx ==2 ){
					negmstrvec.push_back("C");
				}else if(midx ==3 ){
					negmstrvec.push_back("G");
				}else{
					negmstrvec.push_back("-");
				}
				
				//select for reg
				if(Midx > -1){
					if(negmeanvec[j][Midx]>2.0){
						ss.str("");
						ss<< chr_name << "_"<< negposvec[j]<< "_-" << "\t"<< chr_name<< "\t" << negposvec[j] << "\t-\t" << negMstrvec[j] <<  "\t"<< negmeanvec[j][Midx] << "\t" <<negvarvec[j][Midx] <<"\t"<<negvarvec[j][Midx]/(negmeanvec[j][Midx]*negmeanvec[j][Midx])  ;
						regvec.push_back(ss.str());
					}
				}
				if(midx > -1){
					if(negmeanvec[j][midx]>2.0){
						ss.str("");
						ss<< chr_name << "_"<< negposvec[j]<< "_-" << "\t"<< chr_name<< "\t" << negposvec[j] << "\t-\t" << negMstrvec[j] <<  "\t"<< negmeanvec[j][midx] << "\t" <<negvarvec[j][midx] <<"\t"<<negvarvec[j][midx]/(negmeanvec[j][midx]*negmeanvec[j][midx])  ;
						regvec.push_back(ss.str());
					}
				}
				//select hetero SNPs
				if(negflagvec[j] ==RepIDVec.size()){
					
					int cov =0 ;
					ss.str("");
					for(int g =0 ; g< negRawRCvec[j].size(); g++){
						ss<< negRawRCvec[j][g];
				
						if(g< negRawRCvec[j].size() -1){
								ss<<"_";
						} 
				
						vector<string> tmpvec = tok.Tokenize( negRawRCvec[j][g] , ",") ;
						for(int s = 0 ; s< tmpvec.size() ; s++){
							vector<string> tmp1vec = tok.Tokenize( tmpvec[s] , ":") ;
							if(tmp1vec.size()>1){
								cov += atoi(tmp1vec[1].c_str());
							}
						}
					}
					string raw_str = ss.str();
					
					
					ss.str("");
					if(cov >=20){
						ss.str("");
			
						ss<< chr_name << "_"<< negposvec[j]<< "_-" << "\t"<< chr_name<< "\t" << negposvec[j] << "\t-\t"  << raw_str << "\t" << cov << "\t" << negMstrvec[j]<< "\t"<< negmstrvec[j]<< "\t"<< negmeanvec[j][Midx] << "\t" <<negmeanvec[j][midx] ;
						testvec.push_back(ss.str());
					}
						
						
				}
			}
		
	
		}
	
		tio.dlwrite( tmp_dir+"/reg_denovo.rc" , regvec);
    tio.dlwrite( tmp_dir+"/test_denovo.rc" , testvec);
	
	//for testing 
	
	
}

void SNP::SNPtoRC_train(	std::string SNP_dir, 	std::string tmp_dir,  std::vector< std::string > RepIDVec, std::string chr_str, std::map<std::string, double> SizeMap){
	stringstream ss;
	TextIO tio;
	TOK tok;
	vector<string> outVec; 
	
	vector<string> chrVec = tok.Tokenize( chr_str, ",") ; 
	
	//read and normalize raw RC
	
	for(int i =0 ; i<chrVec.size(); i++){
		string chr_name= chrVec[i];
		//read SNP annotation 
		vector<string> SNPVec = tio.dlread_vector(SNP_dir+"/"+chr_name+".txt");
		// ov annotated all SNPs to peak  
	  cout << "read " <<SNP_dir+"/"+chr_name+".txt" <<endl;	  
	  for(int j=0 ; j<RepIDVec.size();j++){
  	
  		string Rep_prefix = RepIDVec[j];  
//			
			listAllOvSNPs_denovo_lite( SNPVec, tmp_dir+"/"+Rep_prefix+"_"+chr_name+"_pos.bed", tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.snp" );			
			cout << "print : "<< tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.snp" << endl;
			listAllOvSNPs_denovo_lite( SNPVec, tmp_dir+"/"+Rep_prefix+"_"+chr_name+"_neg.bed", tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.snp"  );
			cout << "print : "<< tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.snp" << endl;
		}
		
		
		map<string, string> possnpmap,negsnpmap ;
		
//		
	  for(int j=0 ; j<RepIDVec.size();j++){
  	
  		string Rep_prefix = RepIDVec[j];  
			
			vector<string> posvec = tio.dlread_vector(tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.snp");
			vector<string> negvec = tio.dlread_vector(tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.snp");
			
			for(int k =0 ; k<posvec.size() ; k++){
				
				if(possnpmap.find(posvec[k]) == possnpmap.end()){
					possnpmap.insert(make_pair( posvec[k], ""));
				}
			}
			
			for(int k =0 ; k<negvec.size() ; k++){
				
				if(negsnpmap.find(negvec[k]) == negsnpmap.end()){
					negsnpmap.insert(make_pair( negvec[k], ""));
				}
			}
			
		}
////		
		//print the union of SNP 
		vector<string> univec; 
		for(map<string, string>::iterator mitr = possnpmap.begin(); mitr != possnpmap.end(); mitr++){
			univec.push_back(mitr->first);
		}
		
		tio.dlwrite( tmp_dir+"/"+chr_name+"_pos_allovdenovo.snp" , univec);
		cout << "write: "<< tmp_dir+"/"+chr_name+"_pos_allovdenovo.snp"<< endl;
		univec.clear();
		
		for(map<string, string>::iterator mitr = negsnpmap.begin(); mitr != negsnpmap.end(); mitr++){
			univec.push_back(mitr->first);
		}
		
		tio.dlwrite( tmp_dir+"/"+chr_name+"_neg_allovdenovo.snp" , univec);
		cout << "write: "<< tmp_dir+"/"+chr_name+"_neg_allovdenovo.snp"<< endl;
		univec.clear();
		
		
		
		//read read cov
		for(int j=0 ; j<RepIDVec.size();j++){
			string Rep_prefix = RepIDVec[j];
			vector<string> SNPvec = tio.dlread_vector(tmp_dir+"/"+chr_name+"_pos_allovdenovo.snp");
			readRawRC_denovo(SNPvec, tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+".mp" , tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawdenovo.rc");
			cout << "calc raw RC :" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawdenovo.rc" << endl;
		
		//nomalize raw read counts			
			vector<string>  SNP_RepVec= normalize_denovo( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawdenovo.rc",  tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls", tmp_dir , 1 , SizeMap[Rep_prefix]);
			tio.dlwrite( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.rc" , SNP_RepVec);
			cout << "normalize raw RC :" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.rc\t"<< tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls" << "\t"<<SizeMap[Rep_prefix] << endl;
////			
			SNPvec = tio.dlread_vector(tmp_dir+"/"+chr_name+"_neg_allovdenovo.snp");
			readRawRC_denovo(SNPvec, tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+".mp" , tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawdenovo.rc");
			cout << "calc raw RC :" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawdenovo.rc" << endl;
////			
			SNP_RepVec= normalize_denovo( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawdenovo.rc",  tmp_dir+"/"+Rep_prefix+"_neg_2_"+chr_name+".cls", tmp_dir , 0 , SizeMap[Rep_prefix]);
			tio.dlwrite( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.rc" , SNP_RepVec);
			cout << "normalize raw RC :" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.rc\t"<< tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls" << "\t"<<SizeMap[Rep_prefix] << endl;
		}
		
	}
		
	//for regression training 
	
	
		vector<string> regvec ; 
		regvec.push_back("id\tchr\tpos\tstrand\tallele\txmean\tvar\tcv2");
		vector<string> testvec ; 
	  testvec.push_back("id\tchr\tpos\tstrand\traw_str\tcov\tM_str\tm_str\tM\tm");
	   
		for(int i =0 ; i<chrVec.size(); i++){
			
			string chr_name = chrVec[i];
			//merge RCs of pos SNPs
			
			vector<vector<double> > posNRC_Avec,posNRC_Tvec,posNRC_Cvec,posNRC_Gvec ;
			vector<vector<double> > posmeanvec;
			vector<vector<double > > posvarvec;
			vector<vector<string > > posRawRCvec;
			vector<string > posMstrvec,posmstrvec;
			
			vector<string > posposvec; 
			vector<int> posflagvec;
			
			for(int j=0 ; j<RepIDVec.size();j++){
				string Rep_prefix = RepIDVec[j];	
				vector<string>	rcvec = tio.dlread_vector( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.rc" )	;	
				cout << "merge:" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_denovo.rc" << endl;
				if(j==0){
//				//init the xmean vectors
					for(int q =0 ; q<rcvec.size(); q ++){
						vector<double> tmpAvec; 
						posNRC_Avec.push_back(tmpAvec);
						vector<double> tmpTvec; 
						posNRC_Tvec.push_back(tmpTvec);
						vector<double> tmpCvec; 
						posNRC_Cvec.push_back(tmpCvec);
						vector<double> tmpGvec; 
						posNRC_Gvec.push_back(tmpGvec);
						vector<double> tmpvec;
						posmeanvec.push_back(tmpvec);
						vector<double> tmpvarvec;
						posvarvec.push_back(tmpvarvec);
						
						posflagvec.push_back(0);
						
						vector<string> qqvec;
						posRawRCvec.push_back(qqvec);
						
					}	
				}
				for(int k =0 ; k<rcvec.size() ; k++){
					vector<string> tokvec = tok.Tokenize(rcvec[k], "\t" );
					posposvec.push_back(tokvec[0]);
					
					double d_a,d_t,d_c,d_g ;
					d_a = atof(tokvec[7].c_str());
					d_t = atof(tokvec[9].c_str());
					d_c = atof(tokvec[11].c_str());
					d_g = atof(tokvec[13].c_str());
					posNRC_Avec[k].push_back(d_a);
					posNRC_Tvec[k].push_back(d_t);
					posNRC_Cvec[k].push_back(d_c);
					posNRC_Gvec[k].push_back(d_g);
					
					//check if a hetero SNP in the rep
					
					string rawstr = tokvec[3] ;
					posRawRCvec[k].push_back(rawstr);
					vector<string> tmpvec = tok.Tokenize(rawstr, "," );
					vector<double> rawrcvec;
					double dep =0;
					for(int q=0 ; q<4; q++){
						vector<string> tmp1vec = tok.Tokenize(tmpvec[q], ":" );
						rawrcvec.push_back(atof(tmp1vec[1].c_str()));
						dep +=atof(tmp1vec[1].c_str() );
					}
					sort (rawrcvec.begin(), rawrcvec.end());
		    	int dep_idx = dep/10;
		    	if(dep_idx>=9){
		    		dep_idx =9;
		    	}
		    	double M_raw = rawrcvec[3];
		    	double m_raw = rawrcvec[2];
//					if((m_raw/(m_raw+M_raw)) >= cutoff_vec[dep_idx][j] ){
//						posflagvec[k] = posflagvec[k] +1;
//						
//					}
					
				}
				
				
				
				
					
			}
			
			//calc mean of ATCG
			
			for(int j =0 ; j< posNRC_Avec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Avec[j].size(); k++){
					sum = sum + posNRC_Avec[j][k] ;
				}
				posmeanvec[j].push_back(sum/(double)posNRC_Avec[j].size());
			}
			
			for(int j =0 ; j< posNRC_Tvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Tvec[j].size(); k++){
					sum = sum + posNRC_Tvec[j][k] ;
				}
				posmeanvec[j].push_back(sum/(double)posNRC_Tvec[j].size());
			}
			
			for(int j =0 ; j< posNRC_Cvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Cvec[j].size(); k++){
					sum = sum + posNRC_Cvec[j][k] ;
				}
				posmeanvec[j].push_back(sum/(double)posNRC_Cvec[j].size());
			}
			
			for(int j =0 ; j< posNRC_Gvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Gvec[j].size(); k++){
					sum = sum + posNRC_Gvec[j][k] ;
				}
				posmeanvec[j].push_back(sum/(double)posNRC_Gvec[j].size());
			}
			
			//calc var of ATCG
			
			
			for(int j =0 ; j< posNRC_Avec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Avec[j].size(); k++){
					sum = sum + (posNRC_Avec[j][k]-posmeanvec[j][0])*(posNRC_Avec[j][k]-posmeanvec[j][0]) ;
				}
				posvarvec[j].push_back(sum/(double)(posNRC_Avec[j].size()-1) );
			}
			
			for(int j =0 ; j< posNRC_Tvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Tvec[j].size(); k++){
					sum = sum + (posNRC_Tvec[j][k]-posmeanvec[j][1])*(posNRC_Tvec[j][k]-posmeanvec[j][1]) ;
				}
				posvarvec[j].push_back(sum/(double)(posNRC_Tvec[j].size()-1) );
			}
			
			for(int j =0 ; j< posNRC_Cvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Cvec[j].size(); k++){
					sum = sum + (posNRC_Cvec[j][k]-posmeanvec[j][2])*(posNRC_Cvec[j][k]-posmeanvec[j][2]) ;
				}
				posvarvec[j].push_back(sum/(double)(posNRC_Cvec[j].size()-1) );
			}
			
			for(int j =0 ; j< posNRC_Gvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< posNRC_Gvec[j].size(); k++){
					sum = sum + (posNRC_Gvec[j][k]-posmeanvec[j][3])*(posNRC_Gvec[j][k]-posmeanvec[j][3]) ;
				}
				posvarvec[j].push_back(sum/(double)(posNRC_Gvec[j].size()-1) );
			}
			
			//determine major and minor allele
			
			for(int j =0 ; j< posmeanvec.size() ; j++){
						
				vector<double> tmpmeanvec = posmeanvec[j];
				sort (tmpmeanvec.begin(), tmpmeanvec.end());
		    double M_rc = tmpmeanvec[3];
		    double m_rc = tmpmeanvec[2];
				
				//search M idx ;
				int Midx =-1;
				int midx =-1;
				for(int k =0 ; k< posmeanvec[j].size() ; k++){
					if(posmeanvec[j][k] == M_rc){
						Midx = k;
					}
					if(posmeanvec[j][k] == m_rc){
						midx = k;
					}
				}
				
				if(Midx ==0 ){
					posMstrvec.push_back("A");
				}else if(Midx ==1 ){
					posMstrvec.push_back("T");
				}else if(Midx ==2 ){
					posMstrvec.push_back("C");
				}else if(Midx ==3 ){
					posMstrvec.push_back("G");
				}else{
					posMstrvec.push_back("-");
				}
				
				
				if(midx ==0 ){
					posmstrvec.push_back("A");
				}else if(midx ==1 ){
					posmstrvec.push_back("T");
				}else if(midx ==2 ){
					posmstrvec.push_back("C");
				}else if(midx ==3 ){
					posmstrvec.push_back("G");
				}else{
					posmstrvec.push_back("-");
				}
				
				//select for reg
				if(Midx > -1){
					if(posmeanvec[j][Midx]>2.0){
						ss.str("");
						ss<< chr_name << "_"<< posposvec[j]<< "_+" << "\t"<< chr_name<< "\t" << posposvec[j] << "\t+\t" << posMstrvec[j] <<  "\t"<< posmeanvec[j][Midx] << "\t" <<posvarvec[j][Midx] <<"\t"<<posvarvec[j][Midx]/(posmeanvec[j][Midx]*posmeanvec[j][Midx])  ;
						regvec.push_back(ss.str());
					}
				}
				if(midx > -1){
					if(posmeanvec[j][midx]>2.0){
						ss.str("");
						ss<< chr_name << "_"<< posposvec[j]<< "_+" << "\t"<< chr_name<< "\t" << posposvec[j] << "\t+\t" << posmstrvec[j] <<  "\t"<< posmeanvec[j][midx] << "\t" <<posvarvec[j][midx] <<"\t"<<posvarvec[j][midx]/(posmeanvec[j][midx]*posmeanvec[j][midx])  ;
						regvec.push_back(ss.str());
					}
				}
				
				//select hetero SNPs
//				if(posflagvec[j] ==RepIDVec.size()){
//					
//					int cov =0 ;
//					ss.str("");
//					for(int g =0 ; g< posRawRCvec[j].size(); g++){
//						ss<< posRawRCvec[j][g];
//				
//						if(g< posRawRCvec[j].size() -1){
//								ss<<"_";
//						} 
//				
//						vector<string> tmpvec = tok.Tokenize( posRawRCvec[j][g] , ",") ;
//						for(int s = 0 ; s< tmpvec.size() ; s++){
//							vector<string> tmp1vec = tok.Tokenize( tmpvec[s] , ":") ;
//							if(tmp1vec.size()>1){
//								cov += atoi(tmp1vec[1].c_str());
//							}
//						}
//					}
//					string raw_str = ss.str();
//					
//					
//					ss.str("");
//					if(cov >=20){
//						ss.str("");
//			
//						ss<< chr_name << "_"<< posposvec[j]<< "_+" << "\t"<< chr_name<< "\t" << posposvec[j] << "\t+\t"  << raw_str << "\t" << cov << "\t" << posMstrvec[j]<< "\t"<< posmstrvec[j]<< "\t"<< posmeanvec[j][Midx] << "\t" <<posmeanvec[j][midx] ;
//						testvec.push_back(ss.str());
//					}
//						
//						
//				}
				
			}
			
			
			//merge RCs of neg SNPs
			vector<vector<double> > negNRC_Avec,negNRC_Tvec,negNRC_Cvec,negNRC_Gvec ;
			vector<vector<double> > negmeanvec;
			vector<vector<double > > negvarvec;
			vector<vector<string> >negRawRCvec; 
			vector<string > negMstrvec,negmstrvec;
			vector<string > negposvec; 
		
			for(int j=0 ; j<RepIDVec.size();j++){
				string Rep_prefix = RepIDVec[j];	
				vector<string>	rcvec = tio.dlread_vector( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.rc" )	;	
				cout << "merge:" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_denovo.rc" << endl;
				if(j==0){
//				//init the xmean vectors
					for(int q =0 ; q<rcvec.size(); q ++){
						vector<double> tmpAvec; 
						negNRC_Avec.push_back(tmpAvec);
						vector<double> tmpTvec; 
						negNRC_Tvec.push_back(tmpTvec);
						vector<double> tmpCvec; 
						negNRC_Cvec.push_back(tmpCvec);
						vector<double> tmpGvec; 
						negNRC_Gvec.push_back(tmpGvec);
						vector<double> tmpvec;
						negmeanvec.push_back(tmpvec);
						vector<double> tmpvarvec;
						negvarvec.push_back(tmpvarvec);
						
						vector<string> qqvec;
						negRawRCvec.push_back(qqvec);
					}	
				}
				for(int k =0 ; k<rcvec.size() ; k++){
					vector<string> tokvec = tok.Tokenize(rcvec[k], "\t" );
					negposvec.push_back(tokvec[0]);
					
					double d_a,d_t,d_c,d_g ;
					d_a = atof(tokvec[7].c_str());
					d_t = atof(tokvec[9].c_str());
					d_c = atof(tokvec[11].c_str());
					d_g = atof(tokvec[13].c_str());
					negNRC_Avec[k].push_back(d_a);
					negNRC_Tvec[k].push_back(d_t);
					negNRC_Cvec[k].push_back(d_c);
					negNRC_Gvec[k].push_back(d_g);
					
					
					//check if a hetero SNP in the rep
					
					string rawstr = tokvec[3] ;
					negRawRCvec[k].push_back(rawstr);
					vector<string> tmpvec = tok.Tokenize(rawstr, "," );
					vector<double> rawrcvec;
					double dep =0;
					for(int q=0 ; q<4; q++){
						vector<string> tmp1vec = tok.Tokenize(tmpvec[q], ":" );
						rawrcvec.push_back(atof(tmp1vec[1].c_str()));
						dep +=atof(tmp1vec[1].c_str() );
					}
					sort (rawrcvec.begin(), rawrcvec.end());
		    	int dep_idx = dep/10;
		    	if(dep_idx>=9){
		    		dep_idx =9;
		    	}
		    	double M_raw = rawrcvec[3];
		    	double m_raw = rawrcvec[2];
//					if((m_raw/(m_raw+M_raw)) >= cutoff_vec[dep_idx][j] ){
//						negflagvec[k] = negflagvec[k] +1;
//						
//					}
					
				}
					
			}
			
			//calc mean of ATCG
			
			for(int j =0 ; j< negNRC_Avec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Avec[j].size(); k++){
					sum = sum + negNRC_Avec[j][k] ;
				}
				negmeanvec[j].push_back(sum/(double)negNRC_Avec[j].size());
			}
			
			for(int j =0 ; j< negNRC_Tvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Tvec[j].size(); k++){
					sum = sum + negNRC_Tvec[j][k] ;
				}
				negmeanvec[j].push_back(sum/(double)negNRC_Tvec[j].size());
			}
			
			for(int j =0 ; j< negNRC_Cvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Cvec[j].size(); k++){
					sum = sum + negNRC_Cvec[j][k] ;
				}
				negmeanvec[j].push_back(sum/(double)negNRC_Cvec[j].size());
			}
			
			for(int j =0 ; j< negNRC_Gvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Gvec[j].size(); k++){
					sum = sum + negNRC_Gvec[j][k] ;
				}
				negmeanvec[j].push_back(sum/(double)negNRC_Gvec[j].size());
			}
			
			//calc var of ATCG
			
			
			for(int j =0 ; j< negNRC_Avec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Avec[j].size(); k++){
					sum = sum + (negNRC_Avec[j][k]-negmeanvec[j][0])*(negNRC_Avec[j][k]-negmeanvec[j][0]) ;
				}
				negvarvec[j].push_back(sum/(double)(negNRC_Avec[j].size()-1) );
			}
			
			for(int j =0 ; j< negNRC_Tvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Tvec[j].size(); k++){
					sum = sum + (negNRC_Tvec[j][k]-negmeanvec[j][1])*(negNRC_Tvec[j][k]-negmeanvec[j][1]) ;
				}
				negvarvec[j].push_back(sum/(double)(negNRC_Tvec[j].size()-1) );
			}
			
			for(int j =0 ; j< negNRC_Cvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Cvec[j].size(); k++){
					sum = sum + (negNRC_Cvec[j][k]-negmeanvec[j][2])*(negNRC_Cvec[j][k]-negmeanvec[j][2]) ;
				}
				negvarvec[j].push_back(sum/(double)(negNRC_Cvec[j].size()-1) );
			}
			
			for(int j =0 ; j< negNRC_Gvec.size() ; j++){
				double sum =0 ;
				for(int k =0 ; k< negNRC_Gvec[j].size(); k++){
					sum = sum + (negNRC_Gvec[j][k]-negmeanvec[j][3])*(negNRC_Gvec[j][k]-negmeanvec[j][3]) ;
				}
				negvarvec[j].push_back(sum/(double)(negNRC_Gvec[j].size()-1) );
			}
			
			//determine major and minor allele
			
			for(int j =0 ; j< negmeanvec.size() ; j++){
						
				vector<double> tmpmeanvec = negmeanvec[j];
				sort (tmpmeanvec.begin(), tmpmeanvec.end());
		    double M_rc = tmpmeanvec[3];
		    double m_rc = tmpmeanvec[2];
				
				//search M idx ;
				int Midx =-1;
				int midx =-1;
				for(int k =0 ; k< negmeanvec[j].size() ; k++){
					if(negmeanvec[j][k] == M_rc){
						Midx = k;
					}
					if(negmeanvec[j][k] == m_rc){
						midx = k;
					}
				}
				
				if(Midx ==0 ){
					negMstrvec.push_back("A");
				}else if(Midx ==1 ){
					negMstrvec.push_back("T");
				}else if(Midx ==2 ){
					negMstrvec.push_back("C");
				}else if(Midx ==3 ){
					negMstrvec.push_back("G");
				}else{
					negMstrvec.push_back("-");
				}
				
				
				if(midx ==0 ){
					negmstrvec.push_back("A");
				}else if(midx ==1 ){
					negmstrvec.push_back("T");
				}else if(midx ==2 ){
					negmstrvec.push_back("C");
				}else if(midx ==3 ){
					negmstrvec.push_back("G");
				}else{
					negmstrvec.push_back("-");
				}
				
				//select for reg
				if(Midx > -1){
					if(negmeanvec[j][Midx]>2.0){
						ss.str("");
						ss<< chr_name << "_"<< negposvec[j]<< "_-" << "\t"<< chr_name<< "\t" << negposvec[j] << "\t-\t" << negMstrvec[j] <<  "\t"<< negmeanvec[j][Midx] << "\t" <<negvarvec[j][Midx] <<"\t"<<negvarvec[j][Midx]/(negmeanvec[j][Midx]*negmeanvec[j][Midx])  ;
						regvec.push_back(ss.str());
					}
				}
				if(midx > -1){
					if(negmeanvec[j][midx]>2.0){
						ss.str("");
						ss<< chr_name << "_"<< negposvec[j]<< "_-" << "\t"<< chr_name<< "\t" << negposvec[j] << "\t-\t" << negMstrvec[j] <<  "\t"<< negmeanvec[j][midx] << "\t" <<negvarvec[j][midx] <<"\t"<<negvarvec[j][midx]/(negmeanvec[j][midx]*negmeanvec[j][midx])  ;
						regvec.push_back(ss.str());
					}
				}
				//select hetero SNPs
//				if(negflagvec[j] ==RepIDVec.size()){
//					
//					int cov =0 ;
//					ss.str("");
//					for(int g =0 ; g< negRawRCvec[j].size(); g++){
//						ss<< negRawRCvec[j][g];
//				
//						if(g< negRawRCvec[j].size() -1){
//								ss<<"_";
//						} 
//				
//						vector<string> tmpvec = tok.Tokenize( negRawRCvec[j][g] , ",") ;
//						for(int s = 0 ; s< tmpvec.size() ; s++){
//							vector<string> tmp1vec = tok.Tokenize( tmpvec[s] , ":") ;
//							if(tmp1vec.size()>1){
//								cov += atoi(tmp1vec[1].c_str());
//							}
//						}
//					}
//					string raw_str = ss.str();
//					
//					
//					ss.str("");
//					if(cov >=20){
//						ss.str("");
//			
//						ss<< chr_name << "_"<< negposvec[j]<< "_-" << "\t"<< chr_name<< "\t" << negposvec[j] << "\t-\t"  << raw_str << "\t" << cov << "\t" << negMstrvec[j]<< "\t"<< negmstrvec[j]<< "\t"<< negmeanvec[j][Midx] << "\t" <<negmeanvec[j][midx] ;
//						testvec.push_back(ss.str());
//					}
//						
//						
//				}
			}
		
	
		}
	
		tio.dlwrite( tmp_dir+"/reg_denovo.rc" , regvec);
    //tio.dlwrite( tmp_dir+"/test_denovo.rc" , testvec);
	
	//for testing 
	
	
}


void SNP::SNPtoRC_VCF(	std::string VCFfname, 	std::string tmp_dir, std::vector< std::string > RepIDVec , std::string chr_str, std::map<std::string, double> SizeMap){
	stringstream ss;
	TextIO tio;
	TOK tok;
	vector<string> outVec; 
	map<string, vector<string> > vcfmap , heteromap ;
	vector<string> chrVec = tok.Tokenize( chr_str, ",") ; 
	
	//read the VCF file and store by chrms
	
	cout << "reading biallelic SNVs from the VCF file :" << VCFfname << endl; 
	vector<string> vcfvec = tio.dlread_vector(VCFfname);
	for(int i =0 ; i < vcfvec.size() ; i++){
		if(vcfvec[i][0] != '#'){
					vector<string> tokvec = tok.Tokenize( vcfvec[i] , "\t");
					string chr_name = tokvec[0];
					vector<string> tmpvec = tok.Tokenize( tokvec[7] , ";");
					string af_str = tmpvec[1].substr(3); 
					vector<string> afvec = tok.Tokenize( af_str , ",");
					af_str = afvec[0];
					string snp_id = tokvec[0]+"_"+tokvec[1];
					double af = atof(af_str.c_str());
					if(tokvec[3].length()==1 && tokvec[4].length()==1){ // check if biallelic
						
						if(vcfmap.find(chr_name) == vcfmap.end()){
								vector<string> qqvec ;
								qqvec.push_back(vcfvec[i]);
								vcfmap.insert(make_pair( chr_name , qqvec ));
						}
						else{
								vcfmap[chr_name].push_back(vcfvec[i]);
						}
						
						
						
						if(af!=1){
						 
							if(heteromap.find(chr_name) == heteromap.end()){
								vector<string> qqvec ;
								qqvec.push_back(vcfvec[i]);
								heteromap.insert(make_pair( chr_name , qqvec ));
							}
							else{
								heteromap[chr_name].push_back(vcfvec[i]);
							}
						}	
					}
					
		}
	} 

//for regression training 
	
	vector<string> regvec ; 
	regvec.push_back("id\tchr\tpos\tstrand\tallele\txmean\tvar\tcv2");
	for(int i =0 ; i<chrVec.size(); i++){
		string chr_name= chrVec[i];
		
		// ov annotated all SNPs to peak  
	  	  
	  for(int j=0 ; j<RepIDVec.size();j++){
  	
  		string Rep_prefix = RepIDVec[j];  
			
			listAllOvSNPs_VCF_lite( vcfmap[chr_name], tmp_dir+"/"+Rep_prefix+"_"+chr_name+"_pos.bed", tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_vcf.snp" );			cout << "print : "<< tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_vcf.snp" << endl;
			listAllOvSNPs_VCF_lite( vcfmap[chr_name], tmp_dir+"/"+Rep_prefix+"_"+chr_name+"_neg.bed", tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_vcf.snp"  );
			cout << "print : "<< tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_vcf.snp" << endl;
		}
		
		//merge SNPs in rep
		
		map<string, string> possnpmap,negsnpmap ;
		
//		
	  for(int j=0 ; j<RepIDVec.size();j++){
  	
  		string Rep_prefix = RepIDVec[j];  
			
			vector<string> posvec = tio.dlread_vector(tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_vcf.snp");
			vector<string> negvec = tio.dlread_vector(tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_vcf.snp");
			
			for(int k =0 ; k<posvec.size() ; k++){
				
				if(possnpmap.find(posvec[k]) == possnpmap.end()){
					possnpmap.insert(make_pair( posvec[k], ""));
				}
			}
			
			for(int k =0 ; k<negvec.size() ; k++){
				
				if(negsnpmap.find(negvec[k]) == negsnpmap.end()){
					negsnpmap.insert(make_pair( negvec[k], ""));
				}
			}
			
		}
////		
		//print the union of SNP 
		vector<string> univec; 
		for(map<string, string>::iterator mitr = possnpmap.begin(); mitr != possnpmap.end(); mitr++){
			univec.push_back(mitr->first);
		}
		
		tio.dlwrite( tmp_dir+"/"+chr_name+"_pos_allov.snp" , univec);
		cout << "write: "<< tmp_dir+"/"+chr_name+"_pos_allov.snp"<< endl;
		univec.clear();
		
		for(map<string, string>::iterator mitr = negsnpmap.begin(); mitr != negsnpmap.end(); mitr++){
			univec.push_back(mitr->first);
		}
		
		tio.dlwrite( tmp_dir+"/"+chr_name+"_neg_allov.snp" , univec);
		cout << "write: "<< tmp_dir+"/"+chr_name+"_neg_allov.snp"<< endl;
		univec.clear();
//		
		
		
		//read read coverage in reps
		 
		for(int j=0 ; j<RepIDVec.size();j++){
			string Rep_prefix = RepIDVec[j];
			vector<string> SNPvec = tio.dlread_vector(tmp_dir+"/"+chr_name+"_pos_allov.snp");
			readRawRC(SNPvec, tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+".mp" , tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawall.rc");
			cout << "calc raw RC :" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawall.rc" << endl;
			
			vector<string>  SNP_RepVec= normalize_VCF( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawall.rc",  tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls", tmp_dir , 1 , SizeMap[Rep_prefix]);
			tio.dlwrite( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_vcf.rc" , SNP_RepVec);
			cout << "normalize raw RC :" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_vcf.rc\t"<< tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls" << "\t"<<SizeMap[Rep_prefix] << endl;
			
			SNPvec = tio.dlread_vector(tmp_dir+"/"+chr_name+"_neg_allov.snp");
			readRawRC(SNPvec, tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+".mp" , tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawall.rc");
			cout << "calc raw RC :" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawall.rc" << endl;
			
			SNP_RepVec= normalize_VCF( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawall.rc",  tmp_dir+"/"+Rep_prefix+"_neg_2_"+chr_name+".cls", tmp_dir , 0 , SizeMap[Rep_prefix]);
			tio.dlwrite( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_vcf.rc" , SNP_RepVec);
			cout << "normalize raw RC :" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_vcf.rc\t"<< tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls" << "\t"<<SizeMap[Rep_prefix] << endl;
		}
		
		//merge RCs of pos SNPs
		vector<vector<double> > posNRCvec ;
		vector<double > posmeanvec;
		vector<double > posvarvec;
		vector<string > posallelevec;
		vector<string > posposvec; 
		
		for(int j=0 ; j<RepIDVec.size();j++){
			string Rep_prefix = RepIDVec[j];	
			vector<string>	rcvec = tio.dlread_vector( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_vcf.rc" )	;	
			cout << "merge:" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_vcf.rc" << endl;

			if(j==0){
				//init the xmean vectors
				for(int q =0 ; q<2*rcvec.size(); q ++){
					vector<double> tmpvec; 
					posNRCvec.push_back(tmpvec);
				}
			}
			
			for(int k =0 ; k<rcvec.size() ; k++){
				vector<string> tokvec = tok.Tokenize(rcvec[k], "\t" );
				posposvec.push_back(tokvec[0]);
				posposvec.push_back(tokvec[0]);
				posallelevec.push_back(tokvec[1]);
				posallelevec.push_back(tokvec[2]);
				double M_rc, m_rc ;
				M_rc = atof(tokvec[14].c_str());
				m_rc = atof(tokvec[15].c_str());
				posNRCvec[2*k].push_back(M_rc); 
				posNRCvec[2*k+1].push_back(m_rc);
			}
			

		}	
		
		//calc xmean 
		
		for(int q=0 ; q<posNRCvec.size(); q++) {
			double sum =0;
			for(int t = 0 ; t< posNRCvec[q].size(); t++){
				
				sum += posNRCvec[q][t];
				
			}
			posmeanvec.push_back(sum/(double)posNRCvec[q].size() );
			//cout << endl;
		}
		
		//calc var 
		
		for(int q=0 ; q<posNRCvec.size(); q++) {
			double sum =0;
			for(int t = 0 ; t< posNRCvec[q].size(); t++){
				
				sum += (posNRCvec[q][t]-posmeanvec[q])*(posNRCvec[q][t]-posmeanvec[q]);
				
			}
			posvarvec.push_back(sum/(double)(posNRCvec[q].size()-1) );
			//cout << endl;
		}
		
		//select allele for regression
		
		for(int q =0 ; q< posmeanvec.size(); q++){
			if(posmeanvec[q] > 2.0){
				ss.str("");
				ss<< chr_name << "_"<< posposvec[q]<< "_+" << "\t"<< chr_name<< "\t" << posposvec[q] << "\t+\t" << posallelevec[q]<<  "\t"<< posmeanvec[q] << "\t" <<posvarvec[q] <<"\t"<<posvarvec[q]/(posmeanvec[q]*posmeanvec[q])  ;
				regvec.push_back(ss.str());
			}
		}

// merge RCs of neg SNPs
		vector<vector<double> > negNRCvec ;
		vector<double > negmeanvec;
		vector<double > negvarvec;
		vector<string > negallelevec;
		vector<string > negposvec; 
		
		for(int j=0 ; j<RepIDVec.size();j++){
			string Rep_prefix = RepIDVec[j];	
			vector<string>	rcvec = tio.dlread_vector( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_vcf.rc" )	;	
			cout << "merge:" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_vcf.rc" << endl;
			if(j==0){
				//init the xmean vectors
				for(int q =0 ; q<2*rcvec.size(); q ++){
					vector<double> tmpvec; 
					negNRCvec.push_back(tmpvec);
				}
			}
			
			for(int k =0 ; k<rcvec.size() ; k++){
				vector<string> tokvec = tok.Tokenize(rcvec[k], "\t" );
				//cout << "rcvec[k]:"<<  rcvec[k] << endl;
				negposvec.push_back(tokvec[0]);
				negposvec.push_back(tokvec[0]);
				negallelevec.push_back(tokvec[1]);
				negallelevec.push_back(tokvec[2]);
				double M_rc, m_rc ;
				M_rc = atof(tokvec[14].c_str());
				m_rc = atof(tokvec[15].c_str());
				
				//cout << "M_rc :"<<  M_rc << endl;
				//cout << "m_rc :"<<  m_rc << endl;
				
				negNRCvec[2*k].push_back(M_rc); 
				negNRCvec[2*k+1].push_back(m_rc);
				//cout << "read :"<<  negNRCvec[2*k][j] << endl;
				//cout << "read :"<<  negNRCvec[2*k+1][j] << endl;
			}
			

		}	
		
	//		
		//calc xmean 
		
		for(int q=0 ; q<negNRCvec.size(); q++) {
			double sum =0;
			for(int t = 0 ; t< negNRCvec[q].size(); t++){
				
				sum += negNRCvec[q][t];
				
			}
			negmeanvec.push_back(sum/(double)negNRCvec[q].size() );
			//cout << endl;
		}
		
		//calc var 
		
		for(int q=0 ; q<negNRCvec.size(); q++) {
			double sum =0;
			for(int t = 0 ; t< negNRCvec[q].size(); t++){
				
				sum += (negNRCvec[q][t]-negmeanvec[q])*(negNRCvec[q][t]-negmeanvec[q]);
				
			}
			negvarvec.push_back(sum/(double)(negNRCvec[q].size()-1) );
			//cout << endl;
		}
		
		//select allele for regression
		
		for(int q =0 ; q< negmeanvec.size(); q++){
			if(negmeanvec[q] > 2.0){
				//cout << q << ":"<< negNRCvec[q][0] << ","<<negNRCvec[q][1] << endl;
				ss.str("");
				ss<< chr_name << "_"<< negposvec[q]<< "_-" << "\t"<< chr_name<< "\t" << negposvec[q] << "\t-\t" << negallelevec[q]<<  "\t"<< negmeanvec[q] << "\t" <<negvarvec[q] <<"\t"<<negvarvec[q]/(negmeanvec[q]*negmeanvec[q])  ;
				regvec.push_back(ss.str());
			}
		}
			
	}
	
	tio.dlwrite( tmp_dir+"/reg_vcf.rc" , regvec);
	
	
	//for hetero. testing
	 
	cout << "for hetero. testing..... "<< endl; 	
	vector<string> testvec ; 
	testvec.push_back("id\tchr\tpos\tstrand\traw_str\tcov\tM_str\tm_str\tM\tm");
	for(int i =0 ; i<chrVec.size(); i++){
		string chr_name= chrVec[i];
		
		// ov annotated all SNPs to peak  
	  	  
	  for(int j=0 ; j<RepIDVec.size();j++){
  	
  		string Rep_prefix = RepIDVec[j];  
			
			listAllOvSNPs_VCF_lite( heteromap[chr_name], tmp_dir+"/"+Rep_prefix+"_"+chr_name+"_pos.bed", tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.snp" );			
			cout << "print : "<< tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.snp" << endl;
			listAllOvSNPs_VCF_lite( heteromap[chr_name], tmp_dir+"/"+Rep_prefix+"_"+chr_name+"_neg.bed", tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_het.snp"  );
			cout << "print : "<< tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_het.snp" << endl;
		}
		
//		//merge SNPs in rep
//		
		map<string, string> possnpmap,negsnpmap ;
		
//		
	  for(int j=0 ; j<RepIDVec.size();j++){
  	
  		string Rep_prefix = RepIDVec[j];  
			
			vector<string> posvec = tio.dlread_vector(tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.snp");
			vector<string> negvec = tio.dlread_vector(tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_het.snp");
			
			for(int k =0 ; k<posvec.size() ; k++){
				
				if(possnpmap.find(posvec[k]) == possnpmap.end()){
					possnpmap.insert(make_pair( posvec[k], ""));
				}
			}
			
			for(int k =0 ; k<negvec.size() ; k++){
				
				if(negsnpmap.find(negvec[k]) == negsnpmap.end()){
					negsnpmap.insert(make_pair( negvec[k], ""));
				}
			}
			
		}
////		
		//print the union of SNP 
		vector<string> univec; 
		for(map<string, string>::iterator mitr = possnpmap.begin(); mitr != possnpmap.end(); mitr++){
			univec.push_back(mitr->first);
		}
		
		tio.dlwrite( tmp_dir+"/"+chr_name+"_pos_hetov.snp" , univec);
		cout << "write: "<< tmp_dir+"/"+chr_name+"_pos_hetov.snp"<< endl;
		univec.clear();
		
		for(map<string, string>::iterator mitr = negsnpmap.begin(); mitr != negsnpmap.end(); mitr++){
			univec.push_back(mitr->first);
		}
		
		tio.dlwrite( tmp_dir+"/"+chr_name+"_neg_hetov.snp" , univec);
		cout << "write: "<< tmp_dir+"/"+chr_name+"_neg_hetov.snp"<< endl;
		univec.clear();
//		
//		
//		
//		//read read coverage in reps
		 
		for(int j=0 ; j<RepIDVec.size();j++){
			string Rep_prefix = RepIDVec[j];
			vector<string> SNPvec = tio.dlread_vector(tmp_dir+"/"+chr_name+"_pos_hetov.snp");
			readRawRC(SNPvec, tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+".mp" , tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawhet.rc");
			cout << "calc raw RC :" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawhet.rc" << endl;
			
			vector<string>  SNP_RepVec= normalize_VCF( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawhet.rc",  tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls", tmp_dir , 1 , SizeMap[Rep_prefix]);
			tio.dlwrite( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.rc" , SNP_RepVec);
			cout << "normalize raw RC :" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.rc\t"<< tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls" << "\t"<<SizeMap[Rep_prefix] << endl;
			
			SNPvec = tio.dlread_vector(tmp_dir+"/"+chr_name+"_neg_hetov.snp");
			readRawRC(SNPvec, tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+".mp" , tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawhet.rc");
			cout << "calc raw RC :" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawhet.rc" << endl;
			
			SNP_RepVec= normalize_VCF( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawhet.rc",  tmp_dir+"/"+Rep_prefix+"_neg_2_"+chr_name+".cls", tmp_dir , 0 , SizeMap[Rep_prefix]);
			tio.dlwrite( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_het.rc" , SNP_RepVec);
			cout << "normalize raw RC :" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_het.rc\t"<< tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls" << "\t"<<SizeMap[Rep_prefix] << endl;
		}
//		
//		//merge RCs of pos SNPs
		vector<vector<double> > posNRCvec ;
		vector<vector<string> > posrawvec ;
		vector<double > posmeanvec;
		vector<double > posvarvec;
		vector<string > posallelevec;
		vector<string > posposvec; 
		
		int het_size ;
		
		for(int j=0 ; j<RepIDVec.size();j++){
			string Rep_prefix = RepIDVec[j];	
			vector<string>	rcvec = tio.dlread_vector( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.rc" )	;	
			cout << "merge:" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.rc" << endl;
			het_size = rcvec.size();
			if(j==0){
				//init the xmean vectors
				for(int q =0 ; q<2*rcvec.size(); q ++){
					vector<double> tmpvec; 
					posNRCvec.push_back(tmpvec);
				}
				for(int q =0 ; q<rcvec.size(); q ++){
					vector<string> tmpvec; 
					posrawvec.push_back(tmpvec);
				}
			}
			
			for(int k =0 ; k<rcvec.size() ; k++){
				vector<string> tokvec = tok.Tokenize(rcvec[k], "\t" );
				posposvec.push_back(tokvec[0]);
				posposvec.push_back(tokvec[0]);
				posallelevec.push_back(tokvec[1]);
				posallelevec.push_back(tokvec[2]);
				posrawvec[k].push_back(tokvec[3]);
				double M_rc, m_rc ;
				M_rc = atof(tokvec[14].c_str());
				m_rc = atof(tokvec[15].c_str());
				posNRCvec[2*k].push_back(M_rc); 
				posNRCvec[2*k+1].push_back(m_rc);
			}
			

		}	
		
		//calc xmean 
		
		for(int q=0 ; q<posNRCvec.size(); q++) {
			double sum =0;
			for(int t = 0 ; t< posNRCvec[q].size(); t++){
				
				sum += posNRCvec[q][t];
				
			}
			posmeanvec.push_back(sum/(double)posNRCvec[q].size() );
			//cout << endl;
		}
				
		//select allele for regression
		 
		
		for(int q =0 ; q< het_size; q++){
			
			 
			int cov =0 ;
			ss.str("");
			for(int g =0 ; g< posrawvec[q].size(); g++){
				ss<< posrawvec[q][g];
				
				if(g< posrawvec[q].size() -1){
					ss<<"_";
				} 
				
				vector<string> tmpvec = tok.Tokenize( posrawvec[q][g] , ",") ;
				for(int s = 0 ; s< tmpvec.size() ; s++){
					vector<string> tmp1vec = tok.Tokenize( tmpvec[s] , ":") ;
					if(tmp1vec.size()>1){
						cov += atoi(tmp1vec[1].c_str());
					}
				}
			}
			string raw_str = ss.str();
			
			if(cov >=20){
				ss.str("");
				if(posmeanvec[2*q] > posmeanvec[2*q+1]){
					ss<< chr_name << "_"<< posposvec[2*q]<< "_+" << "\t"<< chr_name<< "\t" << posposvec[2*q] << "\t+\t"  << raw_str << "\t" << cov << "\t" << posallelevec[2*q]<< "\t"<< posallelevec[2*q+1]<< "\t"<< posmeanvec[2*q] << "\t" <<posmeanvec[2*q+1] ; 
				}
				else{
					ss<< chr_name << "_"<< posposvec[2*q]<< "_+" << "\t"<< chr_name<< "\t" << posposvec[2*q] << "\t+\t" << raw_str << "\t" << cov << "\t" << posallelevec[2*q+1]<< "\t"<< posallelevec[2*q]<< "\t"<< posmeanvec[2*q+1] << "\t" <<posmeanvec[2*q] ;
				}
				
				testvec.push_back(ss.str());
			}
		}
//
//// merge RCs of neg SNPs
		vector<vector<double> > negNRCvec ;
		vector<vector<string> > negrawvec ;

		vector<double > negmeanvec;
		vector<double > negvarvec;
		vector<string > negallelevec;
		vector<string > negposvec; 
		
		for(int j=0 ; j<RepIDVec.size();j++){
			string Rep_prefix = RepIDVec[j];	
			vector<string>	rcvec = tio.dlread_vector( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_het.rc" )	;	
			het_size = rcvec.size();
			cout << "merge:" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_vcf.rc" << endl;
			if(j==0){
				//init the xmean vectors
				for(int q =0 ; q<2*rcvec.size(); q ++){
					vector<double> tmpvec; 
					negNRCvec.push_back(tmpvec);
				}
				
				for(int q =0 ; q<rcvec.size(); q ++){
					vector<string> tmpvec; 
					negrawvec.push_back(tmpvec);
				}
			}
			
			for(int k =0 ; k<rcvec.size() ; k++){
				vector<string> tokvec = tok.Tokenize(rcvec[k], "\t" );
			
				negposvec.push_back(tokvec[0]);
				negposvec.push_back(tokvec[0]);
				negallelevec.push_back(tokvec[1]);
				negallelevec.push_back(tokvec[2]);
				negrawvec[k].push_back(tokvec[3]);
				double M_rc, m_rc ;
				M_rc = atof(tokvec[14].c_str());
				m_rc = atof(tokvec[15].c_str());
				
				negNRCvec[2*k].push_back(M_rc); 
				negNRCvec[2*k+1].push_back(m_rc);
			
			}
			

		}	
		
	//		
		//calc xmean 
		
		for(int q=0 ; q<negNRCvec.size(); q++) {
			double sum =0;
			for(int t = 0 ; t< negNRCvec[q].size(); t++){
				
				sum += negNRCvec[q][t];
				
			}
			negmeanvec.push_back(sum/(double)negNRCvec[q].size() );
			//cout << endl;
		}
		
		//calc var 
		
		for(int q=0 ; q<negNRCvec.size(); q++) {
			double sum =0;
			for(int t = 0 ; t< negNRCvec[q].size(); t++){
				
				sum += (negNRCvec[q][t]-negmeanvec[q])*(negNRCvec[q][t]-negmeanvec[q]);
				
			}
			negvarvec.push_back(sum/(double)(negNRCvec[q].size()-1) );
			//cout << endl;
		}
		
		//select allele for regression
		
		for(int q =0 ; q< het_size; q++){
			
			int cov =0 ;
			ss.str("");
			for(int g =0 ; g< negrawvec[q].size(); g++){
				ss<< negrawvec[q][g];
				
				if(g< negrawvec[q].size() -1){
					ss<<"_";
				} 
				
				vector<string> tmpvec = tok.Tokenize( negrawvec[q][g] , ",") ;
				for(int s = 0 ; s< tmpvec.size() ; s++){
					vector<string> tmp1vec = tok.Tokenize( tmpvec[s] , ":") ;
					if(tmp1vec.size()>1){
						cov += atoi(tmp1vec[1].c_str());
					}
				}
			}
			string raw_str = ss.str();
			
			if(cov >20 ){
				ss.str("");
				if(negmeanvec[2*q] > negmeanvec[2*q+1]){
					ss<< chr_name << "_"<< negposvec[2*q]<< "_-" << "\t"<< chr_name<< "\t" << negposvec[2*q] << "\t-\t" << raw_str << "\t" << cov << "\t"<< negallelevec[2*q]<< "\t"<< negallelevec[2*q+1]<< "\t"<< negmeanvec[2*q] << "\t" <<negmeanvec[2*q+1] ; 
				}
				else{
					ss<< chr_name << "_"<< negposvec[2*q]<< "_-" << "\t"<< chr_name<< "\t" << negposvec[2*q] << "\t-\t" << raw_str << "\t" << cov << "\t" << negallelevec[2*q+1]<< "\t"<< negallelevec[2*q]<< "\t"<< negmeanvec[2*q+1] << "\t" <<negmeanvec[2*q] ;
				}
				
				testvec.push_back(ss.str());
			}
		}
//			
	}
//	
	tio.dlwrite( tmp_dir+"/test_vcf.rc" , testvec);
//	
	
	
}


void SNP::SNPtoRC_test(	std::string VCFfname, 	std::string tmp_dir, std::vector< std::string > RepIDVec , std::string chr_str, std::map<std::string, double> SizeMap){
	stringstream ss;
	TextIO tio;
	TOK tok;
	vector<string> outVec; 
	map<string, vector<string> > vcfmap , heteromap ;
	vector<string> chrVec = tok.Tokenize( chr_str, ",") ; 
	
	//read the VCF file and store by chrms
	
	cout << "reading biallelic SNVs from the VCF file :" << VCFfname << endl; 
	vector<string> vcfvec = tio.dlread_vector(VCFfname);
	for(int i =0 ; i < vcfvec.size() ; i++){
		if(vcfvec[i][0] != '#'){
					vector<string> tokvec = tok.Tokenize( vcfvec[i] , "\t");
					string chr_name = tokvec[0];
					vector<string> tmpvec = tok.Tokenize( tokvec[7] , ";");
					string af_str = tmpvec[1].substr(3); 
					vector<string> afvec = tok.Tokenize( af_str , ",");
					af_str = afvec[0];
					string snp_id = tokvec[0]+"_"+tokvec[1];
					double af = atof(af_str.c_str());
					if(tokvec[3].length()==1 && tokvec[4].length()==1){ // check if biallelic
						
						if(vcfmap.find(chr_name) == vcfmap.end()){
								vector<string> qqvec ;
								qqvec.push_back(vcfvec[i]);
								vcfmap.insert(make_pair( chr_name , qqvec ));
						}
						else{
								vcfmap[chr_name].push_back(vcfvec[i]);
						}
						
						
						
						if(af!=1){
						 
							if(heteromap.find(chr_name) == heteromap.end()){
								vector<string> qqvec ;
								qqvec.push_back(vcfvec[i]);
								heteromap.insert(make_pair( chr_name , qqvec ));
							}
							else{
								heteromap[chr_name].push_back(vcfvec[i]);
							}
						}	
					}
					
		}
	} 

	//for hetero. testing
	 
	cout << "for hetero. testing..... "<< endl; 	
	vector<string> testvec ; 
	testvec.push_back("id\tchr\tpos\tstrand\traw_str\tcov\tM_str\tm_str\tM\tm");
	for(int i =0 ; i<chrVec.size(); i++){
		string chr_name= chrVec[i];
		
		// ov annotated all SNPs to peak  
	  	  
	  for(int j=0 ; j<RepIDVec.size();j++){
  	
  		string Rep_prefix = RepIDVec[j];  
			
			listAllOvSNPs_VCF_lite( heteromap[chr_name], tmp_dir+"/"+Rep_prefix+"_"+chr_name+"_pos.bed", tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.snp" );			
			cout << "print : "<< tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.snp" << endl;
			listAllOvSNPs_VCF_lite( heteromap[chr_name], tmp_dir+"/"+Rep_prefix+"_"+chr_name+"_neg.bed", tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_het.snp"  );
			cout << "print : "<< tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_het.snp" << endl;
		}
		
//		//merge SNPs in rep
//		
		map<string, string> possnpmap,negsnpmap ;
		
//		
	  for(int j=0 ; j<RepIDVec.size();j++){
  	
  		string Rep_prefix = RepIDVec[j];  
			
			vector<string> posvec = tio.dlread_vector(tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.snp");
			vector<string> negvec = tio.dlread_vector(tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_het.snp");
			
			for(int k =0 ; k<posvec.size() ; k++){
				
				if(possnpmap.find(posvec[k]) == possnpmap.end()){
					possnpmap.insert(make_pair( posvec[k], ""));
				}
			}
			
			for(int k =0 ; k<negvec.size() ; k++){
				
				if(negsnpmap.find(negvec[k]) == negsnpmap.end()){
					negsnpmap.insert(make_pair( negvec[k], ""));
				}
			}
			
		}
////		
		//print the union of SNP 
		vector<string> univec; 
		for(map<string, string>::iterator mitr = possnpmap.begin(); mitr != possnpmap.end(); mitr++){
			univec.push_back(mitr->first);
		}
		
		tio.dlwrite( tmp_dir+"/"+chr_name+"_pos_hetov.snp" , univec);
		cout << "write: "<< tmp_dir+"/"+chr_name+"_pos_hetov.snp"<< endl;
		univec.clear();
		
		for(map<string, string>::iterator mitr = negsnpmap.begin(); mitr != negsnpmap.end(); mitr++){
			univec.push_back(mitr->first);
		}
		
		tio.dlwrite( tmp_dir+"/"+chr_name+"_neg_hetov.snp" , univec);
		cout << "write: "<< tmp_dir+"/"+chr_name+"_neg_hetov.snp"<< endl;
		univec.clear();
//		
//		
//		
//		//read read coverage in reps
		 
		for(int j=0 ; j<RepIDVec.size();j++){
			string Rep_prefix = RepIDVec[j];
			vector<string> SNPvec = tio.dlread_vector(tmp_dir+"/"+chr_name+"_pos_hetov.snp");
			readRawRC(SNPvec, tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+".mp" , tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawhet.rc");
			cout << "calc raw RC :" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawhet.rc" << endl;
			
			vector<string>  SNP_RepVec= normalize_VCF( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_rawhet.rc",  tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls", tmp_dir , 1 , SizeMap[Rep_prefix]);
			tio.dlwrite( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.rc" , SNP_RepVec);
			cout << "normalize raw RC :" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.rc\t"<< tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls" << "\t"<<SizeMap[Rep_prefix] << endl;
			
			SNPvec = tio.dlread_vector(tmp_dir+"/"+chr_name+"_neg_hetov.snp");
			readRawRC(SNPvec, tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+".mp" , tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawhet.rc");
			cout << "calc raw RC :" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawhet.rc" << endl;
			
			SNP_RepVec= normalize_VCF( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_rawhet.rc",  tmp_dir+"/"+Rep_prefix+"_neg_2_"+chr_name+".cls", tmp_dir , 0 , SizeMap[Rep_prefix]);
			tio.dlwrite( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_het.rc" , SNP_RepVec);
			cout << "normalize raw RC :" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_het.rc\t"<< tmp_dir+"/"+Rep_prefix+"_pos_2_"+chr_name+".cls" << "\t"<<SizeMap[Rep_prefix] << endl;
		}
//		
//		//merge RCs of pos SNPs
		vector<vector<double> > posNRCvec ;
		vector<vector<string> > posrawvec ;
		vector<double > posmeanvec;
		vector<double > posvarvec;
		vector<string > posallelevec;
		vector<string > posposvec; 
		
		int het_size ;
		
		for(int j=0 ; j<RepIDVec.size();j++){
			string Rep_prefix = RepIDVec[j];	
			vector<string>	rcvec = tio.dlread_vector( tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.rc" )	;	
			cout << "merge:" << tmp_dir+"/"+Rep_prefix+"_pos_"+chr_name+"_het.rc" << endl;
			het_size = rcvec.size();
			if(j==0){
				//init the xmean vectors
				for(int q =0 ; q<2*rcvec.size(); q ++){
					vector<double> tmpvec; 
					posNRCvec.push_back(tmpvec);
				}
				for(int q =0 ; q<rcvec.size(); q ++){
					vector<string> tmpvec; 
					posrawvec.push_back(tmpvec);
				}
			}
			
			for(int k =0 ; k<rcvec.size() ; k++){
				vector<string> tokvec = tok.Tokenize(rcvec[k], "\t" );
				posposvec.push_back(tokvec[0]);
				posposvec.push_back(tokvec[0]);
				posallelevec.push_back(tokvec[1]);
				posallelevec.push_back(tokvec[2]);
				posrawvec[k].push_back(tokvec[3]);
				double M_rc, m_rc ;
				M_rc = atof(tokvec[14].c_str());
				m_rc = atof(tokvec[15].c_str());
				posNRCvec[2*k].push_back(M_rc); 
				posNRCvec[2*k+1].push_back(m_rc);
			}
			

		}	
		
		//calc xmean 
		
		for(int q=0 ; q<posNRCvec.size(); q++) {
			double sum =0;
			for(int t = 0 ; t< posNRCvec[q].size(); t++){
				
				sum += posNRCvec[q][t];
				
			}
			posmeanvec.push_back(sum/(double)posNRCvec[q].size() );
			//cout << endl;
		}
				
		//select allele for regression
		 
		
		for(int q =0 ; q< het_size; q++){
			
			 
			int cov =0 ;
			ss.str("");
			for(int g =0 ; g< posrawvec[q].size(); g++){
				ss<< posrawvec[q][g];
				
				if(g< posrawvec[q].size() -1){
					ss<<"_";
				} 
				
				vector<string> tmpvec = tok.Tokenize( posrawvec[q][g] , ",") ;
				for(int s = 0 ; s< tmpvec.size() ; s++){
					vector<string> tmp1vec = tok.Tokenize( tmpvec[s] , ":") ;
					if(tmp1vec.size()>1){
						cov += atoi(tmp1vec[1].c_str());
					}
				}
			}
			string raw_str = ss.str();
			
			if(cov >=20){
				ss.str("");
				if(posmeanvec[2*q] > posmeanvec[2*q+1]){
					ss<< chr_name << "_"<< posposvec[2*q]<< "_+" << "\t"<< chr_name<< "\t" << posposvec[2*q] << "\t+\t"  << raw_str << "\t" << cov << "\t" << posallelevec[2*q]<< "\t"<< posallelevec[2*q+1]<< "\t"<< posmeanvec[2*q] << "\t" <<posmeanvec[2*q+1] ; 
				}
				else{
					ss<< chr_name << "_"<< posposvec[2*q]<< "_+" << "\t"<< chr_name<< "\t" << posposvec[2*q] << "\t+\t" << raw_str << "\t" << cov << "\t" << posallelevec[2*q+1]<< "\t"<< posallelevec[2*q]<< "\t"<< posmeanvec[2*q+1] << "\t" <<posmeanvec[2*q] ;
				}
				
				testvec.push_back(ss.str());
			}
		}
//
//// merge RCs of neg SNPs
		vector<vector<double> > negNRCvec ;
		vector<vector<string> > negrawvec ;

		vector<double > negmeanvec;
		vector<double > negvarvec;
		vector<string > negallelevec;
		vector<string > negposvec; 
		
		for(int j=0 ; j<RepIDVec.size();j++){
			string Rep_prefix = RepIDVec[j];	
			vector<string>	rcvec = tio.dlread_vector( tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_het.rc" )	;	
			het_size = rcvec.size();
			cout << "merge:" << tmp_dir+"/"+Rep_prefix+"_neg_"+chr_name+"_het.rc" << endl;
			if(j==0){
				//init the xmean vectors
				for(int q =0 ; q<2*rcvec.size(); q ++){
					vector<double> tmpvec; 
					negNRCvec.push_back(tmpvec);
				}
				
				for(int q =0 ; q<rcvec.size(); q ++){
					vector<string> tmpvec; 
					negrawvec.push_back(tmpvec);
				}
			}
			
			for(int k =0 ; k<rcvec.size() ; k++){
				vector<string> tokvec = tok.Tokenize(rcvec[k], "\t" );
			
				negposvec.push_back(tokvec[0]);
				negposvec.push_back(tokvec[0]);
				negallelevec.push_back(tokvec[1]);
				negallelevec.push_back(tokvec[2]);
				negrawvec[k].push_back(tokvec[3]);
				double M_rc, m_rc ;
				M_rc = atof(tokvec[14].c_str());
				m_rc = atof(tokvec[15].c_str());
				
				negNRCvec[2*k].push_back(M_rc); 
				negNRCvec[2*k+1].push_back(m_rc);
			
			}
			

		}	
		
	//		
		//calc xmean 
		
		for(int q=0 ; q<negNRCvec.size(); q++) {
			double sum =0;
			for(int t = 0 ; t< negNRCvec[q].size(); t++){
				
				sum += negNRCvec[q][t];
				
			}
			negmeanvec.push_back(sum/(double)negNRCvec[q].size() );
			//cout << endl;
		}
		
		//calc var 
		
		for(int q=0 ; q<negNRCvec.size(); q++) {
			double sum =0;
			for(int t = 0 ; t< negNRCvec[q].size(); t++){
				
				sum += (negNRCvec[q][t]-negmeanvec[q])*(negNRCvec[q][t]-negmeanvec[q]);
				
			}
			negvarvec.push_back(sum/(double)(negNRCvec[q].size()-1) );
			//cout << endl;
		}
		
		//select allele for testing
		
		for(int q =0 ; q< het_size; q++){
			
			int cov =0 ;
			ss.str("");
			for(int g =0 ; g< negrawvec[q].size(); g++){
				ss<< negrawvec[q][g];
				
				if(g< negrawvec[q].size() -1){
					ss<<"_";
				} 
				
				vector<string> tmpvec = tok.Tokenize( negrawvec[q][g] , ",") ;
				for(int s = 0 ; s< tmpvec.size() ; s++){
					vector<string> tmp1vec = tok.Tokenize( tmpvec[s] , ":") ;
					if(tmp1vec.size()>1){
						cov += atoi(tmp1vec[1].c_str());
					}
				}
			}
			string raw_str = ss.str();
			
			if(cov >20 ){
				ss.str("");
				if(negmeanvec[2*q] > negmeanvec[2*q+1]){
					ss<< chr_name << "_"<< negposvec[2*q]<< "_-" << "\t"<< chr_name<< "\t" << negposvec[2*q] << "\t-\t" << raw_str << "\t" << cov << "\t"<< negallelevec[2*q]<< "\t"<< negallelevec[2*q+1]<< "\t"<< negmeanvec[2*q] << "\t" <<negmeanvec[2*q+1] ; 
				}
				else{
					ss<< chr_name << "_"<< negposvec[2*q]<< "_-" << "\t"<< chr_name<< "\t" << negposvec[2*q] << "\t-\t" << raw_str << "\t" << cov << "\t" << negallelevec[2*q+1]<< "\t"<< negallelevec[2*q]<< "\t"<< negmeanvec[2*q+1] << "\t" <<negmeanvec[2*q] ;
				}
				
				testvec.push_back(ss.str());
			}
		}
//			
	}
//	
	tio.dlwrite( tmp_dir+"/test_vcf.rc" , testvec);
//	
	
	
}


void SNP::readRawRC_denovo(std::vector<std::string> ChrSNPVec, std::string mp_fname, std::string out_fname ){
	
	map<int,string> MPMap ;
	vector<string> outVec;
	stringstream ss;
	string str;
	fstream filestr;
	filestr.open(mp_fname.c_str(), ios_base::in);
	
	while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 	
			vector<string> tokvec = tok.Tokenize( str,"\t");					
			int mp_pos = atoi(tokvec[1].c_str()); 
			if(atoi(tokvec[3].c_str())>0){
				MPMap.insert( make_pair(mp_pos, str));
			}
		}
	}
	filestr.close();
	
	for( int i =0 ; i< ChrSNPVec.size() ; i++ ){
		//cout << "look at "<< m_itr->first << "\t" << m_itr->second<< endl;
		vector<string> tmptokvec = tok.Tokenize(ChrSNPVec[i]," ");
		int snp_pos = atoi(tmptokvec[1].c_str());
		
		map<string, double> chkSNPMap ; 
		
		chkSNPMap.insert(make_pair("A" , 0.0));
		chkSNPMap.insert(make_pair("T" , 0.0));
		chkSNPMap.insert(make_pair("C" , 0.0));
		chkSNPMap.insert(make_pair("G" , 0.0));
		
		double match_num =0;
		string ref_str = "-";
		if(MPMap.find(snp_pos) != MPMap.end()){
			//cout << MPMap[snp_pos] << endl;
		 	vector<string> tmpvec  = tok.Tokenize( MPMap[snp_pos] , "\t") ; 
		 	ref_str = tmpvec[2]; 
		 	string  base_str = tmpvec[4];
		 	string  score_str = tmpvec[5];
		 	vector<string> baseVec = _sepBaseStr(base_str);
		 	
		 	if(ref_str.compare("a")==0){
		 		ref_str="A";
		 	}
		 	if(ref_str.compare("t")==0){
		 		ref_str="T";
		 	}
		 	if(ref_str.compare("c")==0){
		 		ref_str="C";
		 	}
		 	if(ref_str.compare("g")==0){
		 		ref_str="G";
		 	}
		 	
		 	for( int j =0 ; j<baseVec.size() ; j++){
		 		int score =(int) score_str[j]; 
		 			if(baseVec[j].compare(".")==0 || baseVec[j].compare(",")==0 ){
		 				match_num = match_num +1.0;
		 				chkSNPMap[ref_str] = chkSNPMap[ref_str] + 1.0;
		 			}
		 			if(baseVec[j].compare("A")==0 || baseVec[j].compare("a")==0 ){
		 				if(chkSNPMap.find("A") != chkSNPMap.end()){
		 					chkSNPMap["A"] = chkSNPMap["A"] + 1.0;
		 				}
		 			}
		 			
		 			if(baseVec[j].compare("T")==0 || baseVec[j].compare("t")==0 ){
		 				if(chkSNPMap.find("T") != chkSNPMap.end()){
		 					chkSNPMap["T"] = chkSNPMap["T"] + 1.0;
		 				}
		 			}
		 			
		 			if(baseVec[j].compare("C")==0 || baseVec[j].compare("c")==0 ){
		 				if(chkSNPMap.find("C") != chkSNPMap.end()){
		 					chkSNPMap["C"] = chkSNPMap["C"] + 1.0;
		 				}
		 			}
		 			if(baseVec[j].compare("G")==0 || baseVec[j].compare("g")==0 ){
		 				if(chkSNPMap.find("G") != chkSNPMap.end()){
		 					chkSNPMap["G"] = chkSNPMap["G"] + 1.0;
		 				}
		 			}
		 				
		 		
		 		
		 	}
		 }
		 	ss.str("");
		 	ss<< snp_pos ;
		 	ss<< "\t"<< ref_str ;
		 	ss<< "\t-\t";
		 	double mut_num =0;
		 	for(map<string,double>::iterator m_itr = chkSNPMap.begin(); m_itr != chkSNPMap.end() ; m_itr++){
		 		
		 		ss<<  m_itr->first<< ":" << m_itr->second <<"," ;
		 		mut_num = mut_num+m_itr->second;
		 	}
		 	
		 	outVec.push_back(ss.str());
			
		
		
	}
	
	tio.dlwrite(out_fname, outVec);
	
	
}


void SNP::readRawRC(std::vector<std::string> ChrSNPVec, std::string mp_fname, std::string out_fname ){
	stringstream ss;
	string str;
	//vector<string> ChrSNPVec = tio.dlread_vector(snp_fname);
	vector<int> PosVec;
	vector<string> ASVec;
	
	map<int,string> SNPMap , SNPRefMap ; 
	
	map<int,string> MPMap ;
	
	vector<string> outVec;
	for(int i =0 ; i<ChrSNPVec.size(); i++){
		
		vector<string> tokvec = tok.Tokenize(ChrSNPVec[i], "\t");
		int snp_pos = atoi(tokvec[1].c_str());
		string alt_str = tokvec[4];
		string ref_str = tokvec[3];
		
		if(ref_str.compare("a")==0){
			ref_str= "A" ;
		}
		if(ref_str.compare("t")==0){
			ref_str= "T" ;
		}
		if(ref_str.compare("c")==0){
			ref_str= "C" ;
		}
		if(ref_str.compare("g")==0){
			ref_str= "G" ;
		}
		
		
		if(alt_str.compare("a")==0){
			alt_str= "A" ;
		}
		if(alt_str.compare("t")==0){
			alt_str= "T" ;
		}
		if(alt_str.compare("c")==0){
			alt_str= "C" ;
		}
		if(alt_str.compare("g")==0){
			alt_str= "G" ;
		}
		
		//snp_pos ++;
		
		SNPMap.insert(make_pair( snp_pos ,alt_str ));
		SNPRefMap.insert(make_pair( snp_pos ,ref_str ));	 
	}	
	
	fstream filestr;
	filestr.open(mp_fname.c_str(), ios_base::in);
	
	while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 	
			vector<string> tokvec = tok.Tokenize( str,"\t");					
			int mp_pos = atoi(tokvec[1].c_str()); 
			if(atoi(tokvec[3].c_str())>0){
				MPMap.insert( make_pair(mp_pos, str));
			}
		}
	}
	filestr.close();
	
	for(map<int,string>::iterator m_itr = SNPMap.begin() ; m_itr != SNPMap.end(); m_itr++){
		//cout << "look at "<< m_itr->first << "\t" << m_itr->second<< endl;
		int snp_pos = m_itr->first;
		string alt_str = m_itr->second;
		string ref_str = SNPRefMap[snp_pos];
		
		map<string, double> chkSNPMap ; 
		
		chkSNPMap.insert(make_pair("A" , 0.0));
		chkSNPMap.insert(make_pair("T" , 0.0));
		chkSNPMap.insert(make_pair("C" , 0.0));
		chkSNPMap.insert(make_pair("G" , 0.0));
		
		double match_num =0;
		
		if(MPMap.find(snp_pos) != MPMap.end()){
			//cout << MPMap[snp_pos] << endl;
		 	vector<string> tmpvec  = tok.Tokenize( MPMap[snp_pos] , "\t") ; 
		 	string tmpref_str = tmpvec[2]; 
		 	string  base_str = tmpvec[4];
		 	string  score_str = tmpvec[5];
		 	vector<string> baseVec = _sepBaseStr(base_str);
		 	
		 	if( tmpref_str.compare("a")==0){
		 		tmpref_str="A";
		 	}
		 	if( tmpref_str.compare("t")==0){
		 		tmpref_str="T";
		 	}
		 	if( tmpref_str.compare("c")==0){
		 		tmpref_str="C";
		 	}
		 	if( tmpref_str.compare("g")==0){
		 		tmpref_str="G";
		 	}
		 	
		 	for( int j =0 ; j<baseVec.size() ; j++){
		 		//int score =(int) score_str[j]; 

		 			if(baseVec[j].compare(".")==0 || baseVec[j].compare(",")==0 ){
		 			
		 				chkSNPMap[ tmpref_str] = chkSNPMap[ tmpref_str] + 1.0;
		 			}
		 			else{
		 				string mut_str="$";
		 				if(baseVec[j].compare("A")==0 || baseVec[j].compare("a")==0 ){
		 					mut_str="A";
		 				}
		 			
		 				if(baseVec[j].compare("T")==0 || baseVec[j].compare("t")==0 ){
		 					mut_str="T";
		 				}
		 			
		 				if(baseVec[j].compare("C")==0 || baseVec[j].compare("c")==0 ){
		 					mut_str="C";
		 				}
		 				if(baseVec[j].compare("G")==0 || baseVec[j].compare("g")==0 ){
		 					mut_str="G";
		 				}
		 				
		 				if(mut_str.compare("$") != 0){
		 					chkSNPMap[mut_str] = chkSNPMap[mut_str] + 1.0;
		 				}
		 				
		 			}	
		 		
		 	}
		 	
		}
		ss.str("");
		ss<< snp_pos ;
		ss<< "\t"<< ref_str << "\t" << alt_str;
		ss<< "\t";
		double mut_num =0;
		for(map<string,double>::iterator m_itr = chkSNPMap.begin(); m_itr != chkSNPMap.end() ; m_itr++){
		 		
			ss<<  m_itr->first<< ":" << m_itr->second <<"," ;
			mut_num = mut_num+m_itr->second;
		}
		 	
		outVec.push_back(ss.str());
			
		
		
	}
	
	tio.dlwrite(out_fname, outVec);
	
}

