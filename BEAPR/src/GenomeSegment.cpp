#include "GenomeSegment.h"

using namespace std;

std::map<string, string> GenomeSegment::parse_gid(std::string str){
	map<string, string> outMap ; 
	TOK tok;
	vector<string> tokvec  = tok.Tokenize( str, ";");

	for(int i =0 ; i < tokvec.size() ; i++){
	
		vector<string> tmpvec  = tok.Tokenize( tokvec[i] , " ");
		if(tmpvec.size()>1){
			outMap.insert(make_pair(tmpvec[0] ,tmpvec[1] ));
			//cout << "qq"<< tmpvec[0]<< ","<<tmpvec[1] << endl;
		}
	}
	
	return outMap; 
	
	
}

void GenomeSegment::_debug(){
	
	cout<< "PosVec : " ;
	for(int i=0;i<PosVec.size();i++){
		cout<< PosVec[i] << ",";
	}
	cout << endl;
	
	map<string, map<string,string > >::iterator itr;
	cout<< "ExonMap : " ;	
	for(itr = ExonMap.begin(); itr != ExonMap.end() ; itr ++){
		cout<< itr->first <<" : ";
		map<string, string > gmap = itr->second;
		map<string,string >::iterator mitr;  
		for(mitr = gmap.begin(); mitr != gmap.end(); mitr++){
			cout<< mitr->first <<" , ";
			
		}
		cout << endl;
	}
	
}

std::vector<double> GenomeSegment::readscov(std::string readsfname, int read_len){
	TextIO tio;
	TOK tok;
	vector<double> outDVec;
	map<int,int> OffsetMap;
	stringstream ss;
		
	//init outDVec;
	int offset =0;
	for(int i =0 ; i<PosVec.size()-1; i++){
		ss.str("");
		ss<<PosVec[i] << "_"<< PosVec[i+1];
		
		if(ExonMap.find(ss.str())!=ExonMap.end() ){
			map<string,string> nmap = ExonMap.find(ss.str())->second;
			if(nmap.size()>0){
				int len =PosVec[i+1] - PosVec[i] +1;
				//cout <<PosVec[i] << "_"<< PosVec[i+1] << " " << len << " "<< offset<< endl; 
				OffsetMap.insert(make_pair( PosVec[i], offset ));
				for(int j =0 ; j< len; j++){
					outDVec.push_back(0.0);
				}
			
				offset += len;
			}
		}
	}
	
	vector<vector<string> > readsVec = tio.dlread(readsfname, "\t" );  
	
	for(int i =0 ; i< readsVec.size() ; i++){
		int apos = atoi(readsVec[i][1].c_str())+1;
		string int_str = find_interval(apos) ;
		vector<string> idxvec = tok.Tokenize(int_str,"_");
		int int_aidx= atoi(idxvec[0].c_str());
		int int_bidx= atoi(idxvec[1].c_str());
		int int_apos = PosVec[int_aidx];
		int int_bpos = PosVec[int_bidx];
		
		if(OffsetMap.find(int_apos)!= OffsetMap.end()){
			int r_shift = OffsetMap.find(int_apos)->second;
			int spos = r_shift + apos - int_apos;
		
			//cout << apos << ":" << int_apos << "_" << int_bpos << ":" << r_shift << ":"<< r_shift + apos - int_apos << endl;
			for(int j=0; j< read_len ; j++){
				if(spos +j<outDVec.size()){
					outDVec[spos +j] +=1;
				} 
			}
		}
		 
	}
	
	return outDVec;
}
vector<std::string> GenomeSegment::get_intervals_by_pos(int s_pos, int e_pos){

	vector<std::string> outVec;
 
	string s_interval= find_interval(s_pos);
	string e_interval= find_interval(e_pos);
	
	stringstream ss;
	TOK tok;
	int sa,sb,ea,eb;
	vector<string> tokvec = tok.Tokenize(s_interval, "_");
	sa = atoi(tokvec[0].c_str());
	sb = atoi(tokvec[1].c_str());
	
	tokvec = tok.Tokenize(e_interval, "_");
	ea = atoi(tokvec[0].c_str());
	eb = atoi(tokvec[1].c_str());
	
	map<string, string> addedMap;
	
	for(int i =sa ; i<eb; i++){
		ss.str("");
		ss<< PosVec[i]<< "_" << PosVec[i]; 
		map<string , map<string,string> >::iterator exon_itr = ExonMap.find(ss.str());
		if(exon_itr!=ExonMap.end()){
			if((exon_itr->second).size()>0){
				if(addedMap.find(ss.str()) == addedMap.end()){
					outVec.push_back(ss.str());
					addedMap.insert(make_pair(ss.str(),""));
				}
			}
		}
		ss.str("");
		ss<< PosVec[i]<< "_" << PosVec[i+1]; 
		exon_itr = ExonMap.find(ss.str());
		if(exon_itr!=ExonMap.end()){
			if((exon_itr->second).size()>0){
				if(addedMap.find(ss.str()) == addedMap.end()){
					outVec.push_back(ss.str());
					addedMap.insert(make_pair(ss.str(),""));
				}
			}
		}
		
		ss.str("");
		ss<< PosVec[i+1]<< "_" << PosVec[i+1]; 
		exon_itr = ExonMap.find(ss.str());
		if(exon_itr!=ExonMap.end()){
			if((exon_itr->second).size()>0){
				if(addedMap.find(ss.str()) == addedMap.end()){
					outVec.push_back(ss.str());
					addedMap.insert(make_pair(ss.str(),""));
				}
			}
		}
		
	}
	
	return outVec;
	
}

void GenomeSegment::printReads(string dir_str){
	TextIO tio;
	map<string, map<string, string> >::iterator readmap_itr; 
	for(readmap_itr = ReadMap.begin();readmap_itr != ReadMap.end() ; readmap_itr++)
	{
		map<string, string> & rmap = readmap_itr->second;
		map<string, string>::iterator rmap_itr  ;
		vector<string> rvec;
		for(rmap_itr = rmap.begin();rmap_itr != rmap.end();rmap_itr++ ){
			rvec.push_back(rmap_itr->first);
		}
		tio.dlwrite(dir_str +"/"+readmap_itr -> first +".bed"  ,  rvec);
	}
}

std::vector<std::string> GenomeSegment::segCorr(){
	stringstream ss;
	vector<string> exonvec;
	for(int i =0 ; i<PosVec.size()-1 ;  i++){
		ss.str("");
		ss<< PosVec[i] << "_"<< PosVec[i+1];
		if(ExonMap.find(ss.str()) != ExonMap.end()){
			map<string,string> idmap =  ExonMap.find(ss.str())->second;
			if(idmap.size() > 0 ){
				ss.str("");
				ss<< PosVec[i] << "\t"<< PosVec[i+1];
				exonvec.push_back(ss.str());
				
			}
		}
	}
	
	return exonvec;
}	

void GenomeSegment::initGenomeSegment(string fname,int chr_len){
  vector<string> outVec;
 
	TextIO tio;
	TOK tok;
	stringstream ss;
   
  
  // init segs
  
  PosVec.push_back(0);
  PosVec.push_back(chr_len);
	
	map<string, string> init_gmap;
	ss.str("");
	ss<< 0<< "_"<< chr_len;
	ExonMap.insert(make_pair(ss.str(),init_gmap) );
	
	fstream filestr;
  
	filestr.open(fname.c_str(), ios_base::in);
	//  cout << "init GTF:"<<  fname  <<endl;	
	string str;	 
	while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
			
			vector<string> tokvec = tok.Tokenize(str,"\t");
			
//			cout << tokvec[2]<<endl;;
			
			if(tokvec[2].compare("exon") == 0 ){
			
				
				//cout << tokvec[8] <<endl;
				map<string, string> gidmap = parse_gid(tokvec[8]);
				//cout << gidmap.size() <<endl;
				string gname = (gidmap.find("transcript_id"))->second;
				
				
				//cout << gname << endl;
				
				
				
				int exon_spos = atoi(tokvec[3].c_str());
				int exon_epos = atoi(tokvec[4].c_str()); 
		    // find the 
		    string a_b = find_interval(exon_spos);
				vector<string> posvec = tok.Tokenize(a_b, "_");
				int sa_idx = atoi(posvec[0].c_str());
				int sb_idx = atoi(posvec[1].c_str());
				int sa_pos = PosVec[sa_idx];
				int sb_pos = PosVec[sb_idx];
		    
		    a_b = find_interval(exon_epos);
				posvec = tok.Tokenize(a_b, "_");
				int ea_idx = atoi(posvec[0].c_str());
				int eb_idx = atoi(posvec[1].c_str());
				int ea_pos = PosVec[ea_idx];
				int eb_pos = PosVec[eb_idx];
//initialize readmap

				if(ReadMap.find(gname)== ReadMap.end()){
					map<string,string> initMap;
					ReadMap.insert(make_pair(gname, initMap));
				}
				
//				cout << sa_pos << "," << sb_pos<<endl; 
//				cout << ea_pos << "," << eb_pos<<endl;
//				cout << exon_spos << "," << exon_epos<<endl;
        

				if(ea_idx == sa_idx){
					//split one interval into 3
					  ss.str("");
					  ss<< exon_spos<< "_"<< exon_epos;
					  string exon_name = ss.str();
					  
					  ss.str("");
					  ss<< sa_pos<< "_"<< sb_pos;
					  
					  string int_name = ss.str();
					  
					  map<std::string, map<std::string, string > >::iterator exonmap_itr = ExonMap.find(int_name ); 
						map<string, string> oldgene_map =  exonmap_itr->second;  
						map<string, string> newgene_map = oldgene_map;
						if(newgene_map.find(gname)==newgene_map.end() ){
							newgene_map.insert(make_pair(gname,"") );
						}
						
						//
						ExonMap.erase(int_name);
						//first part 
						
						int a_offset = 0;
						
						if(sa_pos < exon_spos){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+sb_idx,exon_spos);
							a_offset ++;
							if(exon_spos-sa_pos >1){
								vecitr = PosVec.begin();
								PosVec.insert(vecitr+sb_idx,exon_spos-1);
								a_offset ++;
							}
							ss.str("");
							ss<< sa_pos << "_"<< exon_spos-1; 
							ExonMap.insert(make_pair(ss.str(),oldgene_map) );	
						}
						//second part
						
							ExonMap.insert(make_pair(exon_name,newgene_map) );	
						
						//the 3rd part
						if(exon_epos < eb_pos){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+sb_idx+a_offset,exon_epos);
							a_offset ++;
							if(eb_pos - exon_epos >1){
								vecitr = PosVec.begin();
								PosVec.insert(vecitr+sb_idx+a_offset,exon_epos+1);
							}
							ss.str("");
							ss<< exon_epos+1 << "_"<< eb_pos; 
							ExonMap.insert(make_pair(ss.str(),oldgene_map) );	
							
						}
				}
				else{
					//update all internal intervals
					 
					for(int i = sb_idx+1; i<ea_idx; i++ ){
						int int_spos = PosVec[i];
						int int_epos = PosVec[i+1];
			
						ss.str("");
						ss<< int_spos<< "_"<< int_spos;
						
						map<string, map<string, string> >::iterator itr =  ExonMap.find(ss.str());
						if(itr !=  ExonMap.end()){	
							map<string, string> & gmap = itr->second;
							if(gmap.find(ss.str()) == gmap.end()){
								gmap.insert(make_pair( gname, "") );
							}	
						}
									
						
						ss.str("");
						ss<< int_spos<< "_"<< int_epos;
						
						itr =  ExonMap.find(ss.str());
						if(itr !=  ExonMap.end()){	
							map<string, string> & gmap = itr->second;
							if(gmap.find(ss.str()) == gmap.end()){
								gmap.insert(make_pair( gname, "") );
							}	
						}
							 
						ss.str("");
						ss<< int_epos<< "_"<< int_epos;
						
						itr =  ExonMap.find(ss.str());
						if(itr !=  ExonMap.end()){	
							map<string, string> & gmap = itr->second;
							if(gmap.find(ss.str()) == gmap.end()){
								gmap.insert(make_pair( gname, "") );
							}	
						}
			 
					}
					
					
					  ss.str("");
					  
						//create the first interval
						ss<< sa_pos<< "_"<< sb_pos;
					  
					  string int_name = ss.str();
					  
					  map<std::string, map<std::string, string > >::iterator exonmap_itr = ExonMap.find(int_name ); 
						if(exonmap_itr == ExonMap.end()){
							sa_pos = PosVec[sa_idx-1];
							sb_pos = PosVec[sb_idx-1];
							sa_idx--;
							sb_idx--;
							ss.str("");
							ss<< sa_pos<< "_"<< sb_pos;
							exonmap_itr = ExonMap.find(ss.str());
						}
						
						map<string, string> oldgene_map =  exonmap_itr->second;  
						map<string, string> newgene_map = oldgene_map;
						
						if(newgene_map.find(gname)==newgene_map.end() ){
									newgene_map.insert(make_pair(gname,"") );
						}
						
						//
							ExonMap.erase(int_name);
						
						if( sa_pos < exon_spos){
							ss.str("");
					 	  ss<< sa_pos<< "_"<< exon_spos-1;
					 	 
							ExonMap.insert(make_pair(ss.str(), oldgene_map) );
						}
						if( exon_spos <= sb_pos){
							ss.str("");
					 	  ss<< exon_spos<< "_"<< sb_pos;
					 	  ExonMap.insert(make_pair(ss.str(), newgene_map) );
						}
						
						
						
					
					//create the last interval
						
					  ss.str("");
					  
					  ss<< ea_pos<< "_"<< eb_pos;
					  
					  int_name = ss.str();
					  exonmap_itr = ExonMap.find(int_name ); 
					  //cout <<"test" << ss.str()<< endl;
					  if(exonmap_itr == ExonMap.end()){
							ea_pos = PosVec[ea_idx+1];
							eb_pos = PosVec[eb_idx+1];
							ea_idx++;
							eb_idx++;
							ss.str("");
							ss<< ea_pos<< "_"<< eb_pos;
							exonmap_itr = ExonMap.find(ss.str());
							int_name = ss.str();
							//cout << "change to :" << ss.str()<< endl;
						}
					  
						oldgene_map =  exonmap_itr->second;  
						newgene_map = oldgene_map;
						if(newgene_map.find(gname)==newgene_map.end() ){
							newgene_map.insert(make_pair(gname,"") );
						}
						ExonMap.erase(int_name);
						
						if(ea_pos <= exon_epos){
								ss.str("");
					  	  ss<< ea_pos<< "_"<< exon_epos;
					  	  ExonMap.insert(make_pair(ss.str(), newgene_map) );
						}
						if(exon_epos < eb_pos){
								ss.str("");
					  	  ss<< exon_epos+1<< "_"<< eb_pos;
					  	  ExonMap.insert(make_pair(ss.str(), oldgene_map) );
						}
						
						
					  //update posittions	
						int a_offset =0;
						
						if(sa_pos < exon_spos && exon_spos < sb_pos ){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+sb_idx,exon_spos);
							a_offset++;
						}
						
						if(sa_pos < exon_spos-1 && exon_spos-1 < sb_pos ){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+sb_idx,exon_spos-1);
							a_offset ++;
						}
						
						if(ea_pos < exon_epos && exon_epos < eb_pos ){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+eb_idx+a_offset,exon_epos);
							a_offset ++;
						}
						//cout << "test:"<< ea_pos<<"," << exon_epos+1 <<"," << eb_pos<<endl;
						if(ea_pos < exon_epos+1 && exon_epos+1 < eb_pos ){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+eb_idx+a_offset,exon_epos+1);
							//cout << "chk "<<endl;
						}				
				}
					
			}
				
			 
		}
				
	}
	
	//init ExonCountMap
	map<string , map<string, string> >::iterator exon_itr ;
	for(exon_itr = ExonMap.begin();exon_itr!=ExonMap.end() ; exon_itr++){
	//	cout <<exon_itr -> first << endl;
		ExonCountMap.insert(make_pair(exon_itr -> first , 0.0));
		
	}
	
	
}


void GenomeSegment::initGenomeSegment(string fname,int chr_len, string gname_chk){
  vector<string> outVec;
 
	TextIO tio;
	TOK tok;
	stringstream ss;
   
  
  // init segs
  
  PosVec.push_back(0);
  PosVec.push_back(chr_len);
	
	map<string, string> init_gmap;
	ss.str("");
	ss<< 0<< "_"<< chr_len;
	ExonMap.insert(make_pair(ss.str(),init_gmap) );
	
	fstream filestr;
  
	filestr.open(fname.c_str(), ios_base::in);
	//  cout << "init GTF:"<<  fname  <<endl;	
	string str;	 
	while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
			
			vector<string> tokvec = tok.Tokenize(str,"\t");
			
//			cout << tokvec[2]<<endl;;
			
			if(tokvec[2].compare("exon") == 0 ){
			
				
				//cout << tokvec[8] <<endl;
				map<string, string> gidmap = parse_gid(tokvec[8]);
				//cout << gidmap.size() <<endl;
				string gname = (gidmap.find(gname_chk))->second;
				
				gname = gname.substr(1, gname.length()-2);
				//cout << gname << endl;
				
				
				
				int exon_spos = atoi(tokvec[3].c_str());
				int exon_epos = atoi(tokvec[4].c_str()); 
		    // find the 
		    string a_b = find_interval(exon_spos);
				vector<string> posvec = tok.Tokenize(a_b, "_");
				int sa_idx = atoi(posvec[0].c_str());
				int sb_idx = atoi(posvec[1].c_str());
				int sa_pos = PosVec[sa_idx];
				int sb_pos = PosVec[sb_idx];
		    
		    a_b = find_interval(exon_epos);
				posvec = tok.Tokenize(a_b, "_");
				int ea_idx = atoi(posvec[0].c_str());
				int eb_idx = atoi(posvec[1].c_str());
				int ea_pos = PosVec[ea_idx];
				int eb_pos = PosVec[eb_idx];
//initialize readmap

				if(ReadMap.find(gname)== ReadMap.end()){
					map<string,string> initMap;
					ReadMap.insert(make_pair(gname, initMap));
				}
				
//				cout << sa_pos << "," << sb_pos<<endl; 
//				cout << ea_pos << "," << eb_pos<<endl;
//				cout << exon_spos << "," << exon_epos<<endl;
        

				if(ea_idx == sa_idx){
					//split one interval into 3
					  ss.str("");
					  ss<< exon_spos<< "_"<< exon_epos;
					  string exon_name = ss.str();
					  
					  ss.str("");
					  ss<< sa_pos<< "_"<< sb_pos;
					  
					  string int_name = ss.str();
					  
					  map<std::string, map<std::string, string > >::iterator exonmap_itr = ExonMap.find(int_name ); 
						map<string, string> oldgene_map =  exonmap_itr->second;  
						map<string, string> newgene_map = oldgene_map;
						if(newgene_map.find(gname)==newgene_map.end() ){
							newgene_map.insert(make_pair(gname,"") );
						}
						
						//
						ExonMap.erase(int_name);
						//first part 
						
						int a_offset = 0;
						
						if(sa_pos < exon_spos){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+sb_idx,exon_spos);
							a_offset ++;
							if(exon_spos-sa_pos >1){
								vecitr = PosVec.begin();
								PosVec.insert(vecitr+sb_idx,exon_spos-1);
								a_offset ++;
							}
							ss.str("");
							ss<< sa_pos << "_"<< exon_spos-1; 
							ExonMap.insert(make_pair(ss.str(),oldgene_map) );	
						}
						//second part
						
							ExonMap.insert(make_pair(exon_name,newgene_map) );	
						
						//the 3rd part
						if(exon_epos < eb_pos){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+sb_idx+a_offset,exon_epos);
							a_offset ++;
							if(eb_pos - exon_epos >1){
								vecitr = PosVec.begin();
								PosVec.insert(vecitr+sb_idx+a_offset,exon_epos+1);
							}
							ss.str("");
							ss<< exon_epos+1 << "_"<< eb_pos; 
							ExonMap.insert(make_pair(ss.str(),oldgene_map) );	
							
						}
				}
				else{
					//update all internal intervals
					 
					for(int i = sb_idx+1; i<ea_idx; i++ ){
						int int_spos = PosVec[i];
						int int_epos = PosVec[i+1];
			
						ss.str("");
						ss<< int_spos<< "_"<< int_spos;
						
						map<string, map<string, string> >::iterator itr =  ExonMap.find(ss.str());
						if(itr !=  ExonMap.end()){	
							map<string, string> & gmap = itr->second;
							if(gmap.find(ss.str()) == gmap.end()){
								gmap.insert(make_pair( gname, "") );
							}	
						}
									
						
						ss.str("");
						ss<< int_spos<< "_"<< int_epos;
						
						itr =  ExonMap.find(ss.str());
						if(itr !=  ExonMap.end()){	
							map<string, string> & gmap = itr->second;
							if(gmap.find(ss.str()) == gmap.end()){
								gmap.insert(make_pair( gname, "") );
							}	
						}
							 
						ss.str("");
						ss<< int_epos<< "_"<< int_epos;
						
						itr =  ExonMap.find(ss.str());
						if(itr !=  ExonMap.end()){	
							map<string, string> & gmap = itr->second;
							if(gmap.find(ss.str()) == gmap.end()){
								gmap.insert(make_pair( gname, "") );
							}	
						}
			 
					}
					
					
					  ss.str("");
					  
						//create the first interval
						ss<< sa_pos<< "_"<< sb_pos;
					  
					  string int_name = ss.str();
					  
					  map<std::string, map<std::string, string > >::iterator exonmap_itr = ExonMap.find(int_name ); 
						if(exonmap_itr == ExonMap.end()){
							sa_pos = PosVec[sa_idx-1];
							sb_pos = PosVec[sb_idx-1];
							sa_idx--;
							sb_idx--;
							ss.str("");
							ss<< sa_pos<< "_"<< sb_pos;
							exonmap_itr = ExonMap.find(ss.str());
						}
						
						map<string, string> oldgene_map =  exonmap_itr->second;  
						map<string, string> newgene_map = oldgene_map;
						
						if(newgene_map.find(gname)==newgene_map.end() ){
									newgene_map.insert(make_pair(gname,"") );
						}
						
						//
							ExonMap.erase(int_name);
						
						if( sa_pos < exon_spos){
							ss.str("");
					 	  ss<< sa_pos<< "_"<< exon_spos-1;
					 	 
							ExonMap.insert(make_pair(ss.str(), oldgene_map) );
						}
						if( exon_spos <= sb_pos){
							ss.str("");
					 	  ss<< exon_spos<< "_"<< sb_pos;
					 	  ExonMap.insert(make_pair(ss.str(), newgene_map) );
						}
						
						
						
					
					//create the last interval
						
					  ss.str("");
					  
					  ss<< ea_pos<< "_"<< eb_pos;
					  
					  int_name = ss.str();
					  exonmap_itr = ExonMap.find(int_name ); 
					  //cout <<"test" << ss.str()<< endl;
					  if(exonmap_itr == ExonMap.end()){
							ea_pos = PosVec[ea_idx+1];
							eb_pos = PosVec[eb_idx+1];
							ea_idx++;
							eb_idx++;
							ss.str("");
							ss<< ea_pos<< "_"<< eb_pos;
							exonmap_itr = ExonMap.find(ss.str());
							int_name = ss.str();
							//cout << "change to :" << ss.str()<< endl;
						}
					  
						oldgene_map =  exonmap_itr->second;  
						newgene_map = oldgene_map;
						if(newgene_map.find(gname)==newgene_map.end() ){
							newgene_map.insert(make_pair(gname,"") );
						}
						ExonMap.erase(int_name);
						
						if(ea_pos <= exon_epos){
								ss.str("");
					  	  ss<< ea_pos<< "_"<< exon_epos;
					  	  ExonMap.insert(make_pair(ss.str(), newgene_map) );
						}
						if(exon_epos < eb_pos){
								ss.str("");
					  	  ss<< exon_epos+1<< "_"<< eb_pos;
					  	  ExonMap.insert(make_pair(ss.str(), oldgene_map) );
						}
						
						
					  //update positions	
						int a_offset =0;
						
						if(sa_pos < exon_spos && exon_spos < sb_pos ){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+sb_idx,exon_spos);
							a_offset++;
						}
						
						if(sa_pos < exon_spos-1 && exon_spos-1 < sb_pos ){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+sb_idx,exon_spos-1);
							a_offset ++;
						}
						
						if(ea_pos < exon_epos && exon_epos < eb_pos ){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+eb_idx+a_offset,exon_epos);
							a_offset ++;
						}
						//cout << "test:"<< ea_pos<<"," << exon_epos+1 <<"," << eb_pos<<endl;
						if(ea_pos < exon_epos+1 && exon_epos+1 < eb_pos ){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+eb_idx+a_offset,exon_epos+1);
							//cout << "chk "<<endl;
						}				
				}
					
			}
				
			 
		}
				
	}
	
	//init ExonCountMap
	map<string , map<string, string> >::iterator exon_itr ;
	for(exon_itr = ExonMap.begin();exon_itr!=ExonMap.end() ; exon_itr++){
	//	cout <<exon_itr -> first << endl;
		ExonCountMap.insert(make_pair(exon_itr -> first , 0.0));
		
	}
	
	
}


void GenomeSegment::initGenomeSegment_pksbed(string fname,int chr_len){
  vector<string> outVec;
 
	TextIO tio;
	TOK tok;
	stringstream ss;
   
  
  // init segs
  
  PosVec.push_back(0);
  PosVec.push_back(chr_len);
	
	map<string, string> init_gmap;
	ss.str("");
	ss<< 0<< "_"<< chr_len;
	ExonMap.insert(make_pair(ss.str(),init_gmap) );
	
	fstream filestr;
  
	filestr.open(fname.c_str(), ios_base::in);
	//  cout << "init GTF:"<<  fname  <<endl;	
	string str;	 
	while(!filestr.eof()){
	 	getline(filestr,str);
	 	
	 	if(str.length()>0){	 		
			
			vector<string> tokvec = tok.Tokenize(str,"\t");
			
//			cout << tokvec[2]<<endl;;
			
			if(tokvec.size() > 1 ){
			
				
				//cout << tokvec[8] <<endl;
				
				string gname = tokvec[3];
				
				int exon_spos = atoi(tokvec[1].c_str())+1;
				int exon_epos = atoi(tokvec[2].c_str()); 
				
		    // find the 
		    string a_b = find_interval(exon_spos);
				vector<string> posvec = tok.Tokenize(a_b, "_");
				int sa_idx = atoi(posvec[0].c_str());
				int sb_idx = atoi(posvec[1].c_str());
				int sa_pos = PosVec[sa_idx];
				int sb_pos = PosVec[sb_idx];
		    
		    a_b = find_interval(exon_epos);
				posvec = tok.Tokenize(a_b, "_");
				int ea_idx = atoi(posvec[0].c_str());
				int eb_idx = atoi(posvec[1].c_str());
				int ea_pos = PosVec[ea_idx];
				int eb_pos = PosVec[eb_idx];
				//initialize readmap

				if(ReadMap.find(gname)== ReadMap.end()){
					map<string,string> initMap;
					ReadMap.insert(make_pair(gname, initMap));
				}
				

				if(ea_idx == sa_idx){
					
					
					  ss.str("");
					  ss<< exon_spos<< "_"<< exon_epos;
					  string exon_name = ss.str();
					  
					  ss.str("");
					  ss<< sa_pos<< "_"<< sb_pos;
					  
					  string int_name = ss.str();
					  map<std::string, map<std::string, string > >::iterator exonmap_itr = ExonMap.find(int_name ); 
						
						
					if((exonmap_itr->second).size()==0){
						//split one interval into 3
						map<string, string> oldgene_map =  exonmap_itr->second;  
						map<string, string> newgene_map = oldgene_map;
						if(newgene_map.find(gname)==newgene_map.end() ){
							newgene_map.insert(make_pair(gname,"") );
						}
						
						//
						ExonMap.erase(int_name);
						//first part 
						
						int a_offset = 0;
						
						if(sa_pos < exon_spos){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+sb_idx,exon_spos);
							a_offset ++;
							if(exon_spos-sa_pos >1){
								vecitr = PosVec.begin();
								PosVec.insert(vecitr+sb_idx,exon_spos-1);
								a_offset ++;
							}
							ss.str("");
							ss<< sa_pos << "_"<< exon_spos-1; 
							ExonMap.insert(make_pair(ss.str(),oldgene_map) );	
						}
						//second part
						
							ExonMap.insert(make_pair(exon_name,newgene_map) );	
						
						//the 3rd part
						if(exon_epos < eb_pos){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+sb_idx+a_offset,exon_epos);
							a_offset ++;
							if(eb_pos - exon_epos >1){
								vecitr = PosVec.begin();
								PosVec.insert(vecitr+sb_idx+a_offset,exon_epos+1);
							}
							ss.str("");
							ss<< exon_epos+1 << "_"<< eb_pos; 
							ExonMap.insert(make_pair(ss.str(),oldgene_map) );	
							
						}
						
					}
					else{
					//	cout << gname << " is contained." << endl;
					}
						
						
				}
				else{
					//update all internal intervals
					 
					for(int i = sb_idx+1; i<ea_idx; i++ ){
						int int_spos = PosVec[i];
						int int_epos = PosVec[i+1];
			
						ss.str("");
						ss<< int_spos<< "_"<< int_spos;
						
						map<string, map<string, string> >::iterator itr =  ExonMap.find(ss.str());
						if(itr !=  ExonMap.end()){	
							map<string, string> & gmap = itr->second;
							if(gmap.find(ss.str()) == gmap.end()){
								gmap.insert(make_pair( gname, "") );
							}	
						}
									
						
						ss.str("");
						ss<< int_spos<< "_"<< int_epos;
						
						itr =  ExonMap.find(ss.str());
						if(itr !=  ExonMap.end()){	
							map<string, string> & gmap = itr->second;
							if(gmap.find(ss.str()) == gmap.end()){
								gmap.insert(make_pair( gname, "") );
							}	
						}
							 
						ss.str("");
						ss<< int_epos<< "_"<< int_epos;
						
						itr =  ExonMap.find(ss.str());
						if(itr !=  ExonMap.end()){	
							map<string, string> & gmap = itr->second;
							if(gmap.find(ss.str()) == gmap.end()){
								gmap.insert(make_pair( gname, "") );
							}	
						}
			 
					}
					
					
					  ss.str("");
					  
						//create the first interval
						ss<< sa_pos<< "_"<< sb_pos;
					  
					  string int_name = ss.str();
					  
					  map<std::string, map<std::string, string > >::iterator exonmap_itr = ExonMap.find(int_name ); 
						if(exonmap_itr == ExonMap.end()){
							sa_pos = PosVec[sa_idx-1];
							sb_pos = PosVec[sb_idx-1];
							sa_idx--;
							sb_idx--;
							ss.str("");
							ss<< sa_pos<< "_"<< sb_pos;
							exonmap_itr = ExonMap.find(ss.str());
						}
						
						map<string, string> oldgene_map =  exonmap_itr->second;  
						map<string, string> newgene_map = oldgene_map;
						
						if(newgene_map.find(gname)==newgene_map.end() ){
									newgene_map.insert(make_pair(gname,"") );
						}
						
						//
							ExonMap.erase(int_name);
						
						if( sa_pos < exon_spos){
							ss.str("");
					 	  ss<< sa_pos<< "_"<< exon_spos-1;
					 	 
							ExonMap.insert(make_pair(ss.str(), oldgene_map) );
						}
						if( exon_spos <= sb_pos){
							ss.str("");
					 	  ss<< exon_spos<< "_"<< sb_pos;
					 	  ExonMap.insert(make_pair(ss.str(), newgene_map) );
						}
						
						
						
					
					//create the last interval
						
					  ss.str("");
					  
					  ss<< ea_pos<< "_"<< eb_pos;
					  
					  int_name = ss.str();
					  exonmap_itr = ExonMap.find(int_name ); 
					  //cout <<"test" << ss.str()<< endl;
					  if(exonmap_itr == ExonMap.end()){
							ea_pos = PosVec[ea_idx+1];
							eb_pos = PosVec[eb_idx+1];
							ea_idx++;
							eb_idx++;
							ss.str("");
							ss<< ea_pos<< "_"<< eb_pos;
							exonmap_itr = ExonMap.find(ss.str());
							int_name = ss.str();
							//cout << "change to :" << ss.str()<< endl;
						}
					  
						oldgene_map =  exonmap_itr->second;  
						newgene_map = oldgene_map;
						if(newgene_map.find(gname)==newgene_map.end() ){
							newgene_map.insert(make_pair(gname,"") );
						}
						ExonMap.erase(int_name);
						
						if(ea_pos <= exon_epos){
								ss.str("");
					  	  ss<< ea_pos<< "_"<< exon_epos;
					  	  ExonMap.insert(make_pair(ss.str(), newgene_map) );
						}
						if(exon_epos < eb_pos){
								ss.str("");
					  	  ss<< exon_epos+1<< "_"<< eb_pos;
					  	  ExonMap.insert(make_pair(ss.str(), oldgene_map) );
						}
						
						
					  //update posittions	
						int a_offset =0;
						
						if(sa_pos < exon_spos && exon_spos < sb_pos ){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+sb_idx,exon_spos);
							a_offset++;
						}
						
						if(sa_pos < exon_spos-1 && exon_spos-1 < sb_pos ){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+sb_idx,exon_spos-1);
							a_offset ++;
						}
						
						if(ea_pos < exon_epos && exon_epos < eb_pos ){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+eb_idx+a_offset,exon_epos);
							a_offset ++;
						}
						//cout << "test:"<< ea_pos<<"," << exon_epos+1 <<"," << eb_pos<<endl;
						if(ea_pos < exon_epos+1 && exon_epos+1 < eb_pos ){
							vector<int>::iterator vecitr = PosVec.begin();
							PosVec.insert(vecitr+eb_idx+a_offset,exon_epos+1);
							//cout << "chk "<<endl;
						}				
				}
				//cout << gname << ":" << exon_spos<< ","<< exon_epos << ","<< sa_pos<< "_"<< sb_pos<< ","<< ea_pos<< "_"<< eb_pos<< endl;

			}
				
			 
		}
				
	}
	
	//init ExonCountMap 
	/*
	map<string , map<string, string> >::iterator exon_itr ;
	for(exon_itr = ExonMap.begin();exon_itr!=ExonMap.end() ; exon_itr++){
		cout <<exon_itr -> first << endl;
		ExonCountMap.insert(make_pair(exon_itr -> first , 0.0));
		
	}
	*/
	
}

string GenomeSegment::find_interval(int exon_pos){
	
	stringstream ss;
	stringstream int_ss;
	if(PosVec.size() == 0){
		return "";
	}
	
	int a = 0;
	int b = PosVec.size()-1;
	
	while( (b-a)>1 ){
		int mid_idx =(b+a)/2; 
		if(PosVec[mid_idx]<exon_pos){
			a = mid_idx;
		}
		else{
			b= mid_idx;
		}
	}
	
	int c = b+1;
	
	ss.str("");
	int_ss.str("");
	int_ss << PosVec[a]<< "_"<< PosVec[b];
	if(ExonMap.find(int_ss.str())!= ExonMap.end()){
	
		ss.str("");
		ss<< a<< "_"<< b;
	}
	
	int_ss.str("");
	int_ss << PosVec[b]<< "_"<< PosVec[b];
	if(ExonMap.find(int_ss.str())!= ExonMap.end()){
		ss.str("");
		ss<< b<< "_"<< b;
	}
	
	if(c<PosVec.size()){
		int_ss.str("");
		int_ss << PosVec[b]<< "_"<< PosVec[c];
		if(ExonMap.find(int_ss.str())!= ExonMap.end()){
			ss.str("");
			ss<< b<< "_"<< c;
		}
	}
	
	return ss.str();
	
}
//
void GenomeSegment::sep_gtfs_by_genes(std::string fnamestr,std::string out_dir){
	TextIO tio;
	TOK tok;
	map< string, vector<string> > gtfmap;
	stringstream ss;

	fstream filestr;

	filestr.open(fnamestr.c_str(), ios_base::in);
	
	string str;	 
	 while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
	 		vector<string> tmpvec = tok.Tokenize(str,"\t" );
	 		
	 		if(tmpvec[2].compare("exon") == 0){
	 			map<string, string> gidmap = parse_gid(tmpvec[8]);
				//cout << gidmap.size() <<endl;
				string tmpname = (gidmap.find("gene_id"))->second;
				string gname = tmpname.substr(1,tmpname.length()-2);	
				if( gtfmap.find(gname) == gtfmap.end()){
					vector<string> tokvec;
					//cout << tmpvec[1] << endl;
					tokvec.push_back(str);
					gtfmap.insert(make_pair(gname,tokvec));
				}
				else{
				
					(gtfmap.find(gname)->second).push_back(str);
				
				}
			}
		}
	}
	
	filestr.close();
	
	for(map< string, vector<string> >::iterator mapitr = gtfmap.begin();mapitr != gtfmap.end() ; mapitr++){
		tio.dlwrite(out_dir+"/"+mapitr->first,  mapitr->second );
	}
	
	
	
}
