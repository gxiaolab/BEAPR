#include "CrossLinker.h"

using namespace std;


map<int, double> CrossLinker::null_pdf(int len, int site_num, int itr_num){
	map<int, double> distrMap;
	//init
	
	srand(time(NULL));
	
	
	for(int i =0 ; i<=site_num; i++){
		distrMap.insert(make_pair( i, 0.0));
	}
	
	int cl_num=0; 
	for(int i=0; i<itr_num; i++){
		map<int, double> posMap;
		for(int i=0 ; i< len; i++){
			posMap.insert(make_pair( i, 0.0));
		}
		
		//reandom assign positions
		for(int j=0 ; j< site_num; j++){
			int pos =  rand()%len;
			posMap[pos] = posMap[pos]+1;
		}
		
		
		//calc number of site
		for(map<int,double>::iterator m_itr =posMap.begin() ; m_itr != posMap.end(); m_itr++){
			int height = m_itr->second;
			if(height>0){
				distrMap[height] = distrMap[height]+1;
				cl_num++;
			}
		}
		
		
		
		//site distr
	}
	
	for(map<int,double>::iterator m_itr = distrMap.begin() ; m_itr != distrMap.end() ; m_itr++){
			m_itr->second = m_itr->second/(double)cl_num;
	}
	
	return distrMap;
}


void CrossLinker::callCLs(std::string pk_fname, std::string mp_fname, std::string out_fname, double fdr_cut, char delimiter){
	//cout << pk_fname<< endl;
	//cout << mp_fname<< endl;
 	vector<string> PKVec = tio.dlread_vector(pk_fname);
 	vector<string> MP_Vec = tio.dlread_vector(mp_fname);
  vector<string> outVec;
 	map< int, vector<int> > SegPosMap,SegNumMap ;
 	stringstream ss;
 	//get max num 
  	
 	//init SegVec
 	 	
 	//calc all read starting sites
 	for(int i =0 ; i<MP_Vec.size(); i++){
 		vector<string> tokvec = tok.Tokenize(MP_Vec[i],"\t");
 		if(tokvec.size()>4){
 			string mp_str = tokvec[4];
 			int pos = atoi(tokvec[1].c_str());
 			int seq_num = pos/10000;
 			int site_num = 0 ;
 			for(int j =0 ; j< mp_str.length(); j++){
 				if(mp_str[j]== delimiter){
 					site_num++;
 				}
 			}
 			if(site_num>0){
 				if(SegPosMap.find(seq_num)==SegPosMap.end()){
 					vector<int> svec,tvec ;
 					svec.push_back(pos);
 					tvec.push_back(site_num);
 					SegPosMap.insert(make_pair(seq_num , svec));
 					SegNumMap.insert(make_pair(seq_num , tvec));
 					
 				}else{
 					SegPosMap[seq_num].push_back(pos);
 					SegNumMap[seq_num].push_back(site_num);
 				}
 				
 			}
 		}
 	}
 	MP_Vec.clear();
 	
 	
  
 	//
// 	for(map< int, vector<int> >::iterator m_itr = SegPosMap.begin() ; m_itr!=SegPosMap.end(); m_itr++){
// 		int segnum=  m_itr->first;
// 		for(int j =0 ; j<SegPosMap[segnum].size(); j++){
// 			cout << SegPosMap[segnum][j]<< "," << SegNumMap[segnum][j]<< endl;
// 		}  
// 	}
//

 	
 	for(int i =0 ; i< PKVec.size(); i++){
	// the start and end of the cluster
		vector<string> tokvec = tok.Tokenize( PKVec[i], "\t");
		int pk_start, pk_end; 
		map<int,double> obs_distMap, null_distMap;
		pk_start = atoi(tokvec[1].c_str());
		pk_end = atoi(tokvec[2].c_str());
		
		//calc the obs distribution 
		int seg_s = pk_start/10000;
		int seg_e = pk_end/10000;
		
		int max_rnum =0; 
		int total_slots_num =0;
		for(int seg = seg_s  ; seg<=seg_e; seg++){
			
			if(SegPosMap.find(seg)!=SegPosMap.end()){
				for(int j = 0 ; j<SegPosMap[seg].size(); j++){
					if( SegPosMap[seg][j]>=pk_start&& SegPosMap[seg][j]<=pk_end ){
						//cout << SegPosMap[seg][j]<< "," << SegNumMap[seg][j]<< endl;
						
						if(max_rnum< SegNumMap[seg][j]){
							max_rnum = SegNumMap[seg][j]; 
						}
						if(obs_distMap.find(SegNumMap[seg][j]) == obs_distMap.end() ){
							obs_distMap.insert(make_pair( SegNumMap[seg][j], 1.0));
							total_slots_num++;
						}
						else{
							obs_distMap[SegNumMap[seg][j]] = obs_distMap[SegNumMap[seg][j]] +1.0;
							total_slots_num++;
						}
					}
				}
			}
			
		}
//		
		//cout << "pk:" << pk_start <<"<--->" << pk_end << endl;
		//calc the null distribution
		if(max_rnum>3 ){
			//cout << PKVec[i] << endl;
			int len = pk_end - pk_start + 1;
			int site_num =0; 
		 
			for(map<int,double>::iterator m_itr  = obs_distMap.begin() ;m_itr  != obs_distMap.end() ; m_itr ++){
				site_num +=m_itr->first * m_itr->second;
				obs_distMap[m_itr->first] =  obs_distMap[m_itr->first] / (double) total_slots_num;
				 
			}
			
//			for(int j=1; j<=max_rnum; j++){
//				if(obs_distMap.find(j)!=obs_distMap.end()){
//					cout <<  j << ","<< obs_distMap[j] << endl;
//				}
//				else{
//					cout <<  j << ","<< 0 << endl;
//				}
//			}
							
			map<int,double> nullMap =  null_pdf(len,site_num,10);
//		
//			for(map<int,double>::iterator m_itr  = nullMap.begin() ;m_itr  != nullMap.end() ; m_itr ++){
//				
//				cout <<m_itr->first << ","<< nullMap[m_itr->first] << endl; 
//			}		
			
			int min_height = 0 ; 
			
			
			//calc the fdr rate
			
			for(int j =2 ; j<=max_rnum ; j++){
				double null_area =0 ; 
				double obs_area =0 ;
				for(int k=j ; k<=max_rnum; k++){
					if(nullMap.find(k)!=nullMap.end()){
						null_area += nullMap[k];
					}
					if(obs_distMap.find(k)!=obs_distMap.end()){
						obs_area += obs_distMap[k];
					}
				}
				double fdr = null_area/obs_area;
				//cout << "FDR:" << j<< ","<< fdr << endl;
				min_height = j;
				if (fdr< fdr_cut){
					break;
				}
			}
			
			for(int seg = seg_s  ; seg<=seg_e; seg++){
			
				if(SegPosMap.find(seg)!=SegPosMap.end()){
					for(int j = 0 ; j<SegPosMap[seg].size(); j++){
						if( SegPosMap[seg][j]>=pk_start&& SegPosMap[seg][j]<=pk_end ){
						//
							if(SegNumMap[seg][j]>=min_height){
								ss.str("");
								ss << SegPosMap[seg][j]<< "\t" << SegNumMap[seg][j] ;
								outVec.push_back(ss.str());	
							}
	
						}
					}
				}
			
			}
			
		}
	}
	
 	tio.dlwrite(out_fname,outVec ); 
}

std::vector<std::vector<double> > CrossLinker::calcPrior(std::vector< std::string>  CLfFnameVec, std::vector< std::string> SEQfnameVec, std::string out_dir, int half_size ,int pos_flag){
	
   vector<double> Avec,Tvec,Cvec,Gvec;
   vector<vector<double> >outVec;
   stringstream ss;
   
   
		for(int i=0; i<2*half_size+1; i++){
					Avec.push_back(0.0);
					Tvec.push_back(0.0);
					Cvec.push_back(0.0);
					Gvec.push_back(0.0);
		}
		     
	for(int q = 0  ; q< CLfFnameVec.size() ; q++){
		
		CrossLinker clk;
		 
		
		string seq_fname = SEQfnameVec[q];
    fstream filestr;
    string str;
    filestr.clear();
    filestr.open(seq_fname.c_str(), ios_base::in);
    assert(filestr.is_open());
    string seq ="";
        
    while(!filestr.eof()){
        getline(filestr,str);
        if(str.length()>0){
    		   if(str[0]=='>'){
                 	;                	
           }
           else{
             	seq+= str;
           }
                    
                    	
        }
    }
	

    filestr.close();
		
		map<int,int> posmap ; 
		vector<string> clinkvec = tio.dlread_vector(CLfFnameVec[q]);
		for(int i =0 ; i< clinkvec.size(); i++){
				
				vector<string> tokvec =  tok.Tokenize(clinkvec[i] ,"\t"); 
					//string chrname = tokvec[0];
				int pos = atoi(tokvec[0].c_str());
				if(posmap.find(pos)==posmap.end()){
					posmap.insert(make_pair(pos,0));
					string site_str;
					if(pos_flag ==1){
						 site_str = seq.substr(pos-half_size-2,2*half_size+1);
					}
					else{
						site_str = seq.substr(pos-half_size,2*half_size+1);
					}
					for(int j = 0 ; j< site_str.length(); j++){
						if(pos_flag==1){
							if(site_str[j] == 'a'|| site_str[j] == 'A'){
								Avec[j] =Avec[j]+1; 
							}
							if(site_str[j] == 't'|| site_str[j] == 'T'){
								Tvec[j] =Tvec[j]+1; 
							}
							if(site_str[j] == 'c'|| site_str[j] == 'C'){
								Cvec[j] =Cvec[j]+1; 
							}
							if(site_str[j] == 'g'|| site_str[j] == 'G'){
								Gvec[j] =Gvec[j]+1; 
							}
						}
						else{
							if(site_str[j] == 'a'|| site_str[j] == 'A'){
								Tvec[j] =Tvec[j]+1; 
							}
							if(site_str[j] == 't'|| site_str[j] == 'T'){
								Avec[j] =Avec[j]+1; 
							}
							if(site_str[j] == 'c'|| site_str[j] == 'C'){
								Gvec[j] =Gvec[j]+1; 
							}
							if(site_str[j] == 'g'|| site_str[j] == 'G'){
								Cvec[j] =Cvec[j]+1; 
							}
						}
					}
				}
		}
		posmap.clear();

	}
//	
	vector<string> AStrVec, TStrVec, CStrVec, GStrVec; 
		for(int i=0; i<2*half_size+ 1; i++){
			ss.str("");
			ss << i-half_size<<"," << Avec[i]/(double)(Avec[i]+Tvec[i]+Cvec[i]+Gvec[i]) ;
			AStrVec.push_back(ss.str());
		}
		for(int i=0; i<2*half_size+1; i++){
			ss.str("");
			ss << i-half_size<<"," << Tvec[i]/(double)(Avec[i]+Tvec[i]+Cvec[i]+Gvec[i]) ;
			TStrVec.push_back(ss.str());
		}
		for(int i=0; i<2*half_size+1; i++){
			ss.str("");
			ss << i-half_size<<"," << Cvec[i]/(double)(Avec[i]+Tvec[i]+Cvec[i]+Gvec[i]) ;
			CStrVec.push_back(ss.str());
		}
		for(int i=0; i<2*half_size+1; i++){
			ss.str("");
			ss << i-half_size<<"," << Gvec[i]/(double)(Avec[i]+Tvec[i]+Cvec[i]+Gvec[i]) ;
			GStrVec.push_back(ss.str());
		}
	if(pos_flag ==1 ){
	
//cout << (Avec[0]+Tvec[0]+Cvec[0]+Gvec[0]) << endl;
	
		tio.dlwrite(out_dir+"/a_pos.txt", AStrVec);
		tio.dlwrite(out_dir+"/t_pos.txt", TStrVec);
		tio.dlwrite(out_dir+"/c_pos.txt", CStrVec);
		tio.dlwrite(out_dir+"/g_pos.txt", GStrVec);
	}
	else{
		
		tio.dlwrite(out_dir+"/a_neg.txt", AStrVec);
		tio.dlwrite(out_dir+"/t_neg.txt", TStrVec);
		tio.dlwrite(out_dir+"/c_neg.txt", CStrVec);
		tio.dlwrite(out_dir+"/g_neg.txt", GStrVec);
	}
	outVec.push_back(Avec);
	outVec.push_back(Tvec);
	outVec.push_back(Cvec);
	outVec.push_back(Gvec);
	
	return outVec;					
}



std::vector<int> CrossLinker::_mpstr_calc(std::string mpstr){
	vector<int> outVec;
	int match_num = 0;
  int A_num = 0;
  int T_num = 0;
  int C_num = 0;
  int G_num = 0;
  
	for(int i =0 ; i< mpstr.length(); i++){
   		if(mpstr[i]=='.'||mpstr[i]==','){
        			match_num++;
   		}
   		if(mpstr[i]=='a'||mpstr[i]=='A'){
         			A_num++;
      }
      if(mpstr[i]=='t'||mpstr[i]=='T'){
         			T_num++;
      }
      if(mpstr[i]=='c'||mpstr[i]=='C'){
         			C_num++;
      }
      if(mpstr[i]=='g'||mpstr[i]=='G'){
         			G_num++;
      }
      
  }
  
	
	outVec.push_back(match_num);
  outVec.push_back(A_num);
  outVec.push_back(T_num);
  outVec.push_back(C_num);
  outVec.push_back(G_num);
	
	return outVec; 
}


void CrossLinker::predCLs(std::string bed_fname, std::string mp_fname , std::string out_fname, std::string strand_str , std::string FDR_str  ){
	
	double FDR = atof(FDR_str.c_str());
	char end_symbol ;

	if(strand_str.compare("-")==0){
		end_symbol = '$'; 
		
	}
	if(strand_str.compare("+")==0){
		end_symbol = '^';
	}
	callCLs(  bed_fname ,  mp_fname , out_fname, FDR, end_symbol);
}

void CrossLinker::prior(std::string cls_dir,string seq_dir , std::string out_dir, std::string chr_str,std::string size_str){
	vector<string> clsvec ;
	vector<string> seqvec ;
	vector<double> Avec,Tvec,Cvec,Gvec;
  stringstream ss;
  
  string cmd;
  int half_size = atoi(size_str.c_str());
  
   
  for(int i=0; i<2*half_size+1; i++){
				Avec.push_back(0.0);
				Tvec.push_back(0.0);
				Cvec.push_back(0.0);
				Gvec.push_back(0.0);
	}
	
	
	
	vector<string> chrVec = tok.Tokenize( chr_str , ","); 
	
	for(int i =0 ; i<chrVec.size(); i++){
		
		string chr_name= chrVec[i];  
		clsvec.push_back(cls_dir+"/input_pos_2_"+chr_name+".cls");
		cmd= "touch " +  cls_dir+"/input_pos_2_"+chr_name+".cls";
		system(cmd.c_str());
		seqvec.push_back(seq_dir+"/"+chr_name+".fa");
		
	}

	vector<vector<double> > posvec =  calcPrior(clsvec, seqvec, out_dir , half_size, 1);
	
	for(int j=0; j<2*half_size+1; j++){
				Avec[j]+=posvec[0][j];
				Tvec[j]+=posvec[1][j];
				Cvec[j]+=posvec[2][j];
				Gvec[j]+=posvec[3][j];
	}
	
	clsvec.clear();
	seqvec.clear();

	
for(int i =0 ; i<chrVec.size(); i++){
		
		string chr_name= chrVec[i];  
		clsvec.push_back(cls_dir+"/input_neg_2_"+chr_name+".cls");
		cmd= "touch " +  cls_dir+"/input_neg_2_"+chr_name+".cls";
		system(cmd.c_str());
		seqvec.push_back(seq_dir+"/"+chr_name+".fa");
		
	}

	vector<vector<double> > negvec =  calcPrior(clsvec, seqvec, out_dir ,half_size, 0);
	
	for(int j=2*half_size; j>=0; j--){
				Avec[j]+=negvec[0][j];
				Tvec[j]+=negvec[1][j];
				Cvec[j]+=negvec[2][j];
				Gvec[j]+=negvec[3][j];
	}
	
	clsvec.clear();
	seqvec.clear();

	vector<string> AStrVec, TStrVec, CStrVec, GStrVec; 
	
	
	for(int i=0; i<2*half_size+ 1; i++){
		ss.str("");
		ss << i-half_size<<"," << Avec[i]/(double)(Avec[i]+Tvec[i]+Cvec[i]+Gvec[i]) ;
		AStrVec.push_back(ss.str());
	}
	for(int i=0; i<2*half_size+1; i++){
		ss.str("");
		ss << i-half_size<<"," << Tvec[i]/(double)(Avec[i]+Tvec[i]+Cvec[i]+Gvec[i]) ;
		TStrVec.push_back(ss.str());
	}
	for(int i=0; i<2*half_size+1; i++){
		ss.str("");
		ss << i-half_size<<"," << Cvec[i]/(double)(Avec[i]+Tvec[i]+Cvec[i]+Gvec[i]) ;
		CStrVec.push_back(ss.str());
	}
	for(int i=0; i<2*half_size+1; i++){
		ss.str("");
		ss << i-half_size<<"," << Gvec[i]/(double)(Avec[i]+Tvec[i]+Cvec[i]+Gvec[i]) ;
		GStrVec.push_back(ss.str());
	}
//cout << (Avec[0]+Tvec[0]+Cvec[0]+Gvec[0]) << endl;
	
	tio.dlwrite(out_dir+"/a.txt", AStrVec);
	tio.dlwrite(out_dir+"/t.txt", TStrVec);
	tio.dlwrite(out_dir+"/c.txt", CStrVec);
	tio.dlwrite(out_dir+"/g.txt", GStrVec);
}