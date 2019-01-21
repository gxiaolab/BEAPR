#include "Reg.h"

using namespace std;

//std::vector<std::string > Reg::mergeRCs( std::string Rep_str ,std::string tmp_dir, std::string RBP_id, double min_rc ){
//	stringstream ss;
//	vector< string > outVec;
//	map< string, vector<double> > MajMap,MinMap;
//	map<string, double> M_meanMap,m_meanMap, VarMap  ;
//	 
//	map<string,string> idmap; 
//	map<string,string> SNPEndMap;
//	
//	map< string, vector<string> > RawRCMap; 
//	
//	std::vector<std::string > RCVec;
//	vector<string> RCNameVec = tok.Tokenize( Rep_str , ","); 
//	map<string, string> snpidmap ;
//	
//	for( int i=0 ; i<RCNameVec.size()  ; i++){
//		RCVec.push_back(tmp_dir+"/"+RCNameVec[i]+".rc");		
//	}
//		
//	
//	//get all common SNPs 
//	for(int i=0 ; i<RCVec.size(); i++){
//		vector<string> InVec = tio.dlread_vector(RCVec[i]);
//		map<string,string> tmp_map;
//		for(int j =0 ; j<InVec.size(); j++){
//			
//			vector<string> tokvec = tok.Tokenize( InVec[j], "\t");
//			ss.str("");
//			ss<< tokvec[1]<< "_"<< tokvec[2] << "_" << tokvec[5] ;
//			string idstr = ss.str(); 
//			if(i==0){
//				tmp_map.insert(make_pair( idstr, "" ));
//			}
//			else{
//				if(idmap.find(idstr)!= idmap.end()){
//					tmp_map.insert(make_pair( idstr, ""));
//				}
//			}
//			if(snpidmap.find(idstr)==snpidmap.end()){
//				snpidmap.insert(make_pair(idstr , tokvec[0]));
//				
//			}
//		}
//		idmap.clear();
//		idmap = tmp_map;
//	}
//	
//	//initialize MajMap and MinMap ;
//	
//	int rep_num = RCVec.size(); 
//	for(map<string,string>::iterator mitr = idmap.begin() ; mitr!= idmap.end() ; mitr++){
//		string id_str = mitr->first;
//		vector<double> tmpvec ;
//		for(int j=0 ; j<rep_num ; j++){
//			tmpvec.push_back(0.0);
//		}
//		MajMap.insert(make_pair( id_str,tmpvec ));
//		MinMap.insert(make_pair( id_str,tmpvec ));
//		M_meanMap.insert(make_pair( id_str,0.0 ));
//		m_meanMap.insert(make_pair( id_str,0.0 ));
//		VarMap.insert(make_pair( id_str,-1.0 )); 
//		
//		vector<string> tmpstrvec;
//		RawRCMap.insert(make_pair( id_str, tmpstrvec));
//		
//	}
//	
//	//determine the major and the minor allele
// 
//	map< string, vector<double> > ATCGMap ;
//	
//	for(map<string,string>::iterator mitr = idmap.begin() ; mitr!= idmap.end() ; mitr++ ){
//		vector<double> tmpvec ; 
//		for(int j =0 ; j<4;j++ ){
//			tmpvec.push_back(0.0);
//		}
//		ATCGMap.insert(make_pair( mitr->first , tmpvec) );
//		
//	}
//	
//	for(int i=0 ; i<RCVec.size(); i++){
//		vector<string> InVec = tio.dlread_vector(RCVec[i]);
//		for(int j =0 ; j<InVec.size(); j++){
//			vector<string> tokvec = tok.Tokenize( InVec[j], "\t");
//			ss.str("");
//			ss<< tokvec[1]<< "_"<< tokvec[2]<< "_" << tokvec[5];
//			string idstr = ss.str();
//			if(idmap.find(idstr) != idmap.end()){
//				
//				vector<double> & ATCGVec = ATCGMap[idstr];
//				
//				
//				ATCGVec[0] = ATCGVec[0] + atof(tokvec[8].c_str());
//				ATCGVec[1] = ATCGVec[1] + atof(tokvec[10].c_str());
//				ATCGVec[2] = ATCGVec[2] + atof(tokvec[12].c_str());
//				ATCGVec[3] = ATCGVec[3] + atof(tokvec[14].c_str());
//			}
//			
//		}
//		
//	}
//	// MIdxMap, mIdxMap for major , minor allele
//	map<string , int > MIdxMap, mIdxMap; 
//	for(map< string, vector<double> >::iterator m_itr =  ATCGMap.begin() ; m_itr !=  ATCGMap.end() ; m_itr++ ){
//		
//		vector<double> & atcgtmpvec = ATCGMap[m_itr->first]; 
//		double largest =-1;
//		double second =-1;
//		int lidx=0;
//		int sidx=0;
//		for(int q=0 ; q<atcgtmpvec.size() ; q++){
//			double tmpd = atcgtmpvec[q];
//			if (tmpd > largest) {
//  		  second = largest;
//  		  largest = tmpd;
//  		  sidx = lidx;
//  		  lidx = q;
// 			} else if (tmpd > second) {
//  			second = tmpd;
//  			sidx = q;
//  			
// 			}
//		}
//		MIdxMap.insert(make_pair(m_itr->first , lidx ));
//		mIdxMap.insert(make_pair(m_itr->first , sidx ));
//		//cout << m_itr->first << "@" << atcgtmpvec[0] << "," << atcgtmpvec[1]<< "," <<atcgtmpvec[2]  << "," << atcgtmpvec[3]<< "@" << MIdxMap[m_itr->first] << ","<< mIdxMap[m_itr->first] << endl;	
//	}
//	
//		
//	for(int i=0 ; i<RCVec.size(); i++){
//		vector<string> InVec = tio.dlread_vector(RCVec[i]);
//		for(int j =0 ; j<InVec.size(); j++){
//			vector<string> tokvec = tok.Tokenize( InVec[j], "\t");
//			ss.str("");
//			ss<< tokvec[1]<< "_"<< tokvec[2] <<"_" << tokvec[5];
//			string idstr = ss.str();
//			if(idmap.find(idstr) != idmap.end()){
//				
//				vector<double> & MajVec = MajMap[idstr];
//				vector<double> & MinVec = MinMap[idstr];
//				
//				MajVec[i] = atof(tokvec[8+2*MIdxMap[idstr]].c_str());
//				MinVec[i] = atof(tokvec[8+2*mIdxMap[idstr]].c_str());
//				
//				
//				if(i==0){
//					SNPEndMap.insert(make_pair(idstr,tokvec[5]));
//				}
//				RawRCMap[idstr].push_back(tokvec[4]);		
//			}
//			
//		}
//		
//	}
//	
////	
////	//output for regression
////	
//	ss.str("");
//	ss<<"snp_id\tid\tlabel\tallele\txmean\tcv2";
//	outVec.push_back(ss.str());
//	for(map<string,string>::iterator mitr = idmap.begin() ; mitr!= idmap.end() ; mitr++){
//		
//		string id_str = mitr->first;
//		ss.str("");
//		ss<<  snpidmap[id_str]<< "\t" << id_str<< "\t0\t1";
//		
//		//calc M mean 
//		vector<double> rcnumVec= MajMap[id_str];
//		double sum =0;
//		
//		for(int j =0 ; j<rcnumVec.size(); j++ ){
//	
//			sum+=rcnumVec[j];
//		}
//		double M_mean = sum/(double)rep_num;
//		ss<<"\t"<< M_mean; 
//		M_meanMap[id_str] = M_mean;
//		if(M_mean > 0.1 ){
//		//calc M var  
//			double var_sum =0 ; 
//			for(int j =0 ; j<rcnumVec.size(); j++ ){
//				double rc = rcnumVec[j];
//				var_sum+=(rc-M_mean)*(rc-M_mean);
//			}
//			double M_var = var_sum/(double)(rep_num-1);
//			ss<<"\t"<< M_var/(M_mean*M_mean); 
//			//ss<<"\t"<< M_var; 		
//				outVec.push_back(ss.str());
//		}
//		
//		ss.str("");
//		ss<<snpidmap[id_str]<< "\t" << id_str << "\t0\t0";
//		
//		rcnumVec= MinMap[id_str];
//		
//		//calc m mean
//		
//		sum =0;
//		
//		for(int j =0 ; j<rcnumVec.size(); j++ ){
//	
//			sum+=rcnumVec[j];
//		}
//		double m_mean = sum/(double)rep_num;
//		ss<<"\t"<< m_mean; 
//		m_meanMap[id_str] = m_mean;
//		if(m_mean >0.1  ){
//			
//			double var_sum =0 ;
//			
//			for(int j =0 ; j<rcnumVec.size(); j++ ){
//				double rc = rcnumVec[j];
//				
//				var_sum+=(rc-m_mean)*(rc-m_mean);
//			}
//		
//			double m_var = var_sum/(double)(rep_num-1);
//			ss<< "\t"<< m_var/(m_mean*m_mean);
//			//ss<< "\t"<< m_var; 
//			outVec.push_back(ss.str());
//			
//		}
//		
//	}
//	
//	tio.dlwrite(tmp_dir+"/reg.rc",outVec);
//	
//	string cmd = "Rscript Reg.r " + tmp_dir+"/reg.rc "+ tmp_dir +"/reg.rc.out";
////
//	system(cmd.c_str());
//	//read expected std
//	 	
//	vector<string> VarVec = tio.dlread_vector(tmp_dir +"/reg.rc.out");
//	
//	for(int i =1 ; i<VarVec.size(); i++){
//		vector<string> tokvec = tok.Tokenize( VarVec[i], "\t");
//		if(VarMap.find(tokvec[1]) != VarMap.end() ){
//			if(tokvec[3].compare("1")==0){	
//				double std =atof(tokvec[tokvec.size()-1].c_str()); 
//				VarMap[tokvec[1]] = std;
//			} 
//		}
//	}
//	
//	vector<string> predVec; 
//	predVec.push_back("RBP\tsnp_id\tchr\tpos\tstrand\trawstr\tM\tm\tM_mean\tm_mean\texp_std");
//	
//	for(map<string, double>::iterator mitr = M_meanMap.begin() ; mitr!=M_meanMap.end(); mitr++){
//		vector<string> tokvec = tok.Tokenize(mitr->first , "_");
//		ss.str("");
//		ss << RBP_id<< "\t" << snpidmap[mitr->first] <<  "\t"<< tokvec[0]<< "\t"<< tokvec[1] << "\t" << SNPEndMap[mitr->first]<<"\t";
//		
//		vector<string> & rcstrvec  =  RawRCMap[mitr->first];
//		
//		for(int q=0 ; q<rcstrvec.size(); q++){
//			ss<<rcstrvec[q];
//			if(q< rcstrvec.size()-1){
//				ss<<"_";
//			}
//		}
//		string M_allele; 
//		string m_allele;
//		if(MIdxMap[mitr->first] ==0 ){
//			M_allele= "A";
//		}else if(MIdxMap[mitr->first] ==1){
//			M_allele = "T";
//		}else if(MIdxMap[mitr->first] ==2){
//			M_allele = "C";
//		}else{
//			M_allele="G";
//		}
//		ss<< "\t"<< M_allele;
//		
//		if(mIdxMap[mitr->first] ==0 ){
//			m_allele= "A";
//		}else if(mIdxMap[mitr->first] ==1){
//			m_allele = "T";
//		}else if(mIdxMap[mitr->first] ==2){
//			m_allele = "C";
//		}else{
//			m_allele="G";
//		}
//		ss<< "\t"<< m_allele;
//		
//		ss<<"\t" <<mitr->second << "\t"<< m_meanMap[mitr->first] << "\t"<<VarMap[mitr->first] ;
//    
//
//      /*as long as there is a read of the minor allele, the SNP is considered as a hetrogeneous SNP. */
//		if((mitr->second+m_meanMap[mitr->first]) > min_rc && m_meanMap[mitr->first] >0 ){
//			
//			predVec.push_back(ss.str());
//		}
//	}
////	
//	tio.dlwrite( tmp_dir +"/pred.rc", predVec );
//	
////create the mask of filtering
//	map<string, int> posmask; 
//	//deleted by pos
//	map<string, vector<int> > chrsnpmap;
//
//	for(int i =1 ; i<predVec.size() ;i++){
//		
//		vector<string> tokvec = tok.Tokenize( predVec[i] , "\t");
//		
//		string chrname= tokvec[2];
//		int	pos = atoi(tokvec[3].c_str() );
//		
//		if(chrsnpmap.find(chrname) ==chrsnpmap.end() ){
//			vector<int> tmpvec;
//			tmpvec.push_back(pos);
//			chrsnpmap.insert(make_pair( chrname, tmpvec) ); 
//		 
//		}
//		else{
//			chrsnpmap[chrname].push_back(pos);
//			
//		}
//	}
//		
//	for(map< string, vector<int> >::iterator mitr = chrsnpmap.begin() ; mitr != chrsnpmap.end() ; mitr++){
//		
//		vector<int> PosVec = mitr->second;
//	  for(int j =0 ; j<PosVec.size() ; j++){
//	  	int a_pos = PosVec[j];
//	  	for(int k =0 ; k<PosVec.size() ; k++){
//	  		if( PosVec[k] != a_pos){
//	  			if( abs(PosVec[k]-PosVec[j]) < 35  ){
//	  				 ss.str("");
//	  				 ss<< mitr->first << "_"<< a_pos ;
//	  				 string snp_id = ss.str(); 
//						 if( posmask.find(snp_id) == posmask.end() ){
//						 	 
//						 	 posmask.insert(make_pair(snp_id , 0));
//						 	 
//						 }	
//	  			} 
//	  		}
//	  	} 
//	  }
//	   
//	}
//	
//	for(map< string, int >::iterator mitr = posmask.begin() ; mitr != posmask.end() ; mitr++){
//		cout << mitr->first << endl;
//	}
//	
//	
//	cmd = "Rscript Pred.r " + tmp_dir +"/pred.rc "+ tmp_dir +"/seqerr.cut";
//	system(cmd.c_str());
//	
//	vector<string> foutvec ;
//			
//	vector<string> finvec = tio.dlread_vector(tmp_dir +"/pred.rc.out"); 
//	
//	for(int i =0 ; i<finvec.size(); i++){
//		vector<string> tokvec = tok.Tokenize( finvec[i] , "\t");
//		ss.str("");
//		ss << tokvec[2]<< "_" << tokvec[3];
//		
//		string snp_id =  ss.str();
//		if(posmask.find(snp_id) != posmask.end()){
//			ss.str("");
//			for(int q =0 ; q<11 ; q++){
//				ss<< tokvec[q]<< "\t";
//			}
//			ss<< 1.0;
//			foutvec.push_back(ss.str());
//			}
//		else{
//			foutvec.push_back(finvec[i]);
//		}
//				
//	}
//	
//	
//	tio.dlwrite( tmp_dir +"/pred.rc.filtered.out", foutvec );
////	cout <<tmp_dir +"/pred.rc.filtered.out" << endl;
//			
//	return outVec;
//	
//}


std::vector<std::string > Reg::mergeRCs( std::string Rep_str ,std::string tmp_dir, std::string RBP_id, double min_rc ){
	stringstream ss;
	vector< string > outVec;
	map< string, vector<double> > MajMap,MinMap;
	map<string, double> M_meanMap,m_meanMap, VarMap  ;
	 
	map<string,string> idmap; 
	map<string,string> SNPEndMap;
	
	map< string, vector<string> > RawRCMap; 
	
	std::vector<std::string > RCVec;
	vector<string> RCNameVec = tok.Tokenize( Rep_str , ","); 
	map<string, string> snpidmap ;
	
	for( int i=0 ; i<RCNameVec.size()  ; i++){
		RCVec.push_back(tmp_dir+"/"+RCNameVec[i]+".rc");		
	}
		
	
	//get all common SNPs 
	for(int i=0 ; i<RCVec.size(); i++){
		vector<string> InVec = tio.dlread_vector(RCVec[i]);
		map<string,string> tmp_map;
		for(int j =0 ; j<InVec.size(); j++){
			
			vector<string> tokvec = tok.Tokenize( InVec[j], "\t");
			ss.str("");
			ss<< tokvec[1]<< "_"<< tokvec[2] << "_" << tokvec[5] ;
			string idstr = ss.str(); 
			if(i==0){
				tmp_map.insert(make_pair( idstr, "" ));
			}
			else{
				if(idmap.find(idstr)!= idmap.end()){
					tmp_map.insert(make_pair( idstr, ""));
				}
			}
			if(snpidmap.find(idstr)==snpidmap.end()){
				snpidmap.insert(make_pair(idstr , tokvec[0]));
				
			}
		}
		idmap.clear();
		idmap = tmp_map;
	}
	
	//initialize MajMap and MinMap ;
	
	int rep_num = RCVec.size(); 
	for(map<string,string>::iterator mitr = idmap.begin() ; mitr!= idmap.end() ; mitr++){
		string id_str = mitr->first;
		vector<double> tmpvec ;
		for(int j=0 ; j<rep_num ; j++){
			tmpvec.push_back(0.0);
		}
		MajMap.insert(make_pair( id_str,tmpvec ));
		MinMap.insert(make_pair( id_str,tmpvec ));
		M_meanMap.insert(make_pair( id_str,0.0 ));
		m_meanMap.insert(make_pair( id_str,0.0 ));
		VarMap.insert(make_pair( id_str,-1.0 )); 
		
		vector<string> tmpstrvec;
		RawRCMap.insert(make_pair( id_str, tmpstrvec));
		
	}
	
	//determine the major and the minor allele
 
	map< string, vector<double> > ATCGMap ;
	
	for(map<string,string>::iterator mitr = idmap.begin() ; mitr!= idmap.end() ; mitr++ ){
		vector<double> tmpvec ; 
		for(int j =0 ; j<4;j++ ){
			tmpvec.push_back(0.0);
		}
		ATCGMap.insert(make_pair( mitr->first , tmpvec) );
		
	}
	
	for(int i=0 ; i<RCVec.size(); i++){
		vector<string> InVec = tio.dlread_vector(RCVec[i]);
		for(int j =0 ; j<InVec.size(); j++){
			vector<string> tokvec = tok.Tokenize( InVec[j], "\t");
			ss.str("");
			ss<< tokvec[1]<< "_"<< tokvec[2]<< "_" << tokvec[5];
			string idstr = ss.str();
			if(idmap.find(idstr) != idmap.end()){
				
				vector<double> & ATCGVec = ATCGMap[idstr];
				
				
				ATCGVec[0] = ATCGVec[0] + atof(tokvec[8].c_str());
				ATCGVec[1] = ATCGVec[1] + atof(tokvec[10].c_str());
				ATCGVec[2] = ATCGVec[2] + atof(tokvec[12].c_str());
				ATCGVec[3] = ATCGVec[3] + atof(tokvec[14].c_str());
			}
			
		}
		
	}
	// MIdxMap, mIdxMap for major , minor allele
	map<string , int > MIdxMap, mIdxMap; 
	for(map< string, vector<double> >::iterator m_itr =  ATCGMap.begin() ; m_itr !=  ATCGMap.end() ; m_itr++ ){
		
		vector<double> & atcgtmpvec = ATCGMap[m_itr->first]; 
		double largest =-1;
		double second =-1;
		int lidx=0;
		int sidx=0;
		for(int q=0 ; q<atcgtmpvec.size() ; q++){
			double tmpd = atcgtmpvec[q];
			if (tmpd > largest) {
  		  second = largest;
  		  largest = tmpd;
  		  sidx = lidx;
  		  lidx = q;
 			} else if (tmpd > second) {
  			second = tmpd;
  			sidx = q;
  			
 			}
		}
		MIdxMap.insert(make_pair(m_itr->first , lidx ));
		mIdxMap.insert(make_pair(m_itr->first , sidx ));
		//cout << m_itr->first << "@" << atcgtmpvec[0] << "," << atcgtmpvec[1]<< "," <<atcgtmpvec[2]  << "," << atcgtmpvec[3]<< "@" << MIdxMap[m_itr->first] << ","<< mIdxMap[m_itr->first] << endl;	
	}
	
		
	for(int i=0 ; i<RCVec.size(); i++){
		vector<string> InVec = tio.dlread_vector(RCVec[i]);
		for(int j =0 ; j<InVec.size(); j++){
			vector<string> tokvec = tok.Tokenize( InVec[j], "\t");
			ss.str("");
			ss<< tokvec[1]<< "_"<< tokvec[2] <<"_" << tokvec[5];
			string idstr = ss.str();
			if(idmap.find(idstr) != idmap.end()){
				
				vector<double> & MajVec = MajMap[idstr];
				vector<double> & MinVec = MinMap[idstr];
				
				MajVec[i] = atof(tokvec[8+2*MIdxMap[idstr]].c_str());
				MinVec[i] = atof(tokvec[8+2*mIdxMap[idstr]].c_str());
				
				
				if(i==0){
					SNPEndMap.insert(make_pair(idstr,tokvec[5]));
				}
				RawRCMap[idstr].push_back(tokvec[4]);		
			}
			
		}
		
	}
	
//	
//	//output for regression
//	
	ss.str("");
	ss<<"snp_id\tid\tlabel\tallele\txmean\tcv2";
	outVec.push_back(ss.str());
	for(map<string,string>::iterator mitr = idmap.begin() ; mitr!= idmap.end() ; mitr++){
		
		string id_str = mitr->first;
		ss.str("");
		ss<<  snpidmap[id_str]<< "\t" << id_str<< "\t0\t1";
		
		//calc M mean 
		vector<double> rcnumVec= MajMap[id_str];
		double sum =0;
		
		for(int j =0 ; j<rcnumVec.size(); j++ ){
	
			sum+=rcnumVec[j];
		}
		double M_mean = sum/(double)rep_num;
		ss<<"\t"<< M_mean; 
		M_meanMap[id_str] = M_mean;
		
		
		
		//if(M_mean > 2 ){
		//calc M var  
		double var_sum =0 ; 
		for(int j =0 ; j<rcnumVec.size(); j++ ){
			double rc = rcnumVec[j];
			var_sum+=(rc-M_mean)*(rc-M_mean);
		}
		double M_var = var_sum/(double)(rep_num-1);
		ss<<"\t"<< M_var/(M_mean*M_mean); 
			//ss<<"\t"<< M_var; 		
				
		//}
		
		string MRegStr =ss.str() ;
		
		ss.str("");
		ss<<snpidmap[id_str]<< "\t" << id_str << "\t0\t0";
		
		rcnumVec= MinMap[id_str];
		
		//calc m mean
		
		sum =0;
		
		for(int j =0 ; j<rcnumVec.size(); j++ ){
	
			sum+=rcnumVec[j];
		}
		double m_mean = sum/(double)rep_num;
		ss<<"\t"<< m_mean; 
		m_meanMap[id_str] = m_mean;
		//if(m_mean >2  ){
			
		var_sum =0 ;
			
		for(int j =0 ; j<rcnumVec.size(); j++ ){
			double rc = rcnumVec[j];
				
			var_sum+=(rc-m_mean)*(rc-m_mean);
		}
		
		double m_var = var_sum/(double)(rep_num-1);
		ss<< "\t"<< m_var/(m_mean*m_mean);
			//ss<< "\t"<< m_var; 
		
		string mRegStr =ss.str() ;
				
		//}
		
		if(m_mean > 0 && M_mean+m_mean>2.0){
			outVec.push_back(MRegStr);
			if(m_mean > 1.0){
				outVec.push_back(mRegStr);
			}
		}
		
	}
	
	tio.dlwrite(tmp_dir+"/reg.rc",outVec);
	
	string cmd = "Rscript Reg.r " + tmp_dir+"/reg.rc "+ tmp_dir +"/reg.rc.out";
//
	cout<< cmd << endl;
	system(cmd.c_str());
	//read expected std
	 	
	vector<string> VarVec = tio.dlread_vector(tmp_dir +"/reg.rc.out");
	
	for(int i =1 ; i<VarVec.size(); i++){
		vector<string> tokvec = tok.Tokenize( VarVec[i], "\t");
		if(VarMap.find(tokvec[1]) != VarMap.end() ){
			if(tokvec[3].compare("1")==0){	
				double std =atof(tokvec[tokvec.size()-1].c_str()); 
				VarMap[tokvec[1]] = std;
			} 
		}
	}
	
	vector<string> predVec; 
	predVec.push_back("RBP\tsnp_id\tchr\tpos\tstrand\trawstr\tM\tm\tM_mean\tm_mean\texp_std");
	
	for(map<string, double>::iterator mitr = M_meanMap.begin() ; mitr!=M_meanMap.end(); mitr++){
		vector<string> tokvec = tok.Tokenize(mitr->first , "_");
		ss.str("");
		ss << RBP_id<< "\t" << snpidmap[mitr->first] <<  "\t"<< tokvec[0]<< "\t"<< tokvec[1] << "\t" << SNPEndMap[mitr->first]<<"\t";
		
		vector<string> & rcstrvec  =  RawRCMap[mitr->first];
		
		for(int q=0 ; q<rcstrvec.size(); q++){
			ss<<rcstrvec[q];
			if(q< rcstrvec.size()-1){
				ss<<"_";
			}
		}
		string M_allele; 
		string m_allele;
		if(MIdxMap[mitr->first] ==0 ){
			M_allele= "A";
		}else if(MIdxMap[mitr->first] ==1){
			M_allele = "T";
		}else if(MIdxMap[mitr->first] ==2){
			M_allele = "C";
		}else{
			M_allele="G";
		}
		ss<< "\t"<< M_allele;
		
		if(mIdxMap[mitr->first] ==0 ){
			m_allele= "A";
		}else if(mIdxMap[mitr->first] ==1){
			m_allele = "T";
		}else if(mIdxMap[mitr->first] ==2){
			m_allele = "C";
		}else{
			m_allele="G";
		}
		ss<< "\t"<< m_allele;
		
		ss<<"\t" <<mitr->second << "\t"<< m_meanMap[mitr->first] << "\t"<<VarMap[mitr->first] ;
    

      /*as long as there is a read of the minor allele, the SNP is considered as a hetrogeneous SNP. */
		if((mitr->second+m_meanMap[mitr->first]) > min_rc && m_meanMap[mitr->first] >0 ){
			
			predVec.push_back(ss.str());
		}
	}
//	
	tio.dlwrite( tmp_dir +"/pred.rc", predVec );
	
	cmd = "Rscript Pred.r " + tmp_dir +"/pred.rc "+ tmp_dir +"/seqerr.cut";
	system(cmd.c_str());
//			
	return outVec;
	
}



std::vector<std::string > Reg::calcVar(std::string rcfname ){
	
	stringstream ss;
	vector<string> outVec ; 
	vector<string> RCVec = tio.dlread_vector(rcfname);
	
	outVec.push_back("id\tlabel\tallele\txmean\tcv2");
	for(int i=0 ; i< RCVec.size(); i++){
		vector<string> tokvec = tok.Tokenize( RCVec[i],"\t" );
		string id = tokvec[0];
		string label = tokvec[1];
		int size = (tokvec.size() -2)/2;
		//Calc the mean and variance values of M
		
		//mean
		//M mean
		double sum =0;
		
		for(int j =2 ; j<tokvec.size(); j+=2 ){
	
			sum+=atof(tokvec[j].c_str());
		}
		double M_mean = sum/(double)size;
		
		//m mean
		sum =0;
		for(int j =3 ; j<tokvec.size(); j+=2 ){
			 
			sum+=atof(tokvec[j].c_str());
		}
		double m_mean = sum/(double)size;
		//variance
		
		//M variance
		double var_sum =0 ; 
		for(int j =2 ; j<tokvec.size(); j+=2 ){
			double rc = atof(tokvec[j].c_str());
			var_sum+=(rc-M_mean)*(rc-M_mean);
		}
		double M_var = var_sum/(double)size;
			
		//Calc the mean and variance values of m		

		var_sum =0 ; 
		for(int j =3 ; j<tokvec.size(); j+=2 ){
			double rc = atof(tokvec[j].c_str());
			var_sum+=(rc-m_mean)*(rc-m_mean);
		}
		double m_var = var_sum/(double)size;
		
		ss.str("");
		ss << id<< "\t"<< label<< "\t1\t"<< M_mean << "\t"<< M_var/(M_mean*M_mean);
		string M_str = ss.str(); 
		ss.str("");
		ss << id<< "\t"<< label<< "\t0\t"<< m_mean << "\t"<< m_var/(m_mean*m_mean);
		string m_str = ss.str();
		if(m_mean>5 || M_mean>5){
			outVec.push_back(M_str);
			outVec.push_back(m_str);
		}
	}
	
	
	return outVec;
	
}

void Reg::calcFDR(std::string dir_str){
	
	
	vector<string> predvec ;
	predvec.push_back("RBP\tsnp_id\tchr\tpos\tstrand\trawstr\tM\tm\tM_mean\tm_mean\texp_std\tpvals");
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir (dir_str.c_str())) != NULL) {
  /* print all the files and directories within directory */
 	  while ((ent = readdir (dir)) != NULL) {
    
   	 string predfname(ent->d_name) ; 
   	 if(predfname.length() > 5){
   		 cout <<predfname << endl;
   		 
   		 vector<string> invec = tio.dlread_vector( dir_str+"/"+  predfname);
   		 for(int i =1 ; i<invec.size(); i++){
   		 		vector<string> tokvec = tok.Tokenize(invec[i],"\t");
   		 		double p = atof(tokvec[tokvec.size()-1].c_str());
   		 		if(p<1){
   		 			predvec.push_back(invec[i]);
   		 		}
   		 }
   	 }
  	}
  	closedir (dir);
	} 
	else {
  /* could not open directory */
  	perror ("");
  	return;
	}
	
	tio.dlwrite(dir_str+"/RBPS_pred.txt", predvec);
	
	
}



