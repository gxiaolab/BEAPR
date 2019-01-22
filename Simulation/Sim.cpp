#include "Sim.h"

using namespace std;


void Sim::createTestMean(std::string outfname, std::map<int, double> cv2map, std::map<int, double> distmap, double libsize, double allele_bias){
	vector<string> outVec ; 
	stringstream ss ; 
	//allele bias
	double bias = allele_bias;
  
	ss.str("");
	ss << "chrm\tpos\tstrand\tlabel\tM_mean\tm_mean\tMcv2\tmcv2" ; 
	outVec.push_back(ss.str());
	srand (time(NULL));
  int sampled =0; 
	for(int i =1 ; i <10; i++){
		int group_samples ;
		if(i<9){
			group_samples = 5000*distmap[i];
                }	
		else{
			group_samples = 5000-sampled;
		} 
		 
		for(int j =0 ; j<group_samples; j++){
			
			double r;
			r = ((double) rand() / (RAND_MAX))  ;
			double mean;
			double n_mean;
			if(i < 9){		
				n_mean = r*4 + i*4;
				mean = (n_mean)*0.25*libsize;
			}
			if(i==9){
				n_mean = r*40 + i*4;
				mean = (n_mean)*0.25*libsize;
			}	
			ss.str("");
			ss << "chr1\t"<<sampled+1 << "\t+";
			
			if(sampled%10 ==0){
				
				ss<< "\t1";
		
   			ss<<"\t"<< mean * bias << "\t"<< mean*(1-bias) ;
   			int M_idx = n_mean* bias/4.0;
   			int m_idx = n_mean* (1-bias)/4.0;
   			if(M_idx>9){
   				M_idx = 9;
   			}
   			if(m_idx >9){
   				m_idx = 9;
   			}
   			
   			ss<< "\t"<<cv2map[M_idx] << "\t" << cv2map[m_idx];
			}
			else{
				ss<< "\t0";
			
				ss<<"\t"<< mean * 0.5 << "\t"<< mean*(0.5) ;
				int M_idx = n_mean* 0.5/4.0;
				if(M_idx>9){
   				M_idx = 9;
   			}
   			ss<< "\t"<<cv2map[M_idx] << "\t" << cv2map[M_idx];
			}
			outVec.push_back(ss.str());
			sampled++;
		}
	}
	cout << "ok:"<<outfname<<endl;
	tio.dlwrite(outfname , outVec );
	
}


void Sim::createOvTestMean(std::string outfname, std::map<int, double> cv2map, std::map<int, double> distmap, double libsize, double allele_bias){
	vector<string> outVec ; 
	stringstream ss ; 
	//allele bias
	double bias = allele_bias;
  
	ss.str("");
	ss << "chrm\tpos\tstrand\tlabel\tM_mean\tm_mean\tMcv2\tmcv2" ; 
	outVec.push_back(ss.str());
	srand (time(NULL));
  int sampled =0; 
	for(int i =1 ; i <=10 ; i++){
		 
		for(int j =0 ; j<500; j++){
			
			double r;
			r = ((double) rand() / (RAND_MAX))  ;
			double mean;
			double n_mean;
			if(i <= 9){		
				n_mean = r*4 + i*4;
				mean = (n_mean)*0.25*libsize;
			}
			if(i==10){
				n_mean = r*40 + i*4;
				mean = (n_mean)*0.25*libsize;
			}	
			ss.str("");
			ss << "chr1\t"<<sampled+1 << "\t+";
			
			if(j%10 ==0){
				
				ss<< "\t1";
		
   			ss<<"\t"<< mean * bias << "\t"<< mean*(1-bias) ;
   			int M_idx = n_mean* bias/4.0;
   			int m_idx = n_mean* (1-bias)/4.0;
   			if(M_idx>9){
   				M_idx = 9;
   			}
   			if(m_idx >9){
   				m_idx = 9;
   			}
   			
   			ss<< "\t"<<cv2map[M_idx] << "\t" << cv2map[m_idx];
			}
			else{
				ss<< "\t0";
			
				ss<<"\t"<< mean * 0.5 << "\t"<< mean*(0.5) ;
				int M_idx = n_mean* 0.5/4.0;
				if(M_idx>9){
   				M_idx = 9;
   			}
   			ss<< "\t"<<cv2map[M_idx] << "\t" << cv2map[M_idx];
			}
			outVec.push_back(ss.str());
			sampled++;
		}
	}
	cout << "ok:"<<outfname<<endl;
	tio.dlwrite(outfname , outVec );
	
}



std::map<int,double> Sim::learnCV2(std::string regfname){
	map<int, double> OutMap; 
	map<int,double> SumMap;
	map<int,double> NumMap;
	vector<string> RCVec; 
	vector<string> datavec = tio.dlread_vector(regfname);
	for( int i=0; i<10; i++ ){
		SumMap.insert(make_pair( i, 0.0));
		NumMap.insert(make_pair( i, 0.0));
	}
	for(int i =1 ; i< datavec.size() ; i++){
		vector<string> tokvec = tok.Tokenize(datavec[i],"\t");
		double xmean = atof(tokvec[3].c_str());
		if(xmean>2.0){
			double cv2 = atof(tokvec[4].c_str());
			int mean_idx = xmean/4.0;
			if(mean_idx>9){
				mean_idx = 9;
			}
			SumMap[mean_idx] = SumMap[mean_idx] +cv2;
			NumMap[mean_idx] = NumMap[mean_idx] +1.0;
		}
	}
	for(int i=0; i<10; i++){
		OutMap.insert(make_pair(i,SumMap[i] / NumMap[i]));
		cout << i<< "\t"<< SumMap[i] / NumMap[i] << endl;
	}
	
	return OutMap;
}


std::map<int,double> Sim::learnDist(std::string predfname){
		
	map<int,double> NumMap;
	vector<string> RCVec; 
	vector<string> datavec = tio.dlread_vector(predfname);
	for( int i=-1; i<10; i++ ){
		NumMap.insert(make_pair( i, 0.0));
	}
	
	for(int i =1 ; i< datavec.size() ; i++){
		vector<string> tokvec = tok.Tokenize(datavec[i],"\t");
		double xmean = atof(tokvec[6].c_str())+atof(tokvec[7].c_str());
		if(xmean>4.0){
		int mean_idx = xmean/4.0;
			
			if(mean_idx>9){
				mean_idx = 9;
			}
			NumMap[mean_idx] = NumMap[mean_idx] +1.0;
		}
	}
	
	double sum= 0;
	for(int i=0; i<10; i++){
		cout << i<< "\t"<<  NumMap[i] << endl;
		sum= sum+NumMap[i];
	}
	
	for(int i=0; i<10; i++){
		cout << i<< "\t"<<  (int) 100*NumMap[i]/sum << endl;
		NumMap[i] = NumMap[i]/sum;
	}
	
	return NumMap;
}

std::map<int,double> Sim::learnDist_new(std::string predfname){
		
	map<int,double> NumMap;
	vector<string> RCVec; 
	vector<string> datavec = tio.dlread_vector(predfname);
	for( int i=-1; i<10; i++ ){
		NumMap.insert(make_pair( i, 0.0));
	}
	
	for(int i =0 ; i< datavec.size() ; i++){
		vector<string> tokvec = tok.Tokenize(datavec[i],"\t");
		int idx = atoi(tokvec[0].c_str());
		double xmean = atof(tokvec[1].c_str());
		NumMap[idx]=xmean;
		cout << i<< "\t"<<  NumMap[i] << endl;
	}
	
	return NumMap;
}




void Sim::createReg(std::string regfname, std::string predfname, std::string outfname){
	//remove the snp sites to be predicted from the reg file
	vector<string> RegVec = tio.dlread_vector(regfname);
	vector<string> PredVec = tio.dlread_vector(predfname);
	vector<string> OutVec ; 
	stringstream ss; 
	
	map<string, string> IDMap ; 
	
	for(int i =1 ; i< PredVec.size(); i++){
	  vector<string> tokvec = tok.Tokenize( PredVec[i] , "\t");
	  ss.str("");
	  ss<<tokvec[0] << "_"<< tokvec[1];
	  IDMap.insert(make_pair( ss.str(), ""));
  }
  OutVec.push_back(RegVec[0]);
  for(int i =1 ; i< RegVec.size(); i++){
  	vector<string> tokvec =tok.Tokenize( RegVec[i] , "\t") ;
  	if(IDMap.find(tokvec[0]) == IDMap.end()){
  		OutVec.push_back(RegVec[i]);
  	} 
  }
   	
	tio.dlwrite(outfname, OutVec);
}

