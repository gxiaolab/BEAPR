
#include "Motif.h"
using namespace std;


int main(int argc, char **argv){

	Motif mt ; 

	TextIO tio;
	TOK tok;
	stringstream ss; 	
	string cmd;
	string CMDfname= argv[1];

	vector<string> cmdvec = tio.dlread_vector(CMDfname);
	map<string,string> CMDMap ;
//	      
	for(int i =0 ; i<cmdvec.size(); i++){
	  	vector<string> tokvec = tok.Tokenize( cmdvec[i],"=");
	   	CMDMap.insert(make_pair( tokvec[0], tokvec[1]));
	      	
	}

	string asarp_name =  CMDMap["ASARP"];
	string predfname = CMDMap["PREDFILE"];
	string tabfname = CMDMap["TABFILE"];
 	string outfname = CMDMap["OUTFILE"];
  

	vector<string> ASBPredVec ;
	vector<string> outvec;
	vector<string> tmppredvec = tio.dlread_vector(predfname);
	
	//remove from ASE genes

		if(tmppredvec.size()>1){
			outvec.push_back(tmppredvec[0]);
			for(int j =1 ; j<tmppredvec.size() ; j++){
				vector<string> tokvec = tok.Tokenize( tmppredvec[j], "\t");
				//double fdr = atof(tokvec[tokvec.size()-1].c_str());
	 			//if(fdr < 0.1){
	 			ASBPredVec.push_back(tmppredvec[j]);
	 			//} 
	 		}
		
			 	
			ASBPredVec = mt._filter_DGE_SNP( ASBPredVec , asarp_name ,  tabfname );

			
			
		}


	map<string, vector<int> > mask_spos_map , mask_epos_map; 
	vector<string > repeatvec = tio.dlread_vector(CMDMap["REPEATMASK"]);
	

	vector<string > predvec = ASBPredVec;
// //	//read spos of repeat segs
	
	for(int i =1 ; i < repeatvec.size() ;i++){
		vector<string> tokvec = tok.Tokenize(repeatvec[i],"\t");
		string chr_name= tokvec[5];
		int seq_spos = atoi(tokvec[6].c_str());
		int seq_epos = atoi(tokvec[7].c_str());
		string rtype = tokvec[11];
		if(rtype.compare("Simple_repeat")==0){
			if(mask_spos_map.find(chr_name) == mask_spos_map.end()){
				vector<int> tmpvec ; 
				tmpvec.push_back(seq_spos);
				mask_spos_map.insert(make_pair(chr_name ,tmpvec ));
			
			}
			else{
			
				mask_spos_map[chr_name].push_back(seq_spos);
			}
			if(mask_epos_map.find(chr_name) == mask_epos_map.end()){
				vector<int> tmpvec ; 
				tmpvec.push_back(seq_epos);
				mask_epos_map.insert(make_pair(chr_name ,tmpvec ));
			}
			else{
				mask_epos_map[chr_name].push_back(seq_epos);
			}
		
		}
	}
	
	
	map<string , string> seqmap; 
	for(int i =1 ; i<=24; i++){
		
		string chrm_name ;
		
		if(i<=22){
			ss.str("");
			ss<<"chr"<<i;
			chrm_name = ss.str(); 
		}
		if(i==23){
			chrm_name ="chrX";
		}
		if(i==24){
			chrm_name ="chrY";
		}
		
		string seq_fname = CMDMap["FASTADIR"]+"/"+chrm_name+".fa";
		//cout << seq_fname << endl;
		
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
    
    seqmap.insert(make_pair( chrm_name, seq) );
 }

	for(int i =0 ; i < predvec.size() ;i++){
		int flag =0;
		vector<string> tokvec = tok.Tokenize(predvec[i],"\t");
 		string chr_name= tokvec[1];
		int asb_pos = atoi(tokvec[2].c_str());
		vector<int> & sposvec = mask_spos_map[chr_name];
		vector<int> & eposvec = mask_epos_map[chr_name];		
		
		for(int j =0 ; j<sposvec.size() ; j++){
			if( sposvec[j]<= asb_pos && asb_pos <=eposvec[j] ){
				//cout <<predvec[i] << endl;
				flag = 1; 
			}
		}
		
		
		string & fa_seq = seqmap[chr_name];
		
		string asb_seq = fa_seq.substr(asb_pos-5 ,9);
		
		transform(asb_seq.begin(), asb_seq.end(), asb_seq.begin(), ::toupper);
		
		string qstr1 = tokvec[6]+tokvec[6]+tokvec[6]+tokvec[6]+tokvec[6];
		string qstr2 = tokvec[7]+tokvec[7]+tokvec[7]+tokvec[7]+tokvec[7];
		
		if(asb_seq.find(qstr1) !=  string::npos ){
			flag = 1;
		}
		if(asb_seq.find(qstr2) !=  string::npos ){
			
			flag = 1;
		}
		
		if(flag ==0){
			outvec.push_back(predvec[i]);
		}
		
	}

// 	cout << "/home/giovas/eCLIP_download/"+cell_name+"/"+rbp_name+"_test/out/pred_vcf.rc_norepeat.txt" << endl;
	tio.dlwrite( CMDMap["TMPDIR"]+"/pred_vcf.rc_norepeat.txt", outvec);	

	vector<string> qfavec ;
 vector<string> norepeatvec =tio.dlread_vector(CMDMap["TMPDIR"]+"/pred_vcf.rc_norepeat.txt"); 
	for(int i =1 ; i < norepeatvec.size() ;i++){
		vector<string> tokvec = tok.Tokenize(norepeatvec[i],"\t");

		string asb_id = tokvec[0];
		string chr_name = tokvec[1];
		int asb_pos = atoi(tokvec[2].c_str());
		
		string & fa_seq = seqmap[chr_name];
		
		string asb_seq1 = fa_seq.substr(asb_pos-51 ,50);
		string asb_seq2 = fa_seq.substr(asb_pos ,50);
		ss.str("");
		ss << ">"<< asb_id << "_" << tokvec[6];
		qfavec.push_back(ss.str());
   ss.str("");
   ss <<asb_seq1 << tokvec[6] << asb_seq2 ;
   
   string qstr = ss.str(); 
   transform(qstr.begin(), qstr.end(), qstr.begin(), ::toupper);
		//cout <<qstr << endl;
		//cout <<qstr.substr(0,70) <<endl;
		qfavec.push_back(qstr.substr(0,70));
		//cout <<qstr.substr(70) <<endl;
		qfavec.push_back(qstr.substr(70));
		
		ss.str("");
		ss << ">"<< asb_id << "_" << tokvec[7];
		qfavec.push_back(ss.str());
		
		
		ss.str("");
   ss <<asb_seq1 << tokvec[7] << asb_seq2 ;
   qstr = ss.str(); 
   transform(qstr.begin(), qstr.end(), qstr.begin(), ::toupper);
		//cout <<qstr << endl;
		//cout <<qstr.substr(0,70) <<endl;
		qfavec.push_back(qstr.substr(0,70));
		//cout <<qstr.substr(70) <<endl;
		qfavec.push_back(qstr.substr(70));
	}	
	
	tio.dlwrite( CMDMap["TMPDIR"]+"/q.fa", qfavec);	

	ss.str("");
	ss<<"blat -t=dna -q=dna "<<CMDMap["FASTADIR"]<<"/hg19.fa  "<<CMDMap["TMPDIR"] <<"/q.fa "<< CMDMap["TMPDIR"] <<"/out.psl";
	cmd = ss.str();
	system(cmd.c_str());

// ////	
// //	//read & parse out.psl
// ////	
	vector<string> pslvec  =tio.dlread_vector(CMDMap["TMPDIR"]+"/out.psl") ;
//	
	map<string, vector<double> > ASB_idmap; 
	
	for(int i=5; i< pslvec.size(); i++){
		vector<string> tokvec = tok.Tokenize(pslvec[i],"\t");

		double match_num = atof(tokvec[0].c_str());
		string asb_id = tokvec[9]; 
		string strand = tokvec[8];
		if(strand.compare("+") ==0 ){
		  if(ASB_idmap.find(asb_id) == ASB_idmap.end() ){
		  	vector<double> tmpvec;
		  	tmpvec.push_back(match_num);
		  	ASB_idmap.insert(make_pair( asb_id, tmpvec ));
		  }
		  else{
		  	ASB_idmap[asb_id].push_back(match_num);
		  }
		}
	}
	
	map<string, string> tobeRM_Map;
	for(map<string, vector<double> >::iterator mitr= ASB_idmap.begin() ; mitr != ASB_idmap.end(); mitr++){
		
		vector<double> matchvec = mitr->second;
		if(matchvec.size()>1){
			sort(matchvec.begin(), matchvec.end());
			double rate = matchvec[matchvec.size()-2] /  matchvec[matchvec.size()-1];
			if(rate >= 0.95 ){
				//cout << mitr->first<< ","<<matchvec[matchvec.size()-1] << ","<< matchvec[matchvec.size()-2] << endl;
				vector<string> tokvec = tok.Tokenize(mitr->first ,"_");
				ss.str("");
				ss<< tokvec[0]<< "_"<< tokvec[1]<< "_" <<tokvec[2] << "_"<<tokvec[3] ;
				string tobeRM_id =ss.str(); 
				if(tobeRM_Map.find(tobeRM_id) == tobeRM_Map.end() ){
					tobeRM_Map.insert(make_pair( tobeRM_id, ""));
					
				}
			}
		}
	}
	
	vector<string> postoutvec ; 	
		postoutvec.push_back(norepeatvec[0]);
	for(int i =1 ; i < norepeatvec.size() ;i++){
		vector<string> tokvec = tok.Tokenize(norepeatvec[i],"\t");
		string rbp_name =  tokvec[0];
		string asb_id = tokvec[1];
		
		string test_id =rbp_name +"_"+asb_id; 
		if(tobeRM_Map.find(test_id) == tobeRM_Map.end()){
			postoutvec.push_back(norepeatvec[i]);
		}
		else{
			cout << "TBRM:"<< test_id<< endl;
		}
		
	}



map<string,string> flagmap; 

//checking AF (optional)

if(CMDMap.find("HTZONE") != CMDMap.end() ) {
	vector<string> zonevec = tio.dlread_vector(CMDMap["HTZONE"]);

	for(int i =1 ; i< postoutvec.size(); i++){
		vector<string> tokvec = tok.Tokenize( postoutvec[i],"\t");
		string asb_chr= tokvec[1];
		int asb_pos = atoi(tokvec[2].c_str());
		for(int j =0 ; j< zonevec.size() ; j++){
			vector<string> tokzonevec = tok.Tokenize( zonevec[j],"\t");
			if(asb_chr.compare(tokzonevec[0]) ==0 ){
				int apos = atoi(tokzonevec[1].c_str());
				int bpos = atoi(tokzonevec[2].c_str());
				if( apos <= asb_pos &&  asb_pos <= bpos){
					if(flagmap.find(tokvec[0]) == flagmap.end()){
						flagmap.insert(make_pair( tokvec[0],"" ));
						cout << tokvec[0] << " in " <<zonevec[j] <<endl;
					}
				}
			}
		}
	}	


}

// // //checking AF (optional)
if(CMDMap.find("MINAF") != CMDMap.end() ) {
	double minaf = atof(CMDMap["MINAF"].c_str());
	for(int i =1 ; i< postoutvec.size(); i++){
		vector<int> ATCGvec;
		for(int j =0 ; j< 4 ;j++){
			ATCGvec.push_back(0);
		}
		vector<string> tokvec = tok.Tokenize( postoutvec[i],"\t");
		vector<string> tokvec1 =  tok.Tokenize( tokvec[4],"_");
		for(int j =0 ; j<tokvec1.size() ; j++){
			vector<string> tokvec2 =  tok.Tokenize( tokvec1[j],",");
			for( int k =0 ;k <4; k++){
				vector<string> tokvec3 =  tok.Tokenize( tokvec2[k],":");
				int num = atoi(tokvec3[1].c_str());
				ATCGvec[k]=ATCGvec[k]+num;

			}
		}
		int total_num = 0;
		for(int j =0 ; j< 4 ;j++){
			 total_num += ATCGvec[j];
		}
		int idx = 0;
		if(tokvec[6].compare("C")==0)
		{
			idx = 1;
		}
		if(tokvec[6].compare("G")==0)
		{
			idx = 2;
		}
		if(tokvec[6].compare("T")==0)
		{
			idx = 3;
		}	
		int M = ATCGvec[idx];

		
		double af = (double)M/(double)total_num ;
		if(af <= minaf){
			flagmap.insert(make_pair( tokvec[0],"" ));
			
		}
	}
}


vector<string> finaloutvec;
finaloutvec.push_back(postoutvec[0]);
double minfdr = atof( CMDMap["FDR"].c_str() );
for(int i =1 ; i< postoutvec.size(); i++){
	vector<string> tokvec = tok.Tokenize( postoutvec[i],"\t");
	if(flagmap.find(tokvec[0]) == flagmap.end()){
		double fdr = atof(tokvec[tokvec.size()-1].c_str());
		
		if(fdr<=minfdr){
			finaloutvec.push_back(postoutvec[i]);
		}
	}

}

tio.dlwrite( CMDMap["OUTFILE"] , finaloutvec);


// 	tio.dlwrite( "/home/giovas/eCLIP_download/"+cell_name+"/"+rbp_name+"_test/out/pred_vcf.rc_postproc_95.txt", postoutvec);

	return 0; 
}

