#include "Prep.h"

using namespace std;

void Prep::sep_mp_byChrms(std::string  mp_fname , std::string tmpdir , std::string prefix ){
	
//	
  fstream filestr;
  string str;
  vector<string> mpvec ;
  string old_chr="";
        
  filestr.clear();
  filestr.open(mp_fname.c_str(), ios_base::in);
  assert(filestr.is_open());
        
        while(!filestr.eof()){
              getline(filestr,str);
              if(str.length()>0){
                    vector<string> tokvec = tok.Tokenize(str,"\t");
                    if(old_chr.length()==0){
                    	old_chr= tokvec[0];
                    	
                    }
                    if(tokvec[0].compare(old_chr)!=0){
                    	tio.dlwrite(tmpdir+"/" + prefix+"_"+ old_chr+".mp",mpvec );
                    	old_chr= tokvec[0];
                    	mpvec.clear();
                    }
                    
                    mpvec.push_back(str);
                   
              }
        }
				
				tio.dlwrite(tmpdir+"/" + prefix+"_"+ old_chr+".mp",mpvec );
        mpvec.clear();
        filestr.close();
	
	
}
void Prep::sep_beds_by_chrs(std::string inBedfname, std::string tmpdir,std::string prefix){
	      vector<string> PKVec = tio.dlread_vector(inBedfname);
	      
	      map<string,vector<string> > PosBedMap,NegBedMap;
	      for(int i =0; i<PKVec.size(); i++){
	      	vector<string> tokVec = tok.Tokenize( PKVec[i], "\t");
	      	if(tokVec[5].compare("+")==0){
	      		if(PosBedMap.find(tokVec[0]) == PosBedMap.end()){
	      			vector<string> tmpvec ; 
	      			tmpvec.push_back(PKVec[i]);
	      			PosBedMap.insert(make_pair(tokVec[0] , tmpvec));
	      		}
	      		else{
	      			PosBedMap[tokVec[0]].push_back(PKVec[i]);
	      		}	
	      	}
	      	else{
	      		if(NegBedMap.find(tokVec[0]) == NegBedMap.end()){
	      			vector<string> tmpvec ; 
	      			tmpvec.push_back(PKVec[i]);
	      			NegBedMap.insert(make_pair(tokVec[0] , tmpvec));
	      		}
	      		else{
	      			NegBedMap[tokVec[0]].push_back(PKVec[i]);
	      		}
	      	}
	      }
	      
	     for(map<string,vector<string> >::iterator m_itr = PosBedMap.begin() ; m_itr != PosBedMap.end(); m_itr++){
	      	tio.dlwrite( tmpdir+"/"+prefix+"_"+ m_itr->first+"_pos.bed", m_itr->second);
	     }
	      
	     for(map<string,vector<string> >::iterator m_itr = NegBedMap.begin() ; m_itr != NegBedMap.end(); m_itr++){
	      	tio.dlwrite( tmpdir+"/"+prefix+"_"+ m_itr->first+"_neg.bed", m_itr->second);
	     }

}

void Prep::sep_beds_by_chrs(std::string inBedfname, std::string tmpdir,std::string prefix, double rc){
	      vector<string> PKVec = tio.dlread_vector(inBedfname);
	      
	      map<string,vector<string> > PosBedMap,NegBedMap;
	      for(int i =0; i<PKVec.size(); i++){
	      	vector<string> tokVec = tok.Tokenize( PKVec[i], "\t");
	      	double tmp_rc = atof(tokVec[4].c_str());
	      	if(tmp_rc > rc){
	      		if(tokVec[5].compare("+")==0){
	      			if(PosBedMap.find(tokVec[0]) == PosBedMap.end()){
	      				vector<string> tmpvec ; 
	      				tmpvec.push_back(PKVec[i]);
	      				PosBedMap.insert(make_pair(tokVec[0] , tmpvec));
	      			}
	      			else{
	      				PosBedMap[tokVec[0]].push_back(PKVec[i]);
	      			}	
	      		}
	      		else{
	      			if(NegBedMap.find(tokVec[0]) == NegBedMap.end()){
	      				vector<string> tmpvec ; 
	      				tmpvec.push_back(PKVec[i]);
	      				NegBedMap.insert(make_pair(tokVec[0] , tmpvec));
	      			}
	      			else{
	      				NegBedMap[tokVec[0]].push_back(PKVec[i]);
	      			}
	      		}
	      	}
	  	  }
	      
	     for(map<string,vector<string> >::iterator m_itr = PosBedMap.begin() ; m_itr != PosBedMap.end(); m_itr++){
	      	tio.dlwrite( tmpdir+"/"+prefix+"_"+ m_itr->first+"_pos.bed", m_itr->second);
	     }
	      
	     for(map<string,vector<string> >::iterator m_itr = NegBedMap.begin() ; m_itr != NegBedMap.end(); m_itr++){
	      	tio.dlwrite( tmpdir+"/"+prefix+"_"+ m_itr->first+"_neg.bed", m_itr->second);
	     }

}


void Prep::preProc(std::string inBamfname ,std::string Ref_fname ,std::string tmpdir ,std::string prefix){
	
				string cmd;


//for reads in the .bam format
	
//	//extract forward fragments (+)


	cmd = "samtools view -b -f 163 "+ inBamfname + " > " +tmpdir+"/"+prefix+"_163.bam" ;
	cout << cmd<< endl;
  system(cmd.c_str());
	
	cmd = "samtools view -b -f 83 "+ inBamfname + " > " +tmpdir+"/"+prefix+"_83.bam" ;
  cout << cmd<< endl;
  system(cmd.c_str());
	
	
	cmd = "samtools merge "+ tmpdir+"/"+prefix + "_pos.bam " +tmpdir+"/"+prefix+"_83.bam " +tmpdir+"/"+prefix+"_163.bam" ;
  cout << cmd<< endl;
  system(cmd.c_str());
	
	cmd = "samtools mpileup  -f "+Ref_fname+" "+ tmpdir+"/"+prefix + "_pos.bam " + " > " +tmpdir+"/"+prefix+"_pos.mp" ;
	cout << cmd<< endl;
	system(cmd.c_str());
	
	cmd = "samtools mpileup  -f "+Ref_fname+" "+ tmpdir+"/"+prefix+"_163.bam " + " > " +tmpdir+"/"+prefix+"_pos_2.mp" ;
	cout << cmd<< endl;
	system(cmd.c_str());

	
	//extract reverse fragments (-)
	cmd = "samtools view -b -f 99 "+ inBamfname + " > " +tmpdir+"/"+prefix+"_99.bam" ;
  cout << cmd<< endl;
  system(cmd.c_str());

	cmd = "samtools view -b -f 147 "+ inBamfname + " > " +tmpdir+"/"+prefix+"_147.bam" ;
  cout << cmd<< endl;
  system(cmd.c_str());
  
  cmd = "samtools merge "+ tmpdir+"/"+prefix + "_neg.bam " +tmpdir+"/"+prefix+"_99.bam " +tmpdir+"/"+prefix+"_147.bam" ;
  cout << cmd<< endl;
  system(cmd.c_str());
	
	cmd = "samtools mpileup -f "+Ref_fname+" "+ tmpdir+"/"+prefix + "_neg.bam " + " > " +tmpdir+"/"+prefix+"_neg.mp" ;
	cout << cmd<< endl;
	system(cmd.c_str());
	
	cmd = "samtools mpileup -f "+Ref_fname+" "+ tmpdir+"/"+prefix+"_147.bam " + " > " +tmpdir+"/"+prefix+"_neg_2.mp" ;
	cout << cmd<< endl;
	system(cmd.c_str());

	//sep .mp file by chrms
	
  sep_mp_byChrms(tmpdir+"/"+prefix+"_neg.mp" ,  tmpdir , prefix+"_neg" );
  sep_mp_byChrms(tmpdir+"/"+prefix+"_neg_2.mp" ,  tmpdir , prefix+"_neg_2" );
  sep_mp_byChrms(tmpdir+"/"+prefix+"_pos.mp" ,  tmpdir , prefix+"_pos" );
  sep_mp_byChrms(tmpdir+"/"+prefix+"_pos_2.mp" ,  tmpdir , prefix+"_pos_2" );
	
	
	//
	cmd = "rm -rf "+tmpdir+"/"+prefix+"_147.bam"+" "+tmpdir+"/"+prefix+"_99.bam"+" "+tmpdir+"/"+prefix+"_163.bam"+" "+tmpdir+"/"+prefix+"_83.bam";
	cout << cmd<< endl;
	system(cmd.c_str());	 
	cmd = "rm -rf "+tmpdir+"/"+prefix+"_neg.mp"+" "+tmpdir+"/"+prefix+"_neg_2.mp"+" "+tmpdir+"/"+prefix+"_pos.mp"+" "+tmpdir+"/"+prefix+"_pos_2.mp";
	cout << cmd<< endl;
	system(cmd.c_str());
//	
}







void Prep::sep_reads_by_chrs(std::string fnamestr , std::string out_dir ,  std::string prefix){
	 fstream filestr;
  string str;
  vector<string> mpvec ;
  string old_chr="";
        
  filestr.clear();
  filestr.open(fnamestr .c_str(), ios_base::in);
  assert(filestr.is_open());
        
        while(!filestr.eof()){
              getline(filestr,str);
              if(str.length()>0){
                    vector<string> tokvec = tok.Tokenize(str,"\t");
                    if(old_chr.length()==0){
                    	old_chr= tokvec[2];
                    	
                    }
                    if(tokvec[2].compare(old_chr)!=0){
                    	tio.dlwrite(out_dir+"/" + prefix+"_"+ old_chr+".sam",mpvec );
                    	old_chr= tokvec[2];
                    	mpvec.clear();
                    }
                    
                    mpvec.push_back(str);
                   
              }
        }
				
				tio.dlwrite(out_dir+"/" + prefix+"_"+ old_chr+".sam",mpvec );
        mpvec.clear();
        filestr.close();
	
}