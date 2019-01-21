#include <string>
#include <stdlib.h>
#include <stdio.h>
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
#include "Prep.h"
#include "Normalizer.h"
#include "Reg.h"
#include "CrossLinker.h"
#include "SNP.h"


using namespace std;


int main(int argc, char **argv) {
	
        TextIO tio;
        TOK tok;
        stringstream ss; 
	      string cmd="";
	      string CMDfname= argv[1];
	      map<string, double> SizeMap;
	      vector<string> cmdvec = tio.dlread_vector(CMDfname);
	      map<string,string> CMDMap ;
//	      
	      for(int i =0 ; i<cmdvec.size(); i++){
	      	vector<string> tokvec = tok.Tokenize( cmdvec[i],"=");
	      	CMDMap.insert(make_pair( tokvec[0], tokvec[1]));
	      	
	      }
	      
				
				vector<string> RepBamVec = tok.Tokenize( CMDMap["REPBAMFILES"],",");
				vector<string> RepBedVec =tok.Tokenize( CMDMap["REPBEDFILES"],",");
				vector<string> RepIDVec= tok.Tokenize( CMDMap["REPIDS"],",");
				
				vector<string> ChrmVec = tok.Tokenize( CMDMap["CHRMSTR"],",");
				
				string VCFFLAG = CMDMap["VCFFLAG"];
				string VCFfname = CMDMap["VCFNAME"];
////init tmpdir & outdir
//////	  
//				    
//	      cmd="rm -rf "+ CMDMap["TMPDIR"]+"/*input*";
//	      cout<< cmd<<endl;
//	      system(cmd.c_str());
	      
	      cmd="rm -rf "+ CMDMap["OUTDIR"]+"/pred*";
	      cout<< cmd<<endl;
	      system(cmd.c_str());
    
	      cmd="mkdir "+ CMDMap["TMPDIR"];
	      cout<< cmd<<endl;
	      system(cmd.c_str());
	      
	      cmd="mkdir "+ CMDMap["OUTDIR"];
	      cout<< cmd<<endl;
	      system(cmd.c_str());
//////	  		
////
////
//////1. read the lib sizes
//	//	get libsize
//////
	cmd = "samtools flagstat "+ CMDMap["INPUTBAMFILE"] + " > " +CMDMap["TMPDIR"]+"/input.stat" ;
  //system(cmd.c_str());
//
//
	vector<string> tmpvec = tio.dlread_vector( CMDMap["TMPDIR"]+"/input.stat" );
	string numstr;
	for(int i =0 ; i<tmpvec.size(); i++){
		vector<string> tokvec = tok.Tokenize( tmpvec[i], " ");
		int mapflag =0;
		int mateflag =0;
		for(int j =0 ; j<tokvec.size() ; j++){
			
			if(tokvec[j].compare("mapped")==0){
				mapflag =1;
			}
			if(tokvec[j].compare("mate")==0){
				mateflag =1;
			}
			if(mapflag ==1 && mateflag ==0){
				numstr = tmpvec[i];
			}
		}		
	}
	
	cout << numstr << endl; 
	
//	vector<string> tmptokvec = tok.Tokenize(numstr,"+");
//	//cout << atof(tmptokvec[0].c_str()) << endl;
//	SizeMap.insert(make_pair( "input", atof(tmptokvec[0].c_str())/1000000.0 ));
	//cout <<SizeMap["input"] << endl;

	for(int i =0; i<RepBamVec.size(); i++){
		
		cmd = "samtools flagstat "+ RepBamVec[i] + " > " +CMDMap["TMPDIR"]+"/"+RepIDVec[i]+".stat" ;
  	//system(cmd.c_str());
		
		vector<string> rtmpvec = tio.dlread_vector( CMDMap["TMPDIR"]+"/"+RepIDVec[i]+".stat" );
		
		for(int p =0 ; p<rtmpvec.size(); p++){
			vector<string> tokvec = tok.Tokenize( rtmpvec[p], " ");
			int mapflag =0;
			int mateflag =0;
			for(int q =0 ; q<tokvec.size() ; q++){
				if(tokvec[q].compare("mapped")==0){
					mapflag =1;
				}
				if(tokvec[q].compare("mate")==0){
					mateflag =1;
				}
				if(mapflag ==1 && mateflag ==0){
					numstr = rtmpvec[p];
				}
			}		
		}
		cout <<numstr << endl; 
		vector<string> rtmptokvec = tok.Tokenize(numstr,"+");
		//cout << atof(rtmptokvec[0].c_str()) << endl;
		SizeMap.insert(make_pair( RepIDVec[i] , atof(rtmptokvec[0].c_str())/1000000.0 ));
	}

////////
////////// 2. prep the input sample 	  		
////	
		cmd="./preProc "+ CMDMap["INPUTBEDFILE"] +" "+ CMDMap["INPUTBAMFILE"] + " "  +CMDMap["REFSEQ"] +" " + CMDMap["TMPDIR"] + " input";
		cout << cmd<< endl;
    Prep pp;
		//pp.sep_beds_by_chrs(CMDMap["INPUTBEDFILE"], CMDMap["TMPDIR"] ,"input");
		//pp.preProc(CMDMap["INPUTBAMFILE"], CMDMap["REFSEQ"] , CMDMap["TMPDIR"] ,"input" );
//////		
//////		//For Gio's Pipeline
////////		Prep pp;
////////		cmd="./preProc "+ CMDMap["INPUTBAMFILE"] + " "  +CMDMap["REFSEQ"] +" " + CMDMap["TMPDIR"] + " input";
////////		cout << cmd<< endl;
////////		pp.preProc_gio(CMDMap["INPUTBAMFILE"], CMDMap["REFSEQ"] , CMDMap["TMPDIR"] ,"input" );
////////		
////////		
//////////////
////////////// 3. prep all replicates 
//////
  //Prep pp;
//	for(int i =0; i<RepBamVec.size(); i++){
//		cmd="./preProc "+ RepBedVec[i] +" "+ RepBamVec[i] + " "  +CMDMap["REFSEQ"] +" " + CMDMap["TMPDIR"] + " "+RepIDVec[i];
//		cout << cmd<< endl;
//
//		pp.sep_beds_by_chrs(RepBedVec[i] , CMDMap["TMPDIR"] , RepIDVec[i] );
//		pp.preProc( RepBamVec[i] ,CMDMap["REFSEQ"] , CMDMap["TMPDIR"] , RepIDVec[i]);
//	}
////////////
////////////
////////////// 4. filter out PKs by fold-change rates.
//////////////	
//////////////////		
//	if(CMDMap["CALCFCR"].compare("T")==0){
//		for(int i =0; i<RepBamVec.size(); i++){
//		for(int j =0 ;  j< ChrmVec.size(); j++){
//			
//			
//			ss.str("");
//			ss<< "touch " +  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg_raw.bed" +" " + CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_2_"+ChrmVec[j]+".mp"+" "+CMDMap["TMPDIR"]+"/input_neg_2_"+ChrmVec[j]+".mp" ; 
//			cmd = ss.str();
//			system(cmd.c_str());
//			
//			ss.str("");
//			ss<< "filterPKs " +  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg_raw.bed" +" " + CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_2_"+ChrmVec[j]+".mp"+" "+CMDMap["TMPDIR"]+"/input_neg_2_"+ChrmVec[j]+".mp"+" "+CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg.bed"+" "<<SizeMap["input"] << " "<< SizeMap[RepIDVec[i]] << " "<<CMDMap["FILTERFOLDCHANGE"] ; 
//			cmd = ss.str();
//			cout << cmd<< endl;
//			
//			Normalizer nl ; 
//			nl.filterPKs(CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg_raw.bed" , CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_2_"+ChrmVec[j]+".mp",CMDMap["TMPDIR"]+"/input_neg_2_"+ChrmVec[j]+".mp" , CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg.bed" , SizeMap["input"] , SizeMap[RepIDVec[i]] ,atof(CMDMap["FILTERFOLDCHANGE"].c_str()));
//
//			
//			ss.str("");
//			ss<< "touch " +  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos_raw.bed" +" " + CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_2_"+ChrmVec[j]+".mp"+" "+CMDMap["TMPDIR"]+"/input_pos_2_"+ChrmVec[j]+".mp" ; 
//			cmd = ss.str();
//			system(cmd.c_str());
//				
//
//			ss.str("");
//			ss<< "filterPKs " +  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos_raw.bed" +" " + CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_2_"+ChrmVec[j]+".mp"+" "+CMDMap["TMPDIR"]+"/input_pos_2_"+ChrmVec[j]+".mp"+" "+CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos.bed"+" "<<SizeMap["input"] << " "<< SizeMap[RepIDVec[i]] << " "<<CMDMap["FILTERFOLDCHANGE"] ; 
//			cmd = ss.str();
//			cout << cmd<< endl;
//			nl.filterPKs(CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos_raw.bed" , CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_2_"+ChrmVec[j]+".mp",CMDMap["TMPDIR"]+"/input_pos_2_"+ChrmVec[j]+".mp",CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos.bed", SizeMap["input"] , SizeMap[RepIDVec[i]] ,atof(CMDMap["FILTERFOLDCHANGE"].c_str()));
//			
//						
//		}
//		}
//	}
//	else{
//		for(int i =0; i<RepBamVec.size() ; i++){
//			for(int j =0 ;  j< ChrmVec.size(); j++){
//			
//				ss.str("");
//				ss<< "touch " +  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg_raw.bed";
//				cmd = ss.str();
//				system(cmd.c_str());
//			
//				ss.str("");
//				ss<< "filterPKs " +  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg_raw.bed" +" "+CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg.bed"+ " "<< SizeMap[RepIDVec[i]] << " "<<CMDMap["FILTERFOLDCHANGE"] ; 
//				cmd = ss.str();
//				cout << cmd<< endl;
//			
//				Normalizer nl ; 
//				nl.filterPKs(CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg_raw.bed" , CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg.bed" , log2( atof(CMDMap["FILTERFOLDCHANGE"].c_str()) ) );
//
//			
//				ss.str("");
//				ss<< "touch " +  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos_raw.bed"; 
//				cmd = ss.str();
//				system(cmd.c_str());
//				
//
//				ss.str("");
//				ss<< "filterPKs " +  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos_raw.bed" +" "+CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos.bed"+" "<<CMDMap["FILTERFOLDCHANGE"] ; 
//				cmd = ss.str();
//				cout << cmd<< endl;
//				nl.filterPKs(CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos_raw.bed" , CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos.bed", log2(atof(CMDMap["FILTERFOLDCHANGE"].c_str()) ) );
//			
//						
//			}
//		}
//	}
	
////////////	
//////////
/////////////// 5. estimate seq error rates 
////////////
////////
//	double fdr_cutoff = atof(CMDMap["SEQERRFDR"].c_str());
//	vector<string> outcutvec; 
//	
//	for(int i =0; i<10; i++){
//		outcutvec.push_back("");
//	}
//	for(int i =0; i<RepBamVec.size(); i++){
//		map<int, vector<double> > SeqErrDistMap;
//
//		for(int p =0 ; p< 10; p++){
//			vector<double> tmpvec ;
//			for(int q =0 ; q<20; q++){
//				tmpvec.push_back(0.0); 
//			}
//			SeqErrDistMap.insert(make_pair( p, tmpvec ));
//		}
//		
//		for(int j =0 ;  j< ChrmVec.size(); j++){
//			
//			ss.str("");
//			ss<< "estimateSeqError " +  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_"+ChrmVec[j]+".mp" ;
//			cmd = ss.str();
//			cout << cmd<< endl;
//   		map<int, vector<double> > tmpdistmap;	
//			Normalizer nl ; 
//			 //tmpdistmap = nl.estimateSeqError( CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_"+ChrmVec[j]+".mp", CMDMap["SNPDIR"]+"/"+ChrmVec[j]+".txt", SizeMap[RepIDVec[i]],1);
//			 //cout << CMDMap["SNPDIR"]+"/"+ChrmVec[j]+".txt" << endl;
//			 tmpdistmap = nl.estimateSeqError( CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_"+ChrmVec[j]+".mp", CMDMap["SNPDIR"]+"/"+ChrmVec[j]+".txt", SizeMap[RepIDVec[i]],1);
//			 for(int p =0 ; p< 10; p++){
//				
//				for(int q =0 ; q<20; q++){
//					 SeqErrDistMap[p][q] = SeqErrDistMap[p][q] + tmpdistmap[p][q]; 
//				}
//				
//			 }
//			 
//			 
//			ss.str("");
//			ss<< "estimateSeqError " +  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_"+ChrmVec[j]+".mp" ;
//			cmd = ss.str();
//			cout << cmd<< endl;
//			
//			 tmpdistmap = nl.estimateSeqError( CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_"+ChrmVec[j]+".mp", CMDMap["SNPDIR"]+"/"+ChrmVec[j]+".txt",SizeMap[RepIDVec[i]],0);
//			 for(int p=0 ; p< 10; p++){
//				
//				for(int q=0 ; q<20; q++){
//					 SeqErrDistMap[p][q] = SeqErrDistMap[p][q] + tmpdistmap[p][q]; 
//				}
//				
//			 }
//			 			
//		}
//		
//		
//    
//    vector<string> outrawvec;
//		for(int p =0 ; p< 10; p++){
//			ss.str("");
//			for(int q =0 ; q<20; q++){
//				ss << SeqErrDistMap[p][q];
//				if(q<19){
//					ss<<"\t";
//				}
//			}
//			outrawvec.push_back(ss.str());
//		}
//		tio.dlwrite( CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_seqerr.raw", outrawvec);
//    
//		
//		
//		for(int p =0 ; p< 10; p++){
//			double row_sum =0;
//			for(int q=0 ; q<20; q++){
//				 row_sum += SeqErrDistMap[p][q];
//			}
//			
//			for(int q=0 ; q<20; q++){
//				 SeqErrDistMap[p][q] = SeqErrDistMap[p][q] /row_sum ;
//			}
//		}
//		
//		
//		//output seq err pmf
//		vector<string> outpmfvec;
//		for(int p =0 ; p< 10; p++){
//			ss.str("");
//			for(int q =0 ; q<20; q++){
//				ss << SeqErrDistMap[p][q];
//				if(q<19){
//					ss<<"\t";
//				}
//			}
//			outpmfvec.push_back(ss.str());
//		}
//		tio.dlwrite( CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_seqerr.pmf", outpmfvec);
//		
//		vector<string> outcdfvec;
//		
//		for(int p =0 ; p< 10; p++){
//			double sum =0 ; 
//			ss.str("");
//			
//			for(int q =0 ; q<20; q++){
//				sum+=SeqErrDistMap[p][q];
//				ss << sum;
//				if(q<19){
//					ss<<"\t";
//				}
//			}
//			outcdfvec.push_back(ss.str());
//		}
//		tio.dlwrite( CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_seqerr.cdf", outcdfvec);
//		
//		//write the cut-off vector
//		for(int p =0 ; p< 10; p++){
//			double sum =0 ; 
//			
//			int cutidx = 0;
//			for(int q =0 ; q<20; q++){
//				sum+=SeqErrDistMap[p][q];
//				if(sum>0.9){
//					cutidx=q;
//					
//					break;
//				}
//			}
//			ss.str("");
//			ss<< (cutidx +1)*0.05;
//			if(i<RepBamVec.size()-1){
//				ss<<"\t";
//			}
//			outcutvec[p]=outcutvec[p]+ss.str();
//		}
//		
//	}
//	
//	tio.dlwrite( CMDMap["TMPDIR"]+"/seqerr.cut", outcutvec);
//	
//////	
//////////	
////////
//////////	
//////////////	
////////////
//////////////	
//////////////	
////////////////////6. predict CLs in the input sample
//////////////////
////////
////////	for(int j =0 ;  j< ChrmVec.size(); j++){
////////			
////////			cmd = "touch " + CMDMap["TMPDIR"]+"/input_"+ChrmVec[j]+"_pos_raw.bed " + CMDMap["TMPDIR"]+"/input_pos_2_"+ChrmVec[j]+".mp";
////////			system(cmd.c_str());
////////					
////////			ss.str("");
////////			ss<< "./predCLs " <<  CMDMap["TMPDIR"]+"/input_"+ChrmVec[j]+"_pos_raw.bed" << " "<<CMDMap["TMPDIR"]+"/input_pos_2_"+ChrmVec[j]+".mp" << " " << CMDMap["TMPDIR"]+"/input_pos_2_"+ChrmVec[j]+".cls" << " + "<< CMDMap["CLSFDR"];
////////			cmd = ss.str();
////////			cout<<cmd<<endl;
//////////			system(cmd.c_str());
////////			CrossLinker cl ; 
////////			cl.predCLs(CMDMap["TMPDIR"]+"/input_"+ChrmVec[j]+"_pos_raw.bed", CMDMap["TMPDIR"]+"/input_pos_2_"+ChrmVec[j]+".mp" , CMDMap["TMPDIR"]+"/input_pos_2_"+ChrmVec[j]+".cls", "+" , CMDMap["CLSFDR"]  );							
////////			
////////			
////////			cmd = "touch " + CMDMap["TMPDIR"]+"/input_"+ChrmVec[j]+"_neg_raw.bed " + CMDMap["TMPDIR"]+"/input_neg_2_"+ChrmVec[j]+".mp";
////////			system(cmd.c_str());
////////			
////////			ss.str("");
////////			ss<< "./predCLs " <<  CMDMap["TMPDIR"]+"/input_"+ChrmVec[j]+"_neg_raw.bed" << " "<<CMDMap["TMPDIR"]+"/input_neg_2_"+ChrmVec[j]+".mp" << " " << CMDMap["TMPDIR"]+"/input_neg_2_"+ChrmVec[j]+".cls" << " - "<< CMDMap["CLSFDR"];
////////			cmd = ss.str();
////////			cout<<cmd<<endl;
////////			cl.predCLs(CMDMap["TMPDIR"]+"/input_"+ChrmVec[j]+"_neg_raw.bed", CMDMap["TMPDIR"]+"/input_neg_2_"+ChrmVec[j]+".mp" , CMDMap["TMPDIR"]+"/input_neg_2_"+ChrmVec[j]+".cls", "-" , CMDMap["CLSFDR"]  );							
////////		
////////	}
//////
////////////	
//////////////////
////////////////  7. predict CLs in the Rep samples
////
//	for(int i=0 ; i<RepIDVec.size();i++ ){
//		for(int j =0 ;  j< ChrmVec.size(); j++){
//			
//			cmd = "touch " + CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos.bed " + CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_2_"+ChrmVec[j]+".mp";
//			system(cmd.c_str());
//					
//			ss.str("");
//			ss<< "./predCLs " <<  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos.bed" << " "<<CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_2_"+ChrmVec[j]+".mp" << " " << CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_2_"+ChrmVec[j]+".cls" << " + "<< CMDMap["CLSFDR"];
//			cmd = ss.str();
//			cout<<cmd<<endl;
//
//			CrossLinker cl ; 
//			cl.predCLs(CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos.bed", CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_2_"+ChrmVec[j]+".mp" , CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_2_"+ChrmVec[j]+".cls", "+" , CMDMap["CLSFDR"]  );							
//
//			
//			
//			cmd = "touch " + CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg.bed " + CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_2_"+ChrmVec[j]+".mp";
//			system(cmd.c_str());
//			
//			ss.str("");
//			ss<< "./predCLs " <<  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg.bed" << " "<<CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_2_"+ChrmVec[j]+".mp" << " " << CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_2_"+ChrmVec[j]+".cls" << " - "<< CMDMap["CLSFDR"];
//			cmd = ss.str();
//			cout<<cmd<<endl;
//			//system(cmd.c_str());
//			cl.predCLs(CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg.bed", CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_2_"+ChrmVec[j]+".mp" , CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_2_"+ChrmVec[j]+".cls", "-" , CMDMap["CLSFDR"]  );							
//		}
//	}
////////	
//////////	
////////	
//////////////
////////////////// 8. calc Prior of Seq
////////////
////////////////	//get refseq path
//////
//////	string seqfname = CMDMap["REFSEQ"] ;
//////	
//////	size_t found;
//////
//////  found=seqfname.find_last_of("/\\");
//////	string seqpath = seqfname.substr(0,found);
//////		
//////	cmd = "calcPrior " +CMDMap["TMPDIR"] + " "+ seqpath +" " + CMDMap["TMPDIR"] +" "+CMDMap["CHRMSTR"] + " "+ CMDMap["PRIORWINDOW"];
//////	cout << cmd << endl;
//////	CrossLinker cl; 
//////	cl.prior(CMDMap["TMPDIR"] , seqpath ,  CMDMap["TMPDIR"] , CMDMap["CHRMSTR"] , CMDMap["PRIORWINDOW"]);
//////
////////////
////////////
////////////// 9(a).  for each rep sample, calc the RCs of SNPs
//////
  if(VCFFLAG.compare("T") ==0 ){
  	
  	//for traning 
  	
  	SNP snp;
    snp.SNPtoRC_train(CMDMap["SNPDIR"] , CMDMap["TMPDIR"] , RepIDVec , CMDMap["CHRMSTR"], SizeMap);
		
		cout << "training ..... done"<< endl;
		
		//for testing
		
		//snp.SNPtoRC_test(VCFfname , CMDMap["TMPDIR"] , RepIDVec  , CMDMap["CHRMSTR"], SizeMap);
		
		
  }
	else{
		
// old version of de novo prediction 	
//		for(int i=0 ; i<RepIDVec.size();i++){
//			ss.str("");
//			ss<< "SNPtoRC " + CMDMap["SNPDIR"] + " "+ CMDMap["TMPDIR"]+" "+ RepIDVec[i]+" "+CMDMap["CHRMSTR"]+ " "<<SizeMap[RepIDVec[i]]; 
//			cmd= ss.str();
//			cout << cmd<< endl;
//		
//			SNP snp;
//			snp.SNPtoRC(CMDMap["SNPDIR"] , CMDMap["TMPDIR"] , RepIDVec[i] , CMDMap["CHRMSTR"], SizeMap[RepIDVec[i]]);
////		
//		}

// new version of de novo prediction 
//	
//			ss<< "SNPtoRC " + CMDMap["SNPDIR"] + " "+ CMDMap["TMPDIR"]+" "+CMDMap["CHRMSTR"]; 
//			cmd= ss.str();
//			cout << cmd<< endl;
////		
//			SNP snp;
//			snp.SNPtoRC_denovo(CMDMap["SNPDIR"] , CMDMap["TMPDIR"] , RepIDVec , CMDMap["CHRMSTR"], SizeMap);
//			
		
		
	}
////	
////////////////	
//////////
//////////
//////////
////////////////
////////////////////10.  merge RCs and predicts the ASBs
////////
//	cmd="mergeRCs&Pred " + CMDMap["TMPDIR"] + " "+ CMDMap["REPIDS"];
////	
//	cout << cmd<< endl;
////	//system(cmd.c_str());
////
////
//	Reg reg;
//	if(VCFFLAG.compare("T") ==0 ){
//		cmd="Rscript Reg_vcf.r "+ CMDMap["TMPDIR"]+"/reg_vcf.rc "+ CMDMap["TMPDIR"]+"/test_vcf.rc "+CMDMap["TMPDIR"] +"/pred_vcf.rc 2 ";
//		cout << cmd<< endl;
//	  system(cmd.c_str());;
//	}
//	else{
//		//old version
//		//reg.mergeRCs(CMDMap["REPIDS"] , CMDMap["TMPDIR"] , CMDMap["RBPID"], atof(CMDMap["MINRC"].c_str()) );
//		
//		//new vesion denovo 1/15
//		cmd="Rscript Reg_denovo.r "+ CMDMap["TMPDIR"]+"/reg_denovo.rc "+ CMDMap["TMPDIR"]+"/test_denovo.rc "+CMDMap["TMPDIR"] +"/pred_denovo.rc 2 ";
//		cout << cmd<< endl;
//	  system(cmd.c_str());
//	}
////
////
//////		
////	11. calc FDR values and list ASB SNPs
//
//	Reg reg;
//
//	string dir_path(argv[1]) ; 
//	cout << "open  "<< dir_path<< endl; 
//	
//	reg.calcFDR(dir_path );
//	
	// here run Rscript FDR.r the predict file to convert p-values to the FDR values
////	

	

	  
	return 0;
}