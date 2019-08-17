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
     	double PK_rc = atof(CMDMap["FILTERFOLDCHANGE"].c_str());
////init tmpdir & outdir

	   	   
	    cmd="rm -rf "+ CMDMap["OUTDIR"]+"/pred*";
	    cout<< cmd<<endl;
	    system(cmd.c_str());
    
    	cmd="rm -rf "+ CMDMap["TMPDIR"];
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
	cout<< cmd<<endl;
    system(cmd.c_str());
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
		cout<< cmd<<endl;
  	    system(cmd.c_str());
		
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
		
		vector<string> rtmptokvec = tok.Tokenize(numstr,"+");
		cout << atof(rtmptokvec[0].c_str()) << endl;
		SizeMap.insert(make_pair( RepIDVec[i] , atof(rtmptokvec[0].c_str())/1000000.0 ));
	}

////////
////////// 2. prep the input sample 	  		
////	
		cmd="./preProc "+ CMDMap["INPUTBEDFILE"] +" "+ CMDMap["INPUTBAMFILE"] + " "  +CMDMap["REFSEQ"] +" " + CMDMap["TMPDIR"] + " input";
		cout << cmd<< endl;
    	Prep pp;
		pp.sep_beds_by_chrs(CMDMap["INPUTBEDFILE"], CMDMap["TMPDIR"] ,"input");
		pp.preProc(CMDMap["INPUTBAMFILE"], CMDMap["REFSEQ"] , CMDMap["TMPDIR"] ,"input" );
//////		

///////// 3. prep all replicates 
//////
  //Prep pp;
	for(int i =0; i<RepBamVec.size(); i++){
		cmd="./preProc "+ RepBedVec[i] +" "+ RepBamVec[i] + " "  +CMDMap["REFSEQ"] +" " + CMDMap["TMPDIR"] + " "+RepIDVec[i];
		cout << cmd<< endl;

		pp.sep_beds_by_chrs(RepBedVec[i] , CMDMap["TMPDIR"] , RepIDVec[i], log2(PK_rc) );
		pp.preProc( RepBamVec[i] ,CMDMap["REFSEQ"] , CMDMap["TMPDIR"] , RepIDVec[i]);
	}

////////////

	
////////////////////4. predict CLs in the input sample
//////////////////

	for(int j =0 ;  j< ChrmVec.size(); j++){
			
			cmd = "touch " + CMDMap["TMPDIR"]+"/input_"+ChrmVec[j]+"_pos.bed " + CMDMap["TMPDIR"]+"/input_pos_2_"+ChrmVec[j]+".mp";
			system(cmd.c_str());
					
			ss.str("");
			ss<< "./predCLs " <<  CMDMap["TMPDIR"]+"/input_"+ChrmVec[j]+"_pos.bed" << " "<<CMDMap["TMPDIR"]+"/input_pos_2_"+ChrmVec[j]+".mp" << " " << CMDMap["TMPDIR"]+"/input_pos_2_"+ChrmVec[j]+".cls" << " + "<< CMDMap["CLSFDR"];
			cmd = ss.str();
			cout<<cmd<<endl;
//			system(cmd.c_str());
			CrossLinker cl ; 
			cl.predCLs(CMDMap["TMPDIR"]+"/input_"+ChrmVec[j]+"_pos.bed", CMDMap["TMPDIR"]+"/input_pos_2_"+ChrmVec[j]+".mp" , CMDMap["TMPDIR"]+"/input_pos_2_"+ChrmVec[j]+".cls", "+" , CMDMap["CLSFDR"]  );							
			
			
			cmd = "touch " + CMDMap["TMPDIR"]+"/input_"+ChrmVec[j]+"_neg.bed " + CMDMap["TMPDIR"]+"/input_neg_2_"+ChrmVec[j]+".mp";
			system(cmd.c_str());
			
			ss.str("");
			ss<< "./predCLs " <<  CMDMap["TMPDIR"]+"/input_"+ChrmVec[j]+"_neg.bed" << " "<<CMDMap["TMPDIR"]+"/input_neg_2_"+ChrmVec[j]+".mp" << " " << CMDMap["TMPDIR"]+"/input_neg_2_"+ChrmVec[j]+".cls" << " - "<< CMDMap["CLSFDR"];
			cmd = ss.str();
			cout<<cmd<<endl;
			cl.predCLs(CMDMap["TMPDIR"]+"/input_"+ChrmVec[j]+"_neg.bed", CMDMap["TMPDIR"]+"/input_neg_2_"+ChrmVec[j]+".mp" , CMDMap["TMPDIR"]+"/input_neg_2_"+ChrmVec[j]+".cls", "-" , CMDMap["CLSFDR"]  );							
		
	}
//////
////////////	
//////////////////
////////////////  5. predict CLs in the Rep samples
// ////
	for(int i=0 ; i<RepIDVec.size();i++ ){
		for(int j =0 ;  j< ChrmVec.size(); j++){
			
			cmd = "touch " + CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos.bed " + CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_2_"+ChrmVec[j]+".mp";
			system(cmd.c_str());
					
			ss.str("");
			ss<< "./predCLs " <<  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos.bed" << " "<<CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_2_"+ChrmVec[j]+".mp" << " " << CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_2_"+ChrmVec[j]+".cls" << " + "<< CMDMap["CLSFDR"];
			cmd = ss.str();
			cout<<cmd<<endl;

			CrossLinker cl ; 
			cl.predCLs(CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_pos.bed", CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_2_"+ChrmVec[j]+".mp" , CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_pos_2_"+ChrmVec[j]+".cls", "+" , CMDMap["CLSFDR"]  );							

			
			
			cmd = "touch " + CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg.bed " + CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_2_"+ChrmVec[j]+".mp";
			system(cmd.c_str());
			
			ss.str("");
			ss<< "./predCLs " <<  CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg.bed" << " "<<CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_2_"+ChrmVec[j]+".mp" << " " << CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_2_"+ChrmVec[j]+".cls" << " - "<< CMDMap["CLSFDR"];
			cmd = ss.str();
			cout<<cmd<<endl;
			//system(cmd.c_str());
			cl.predCLs(CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_"+ChrmVec[j]+"_neg.bed", CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_2_"+ChrmVec[j]+".mp" , CMDMap["TMPDIR"]+"/"+RepIDVec[i]+"_neg_2_"+ChrmVec[j]+".cls", "-" , CMDMap["CLSFDR"]  );							
		}
	}
////////	
//////////	
////////	
//////////////
////////////////// 6. calc Prior of Seq
////////////
////////////////	//get refseq path
// //////
	string seqfname = CMDMap["REFSEQ"] ;
	
	size_t found;

 	found=seqfname.find_last_of("/\\");
	string seqpath = seqfname.substr(0,found);
		
	cmd = "calcPrior " +CMDMap["TMPDIR"] + " "+ seqpath +" " + CMDMap["TMPDIR"] +" "+CMDMap["CHRMSTR"] + " "+ CMDMap["PRIORWINDOW"];
	cout << cmd << endl;
	CrossLinker cl; 
	cl.prior(CMDMap["TMPDIR"] , seqpath ,  CMDMap["TMPDIR"] , CMDMap["CHRMSTR"] , CMDMap["PRIORWINDOW"]);

// ////////////
// ////////////
// ////////////// 7.  predic candidates of ASB SNPs
// //////
// 
  	
//   	//for traning 
  	
  	SNP snp;
    cout << "training ..... start"<< endl;
    snp.SNPtoRC_train(CMDMap["SNPDIR"] , CMDMap["TMPDIR"] , RepIDVec , CMDMap["CHRMSTR"], SizeMap);
		
	cout << "training ..... done"<< endl;
		
		//for testing
	cout << "testing ..... start"<< endl;	
	snp.SNPtoRC_test(VCFfname , CMDMap["TMPDIR"] , RepIDVec  , CMDMap["CHRMSTR"], SizeMap);
	cout << "testing ..... done"<< endl;
		

	ss.str("");
	ss<<"Rscript Reg_denovo.r "+ CMDMap["TMPDIR"]+"/reg_denovo.rc "+ CMDMap["TMPDIR"]+"/test_vcf.rc "+CMDMap["OUTDIR"] +"/asb_candidates.txt "<< RepIDVec.size(); 
	cmd=ss.str();
	cout << cmd<< endl;	
        system(cmd.c_str());	
		

//	
	// here run Rscript FDR.r the predict file to convert p-values to the FDR values
////	

	

	  
	return 0;
}
