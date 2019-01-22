#include "Sim.h"

#include <string>
#include <stdlib.h>
#include <stdio.h>


using namespace std;

int main(int argc, char **argv) {
	
	//Sim sm ; 
	TextIO tio;
	TOK tok;
	stringstream ss; 	
	int selected =0;
	double AR =atof(argv[1]);

	vector<string>  auc_statvec ;
	auc_statvec.push_back("ID\tBEAPR_AUC\tBEAPR_SEN\tBEAPR_SPE\tCHI_AUC\tCHI_SEN\tCHI_SPE\tFET_AUC\tFET_SEN\tFET_SPE");

	for(int q =0 ; q<500; q++){
	stringstream ss; 
	string cmd; 
	ss.str("");
	ss<<"./genMean 10 " << AR;
	cmd = ss.str();
	system(cmd.c_str());

	ss.str("");
	ss<<"Rscript Sim.r /home/yyang/data/eclip/K562/SRSF1_test/sim/mean 2 10";
	cmd = ss.str();
	system(cmd.c_str());

	ss.str("");
	ss<<"Rscript Reg.r /home/yyang/data/eclip/K562/SRSF1_test/tmp/reg.rc /home/yyang/data/eclip/K562/SRSF1_test/sim/mean.out /home/yyang/data/eclip/K562/SRSF1_test/sim/asbclip.in";
	cmd = ss.str();
	system(cmd.c_str());	


	ss.str("");
	ss<<"Rscript Pred.r /home/yyang/data/eclip/K562/SRSF1_test/sim/asbclip.in /home/yyang/data/eclip/K562/SRSF1_test/sim/asbclip.out 2";
	cmd = ss.str();
	system(cmd.c_str());

	ss.str("");
	ss<<"Rscript chi.r /home/yyang/data/eclip/K562/SRSF1_test/sim/mean.out /home/yyang/data/eclip/K562/SRSF1_test/sim/chi.out";
	cmd = ss.str();
	system(cmd.c_str());

		ss.str("");
		ss<<"Rscript fisher.r /home/yyang/data/eclip/K562/SRSF1_test/sim/mean.out /home/yyang/data/eclip/K562/SRSF1_test/sim/fisher.out";
		cmd = ss.str();
		system(cmd.c_str());

		 ss.str("");
		 ss<<"Rscript pROC_new.r /home/yyang/data/eclip/K562/SRSF1_test/sim/asbclip.out";
		 cmd = ss.str();
		 system(cmd.c_str());

		 ss.str("");
		 ss<<"Rscript pROC_new.r /home/yyang/data/eclip/K562/SRSF1_test/sim/chi.out";
		 cmd = ss.str();
		 system(cmd.c_str());

		 ss.str("");
		 ss<<"Rscript pROC_new.r /home/yyang/data/eclip/K562/SRSF1_test/sim/fisher.out";
		 cmd = ss.str();
		 system(cmd.c_str());

		 ss.str("");
		 ss<<"Rscript PRROC.r /home/yyang/data/eclip/K562/SRSF1_test/sim/asbclip.out  /home/yyang/data/eclip/K562/SRSF1_test/sim/chi.out /home/yyang/data/eclip/K562/SRSF1_test/sim/fisher.out";
		 cmd = ss.str();
                 system(cmd.c_str());		

		
		

	//
    vector<string> PRROCvec =tio.dlread_vector("/home/yyang/data/eclip/K562/SRSF1_test/sim/asbclip.out.pr"); 
    
    double asb_auc = atof(PRROCvec[1].c_str());
    double chi_auc = atof(PRROCvec[2].c_str());
    double fet_auc = atof(PRROCvec[3].c_str());

    vector<string> pROCvec =tio.dlread_vector("/home/yyang/data/eclip/K562/SRSF1_test/sim/asbclip.out.matrices"); 
    
    vector<string> tmpvec =  tok.Tokenize(pROCvec[1] , "\t");;
    double asb_sen = atof(tmpvec[1].c_str());
    double asb_spe = atof(tmpvec[2].c_str());

    pROCvec =tio.dlread_vector("/home/yyang/data/eclip/K562/SRSF1_test/sim/chi.out.matrices"); 
  
    tmpvec =   tok.Tokenize(pROCvec[1] , "\t");;
    double chi_sen = atof(tmpvec[1].c_str());
    double chi_spe = atof(tmpvec[2].c_str());

    pROCvec =tio.dlread_vector("/home/yyang/data/eclip/K562/SRSF1_test/sim/fisher.out.matrices"); 
  
    tmpvec =   tok.Tokenize(pROCvec[1] , "\t");;
    double fet_sen = atof(tmpvec[1].c_str());
    double fet_spe = atof(tmpvec[2].c_str());



    cout <<"================"<< asb_sen - fet_sen << endl;
    if(asb_auc - fet_auc > 0.03 ){
    	ss.str("") ;
		ss<< q;    	
    	ss<< "\t"<< asb_auc<<"\t";
    	ss<< "\t"<< asb_sen<<"\t";
    	ss<< "\t"<< asb_spe<<"\t";
    	ss<< "\t"<< chi_auc<<"\t";
    	ss<< "\t"<< chi_sen<<"\t";
    	ss<< "\t"<< chi_spe<<"\t";
    	ss<< "\t"<< fet_auc<<"\t";
    	ss<< "\t"<< fet_sen<<"\t";
    	ss<< "\t"<< fet_spe<<"\t";

//     	string auc_str = ss.str(); 

//         vector<string> outvec ;

//   vector<string> asb_predvec = tio.dlread_vector("/home/yyang/data/eclip/K562/SRSF1_test/sim/asbclip.out");;
//   vector<string> chi_predvec = tio.dlread_vector("/home/yyang/data/eclip/K562/SRSF1_test/sim/chi.out");;
//   vector<string> fet_predvec = tio.dlread_vector("/home/yyang/data/eclip/K562/SRSF1_test/sim/fisher.out");; 
 
//  double fdr_cut = 0.1;

//  int asb_tp,chi_tp,fet_tp;
//  asb_tp=chi_tp=fet_tp = 0 ;
//  int asb_fp,chi_fp,fet_fp;
//  asb_fp=chi_fp=fet_fp = 0 ;
//  int asb_tn,chi_tn,fet_tn;
//  asb_tn=chi_tn=fet_tn = 0 ;
//  int asb_fn,chi_fn,fet_fn;
//  asb_fn=chi_fn=fet_fn = 0 ;


//   // vector<string> chivec ; 
//   // chivec.push_back("label\tfdr");
//   // vector<string> asbvec ;
//   // asbvec.push_back("label\tfdr");
//   // vector<string> fetvec ;
//   // fetvec.push_back("label\tfdr");




//  for(int i =0; i<asb_predvec.size(); i++){

//   	vector<string> chivec  = tok.Tokenize(chi_predvec[i] , "\t"); 
//   	string label = chivec[3];

  	
    
//   	vector<string> asbvec = tok.Tokenize(asb_predvec[i] , "\t"); 
//     //asbvec.push_back(label +"\t"+tokvec[tokvec.size()-1] );

//     vector<string> fetvec = tok.Tokenize(fet_predvec[i] , "\t"); 
// 	//fetvec.push_back(label +"\t"+tokvec[tokvec.size()-1]);

// 	double asb_fdr = atof(asbvec[asbvec.size()-1].c_str()); 
// 	double chi_fdr = atof(chivec[chivec.size()-1].c_str());
// 	double fet_fdr = atof(fetvec[fetvec.size()-1].c_str());

// 	if(chivec[3].compare("0") == 0){

// 		if(asb_fdr< fdr_cut) {
// 			asb_fp ++;
// 		}
// 		else{
// 			asb_tn++;

// 		}

// 		if(chi_fdr < fdr_cut){
// 			chi_fp++;
// 		} 
// 		else{
// 			chi_tn++;
// 		}
// 		if(fet_fdr < fdr_cut){
// 			fet_fp++;
// 		} 
// 		else{
// 			fet_tn++;	
// 		}

// 	}
// 	else{
// 		if(asb_fdr< fdr_cut) {
// 			asb_tp ++;
// 		}
// 		else{
// 			asb_fn++;
// 		}

// 		if(chi_fdr < fdr_cut){
// 			chi_tp++;
// 		} 
// 		else{
// 			chi_fn++;
// 		}
// 		if(fet_fdr < fdr_cut){
// 			fet_tp++;
// 		} 
// 		else{
// 			fet_fn++;
// 		}
// 	}

//    }







//  int asb_tp5,chi_tp5,fet_tp5;
//  asb_tp5=chi_tp5=fet_tp5 = 0 ;
//  int asb_fp5,chi_fp5,fet_fp5;
//  asb_fp5=chi_fp5=fet_fp5 = 0 ;
//  int asb_tn5,chi_tn5,fet_tn5;
//  asb_tn5=chi_tn5=fet_tn5 = 0 ;
//  int asb_fn5,chi_fn5,fet_fn5;
//  asb_fn5=chi_fn5=fet_fn5 = 0 ;
  
//   for(int i =0; i<asb_predvec.size(); i++){

//   	vector<string> chivec  = tok.Tokenize(chi_predvec[i] , "\t"); 
//   	string label = chivec[3];

  	
    
//   	vector<string> asbvec = tok.Tokenize(asb_predvec[i] , "\t"); 
//     //asbvec.push_back(label +"\t"+tokvec[tokvec.size()-1] );

//     vector<string> fetvec = tok.Tokenize(fet_predvec[i] , "\t"); 
// 	//fetvec.push_back(label +"\t"+tokvec[tokvec.size()-1]);

// 	double asb_fdr = atof(asbvec[asbvec.size()-1].c_str()); 
// 	double chi_fdr = atof(chivec[chivec.size()-1].c_str());
// 	double fet_fdr = atof(fetvec[fetvec.size()-1].c_str());

// 	if(chivec[3].compare("0") == 0){

// 		if(asb_fdr< 0.05) {
// 			asb_fp5 ++;
// 		}
// 		else{
// 			asb_tn5++;

// 		}

// 		if(chi_fdr < 0.05){
// 			chi_fp5++;
// 		} 
// 		else{
// 			chi_tn5++;
// 		}
// 		if(fet_fdr < 0.05){
// 			fet_fp5++;
// 		} 
// 		else{
// 			fet_tn5++;	
// 		}

// 	}
// 	else{
// 		if(asb_fdr< 0.05) {
// 			asb_tp5 ++;
// 		}
// 		else{
// 			asb_fn5++;
// 		}

// 		if(chi_fdr < 0.05){
// 			chi_tp5++;
// 		} 
// 		else{
// 			chi_fn5++;
// 		}
// 		if(fet_fdr < 0.05){
// 			fet_tp5++;
// 		} 
// 		else{
// 			fet_fn5++;
// 		}
// 	}

//    } 
// // // tio.dlwrite( "/home/yyang/projects/eCLIP/motif/Fig.2/Sim/result/08222018/0.8/chi.out" ,  chivec );
// // // tio.dlwrite( "/home/yyang/projects/eCLIP/motif/Fig.2/Sim/result/08222018/0.8/asb.out" ,  asbvec );
// // // tio.dlwrite( "/home/yyang/projects/eCLIP/motif/Fig.2/Sim/result/08222018/0.8/fet.out" ,  fetvec );
// //   double pre,rec ;


//    double asb_pre = asb_tp /(double)(asb_tp+asb_fp); 
//    double asb_rec = asb_tp /(double)(asb_tp+asb_fn);
//    double asb_fpr = asb_fp /(double)(asb_tn+asb_fp);
//    double asb_f = 2*(asb_pre*asb_rec) /(double)(asb_pre+asb_rec);
//    ss<< asb_pre<< "\t"<< asb_rec<< "\t"<< asb_fpr<< "\t"<< asb_f <<"\t";
   
//    double chi_pre = chi_tp /(double)(chi_tp+chi_fp); 
//    double chi_rec = chi_tp /(double)(chi_tp+chi_fn);
//    double chi_fpr = chi_fp /(double)(chi_tn+chi_fp);
//    double chi_f = 2*(chi_pre*chi_rec) /(double)(chi_pre+chi_rec);
//    ss<< chi_pre<< "\t"<< chi_rec<< "\t"<< chi_fpr<< "\t"<< chi_f <<"\t";
   
//    double fet_pre = fet_tp /(double)(fet_tp+fet_fp); 
//    double fet_rec = fet_tp /(double)(fet_tp+fet_fn);
//    double fet_fpr = fet_fp /(double)(fet_tn+fet_fp);
//    double fet_f = 2*(fet_pre*fet_rec) /(double)(fet_pre+fet_rec);
//     ss<< fet_pre<< "\t"<< fet_rec<< "\t"<< fet_fpr<< "\t"<< fet_f <<"\t";
	
    

//    	cout << "asb_pre : "<< asb_tp /(double)(asb_tp+asb_fp) << endl;

//    		cout << "asb rec : "<< asb_tp /(double)(asb_tp+asb_fn) << endl;
//    				cout << "asb fpr : "<< asb_fp /(double)(asb_tn+asb_fp) << endl;

//    				cout << "chi pre : "<< chi_tp /(double)(chi_tp+chi_fp) << endl;
//    				cout << "chi rec : "<< chi_tp /(double)(chi_tp+chi_fn) << endl;
//    				cout << "chi fpr : "<< chi_fp /(double)(chi_tn+chi_fn) << endl;
 
//    				cout << "fet pre : "<< fet_tp /(double)(fet_tp+fet_fp) << endl;
//   				cout << "fet rec : "<< fet_tp /(double)(fet_tp+fet_fn) << endl;
//    				cout << "fet fpr : "<< fet_fp /(double)(fet_tn+fet_fp) << endl;


//    	 double asb_pre5 = asb_tp5 /(double)(asb_tp5+asb_fp5); 
//    double asb_rec5 = asb_tp5 /(double)(asb_tp5+asb_fn5);
//    double asb_fpr5 = asb_fp5 /(double)(asb_tn5+asb_fp5);
 
//    double asb_f5 = 2*(asb_pre5*asb_rec5) /(double)(asb_pre5+asb_rec5);
//    ss<< asb_pre5<< "\t"<< asb_rec5<< "\t"<< asb_fpr5<< "\t"<< asb_f5 <<"\t";
   
//    double chi_pre5 = chi_tp5 /(double)(chi_tp5+chi_fp5); 
//    double chi_rec5 = chi_tp5 /(double)(chi_tp5+chi_fn5);
//    double chi_fpr5 = chi_fp5 /(double)(chi_tn5+chi_fp5);
//     double chi_f5 = 2*(chi_pre5*chi_rec5) /(double)(chi_pre5+chi_rec5);
//    ss<< chi_pre5<< "\t"<< chi_rec5<< "\t"<< chi_fpr5<< "\t"<< chi_f5 <<"\t";


//    double fet_pre5 = fet_tp5 /(double)(fet_tp5+fet_fp5); 
//    double fet_rec5 = fet_tp5 /(double)(fet_tp5+fet_fn5);
//    double fet_fpr5 = fet_fp5 /(double)(fet_tn5+fet_fp5);
// 	double fet_f5 = 2*(fet_pre5*fet_rec5) /(double)(fet_pre5+fet_rec5);
//    ss<< fet_pre5<< "\t"<< fet_rec5<< "\t"<< fet_fpr5<< "\t"<< fet_f5 <<"\t";   

   		auc_statvec.push_back(ss.str());


   		ss.str("");
  		ss<< "cp /home/yyang/data/eclip/K562/SRSF1_test/sim/asbclip.out " << "/home/yyang/projects/eCLIP/motif/chronical/08232018_simulation_rc/rc/"<< AR<< "/asbclip_"<<selected<<".txt"; 
  		cmd = ss.str();
  		system(cmd.c_str());

  		ss.str("");
  		ss<< "cp /home/yyang/data/eclip/K562/SRSF1_test/sim/chi.out " << "/home/yyang/projects/eCLIP/motif/chronical/08232018_simulation_rc/rc/"<< AR<< "/chi_"<<selected<<".txt"; 
  		cmd = ss.str();
  		system(cmd.c_str());

  		ss.str("");
  		ss<< "cp /home/yyang/data/eclip/K562/SRSF1_test/sim/fisher.out " << "/home/yyang/projects/eCLIP/motif/chronical/08232018_simulation_rc/rc/"<< AR<< "/fisher_"<<selected<<".txt"; 
  		cmd = ss.str();
  		system(cmd.c_str());

   		selected++;
	}


	}

	tio.dlwrite( "res.txt" ,  auc_statvec );
	return 0;
	
}
