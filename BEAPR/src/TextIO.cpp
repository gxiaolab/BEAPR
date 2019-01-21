#include "TextIO.h"

using namespace std;

TextIO::TextIO(){
;
}


std::vector< std::vector<std::string> > TextIO::dlread(std::string fname,  std::string delimiters ){
		vector<vector<string> > outvec;
		
		fstream filestr;
  	stringstream ss;
  	string str;
  	filestr.clear();
 		TOK tok;
		filestr.open(fname.c_str(), ios_base::in);
		assert(filestr.is_open());
		while(!filestr.eof()){
	 		getline(filestr,str);
			if(str.length()>0){
				vector<string> tokvec = tok.Tokenize(str,delimiters);
				outvec.push_back(tokvec);
			}
		}
	
	filestr.close();
	return outvec;
}

std::vector<std::string>  TextIO::dlread_vector(std::string fname ){
	vector<string>  outvec;
		
	fstream filestr;
  	stringstream ss;
  	string str;
  	filestr.clear();
 		TOK tok;
		filestr.open(fname.c_str(), ios_base::in);
		assert(filestr.is_open());
		while(!filestr.eof()){
	 		getline(filestr,str);
			if(str.length()>0){
				//vector<string> tokvec = tok.Tokenize(str,delimiters);
				outvec.push_back(str);
			}
		}
	
	filestr.close();
	return outvec;
}

std::vector< std::vector<double> > TextIO::dlread_double(std::string fname,  std::string delimiters ){
		vector<vector<double> > outvec;
		
		fstream filestr;
  	stringstream ss;
  	string str;
  	filestr.clear();
 		TOK tok;
		filestr.open(fname.c_str(), ios_base::in);
		assert(filestr.is_open());
		while(!filestr.eof()){
	 		getline(filestr,str);
			if(str.length()>0){
				vector<string> tokvec = tok.Tokenize(str,delimiters);
				vector<double> tmpvec;
				for(int i=0;i<tokvec.size();i++){
					tmpvec.push_back(atof(tokvec[i].c_str()));
				}
				outvec.push_back(tmpvec);
			}
		}
	
	filestr.close();
//	cout<< "ok" <<endl;
	return outvec;
}

void TextIO::dlwrite(std::string fname,  std::vector<std::string> tokvec ){
	fstream filestr;
  	stringstream ss;
  	string str;
  	filestr.clear();
 		
		filestr.open(fname.c_str(), ios_base::out);
		assert(filestr.is_open());
	
		for(int i =0;i< tokvec.size();i++){
			filestr<<tokvec[i];
			filestr<<endl;
		}
	
		filestr.close();
}

void TextIO::dlwrite(std::string fname, std::vector< std::vector<std::string> > tokvec , std::string delimiters ){
		fstream filestr;
  	stringstream ss;
  	string str;
  	filestr.clear();
 		
		filestr.open(fname.c_str(), ios_base::out);
		assert(filestr.is_open());
	
		for(int i =0;i< tokvec.size();i++){
			vector<string> tmpvec = tokvec[i];
			
			for(int j=0;j<tmpvec.size();j++){

				filestr<<tmpvec[j];
				if(j<tmpvec.size()-1){
				  filestr<<delimiters;
				}
			}
			filestr<<endl;
		}
	
		filestr.close();
}


void TextIO::dlwrite_double(std::string fname, std::vector< std::vector<double> > tokvec , std::string delimiters ){
		fstream filestr;
  	stringstream ss;
  	string str;
  	filestr.clear();
 		
		filestr.open(fname.c_str(), ios_base::out);
		assert(filestr.is_open());
	
		for(int i =0;i< tokvec.size();i++){
			vector<double> tmpvec = tokvec[i];
			
			for(int j=0;j<tmpvec.size();j++){
				filestr<<tmpvec[j]<<delimiters;
			}
			filestr<<endl;
		}
	
		filestr.close();
}

void TextIO::dlwrite_withHead(std::string fname, std::vector<std::string> HeadVec, std::vector< std::vector<std::string> > tokvec , std::string delimiters ){
		fstream filestr;
  	stringstream ss;
  	string str;
  	filestr.clear();
 		
		filestr.open(fname.c_str(), ios_base::out);
		assert(filestr.is_open());
	
		for(int i =0;i< tokvec.size();i++){
			filestr<<HeadVec[i]<<delimiters;
			vector<string> tmpvec = tokvec[i];
			
			for(int j=0;j<tmpvec.size();j++){
				filestr<<tmpvec[j]<<delimiters;
			}
			filestr<<endl;
		}
	
		filestr.close();

}

void TextIO::dlwrite_withHead(std::string fname, std::vector<std::string> HeadVec, std::vector< std::vector<double> > tokvec , std::string delimiters ){
		fstream filestr;
  	stringstream ss;
  	string str;
  	filestr.clear();
 		
		filestr.open(fname.c_str(), ios_base::out);
		assert(filestr.is_open());
	
		for(int i =0;i< tokvec.size();i++){
			filestr<<HeadVec[i]<<delimiters;
			vector<double> tmpvec = tokvec[i];
			
			for(int j=0;j<tmpvec.size();j++){
				filestr<<tmpvec[j]<<delimiters;
			}
			filestr<<endl;
		}
	
		filestr.close();

}

void TextIO::dlwrite_sparse(std::string fname, std::vector< std::vector<double> > tokvec ){
		fstream filestr;
  	stringstream ss;
  	string str;
  	filestr.clear();
 		
		filestr.open(fname.c_str(), ios_base::out);
		assert(filestr.is_open());
		int row_size = tokvec.size();
		int col_size = tokvec[0].size();
	  filestr<< row_size<<endl;
	  filestr<< col_size<<endl;
		for(int i =0;i< tokvec.size();i++){
			vector<double> tmpvec = tokvec[i];
			
			for(int j=0;j<tmpvec.size();j++){
				if(tmpvec[j]!=0){
					filestr<< i<<"_" << j<< "_"<<tmpvec[j]<<endl;
				}
			}
			
		}
	
		filestr.close();
}

std::vector< std::vector<double> > TextIO::dlread_sparse(std::string fname ){
		vector<vector<double> > outvec;
		
		fstream filestr;
  	stringstream ss;
  	string str;
  	filestr.clear();
 		TOK tok;
		filestr.open(fname.c_str(), ios_base::in);
		assert(filestr.is_open());
		getline(filestr,str);
		int row_size = atoi(str.c_str());
		getline(filestr,str);
		int col_size = atoi(str.c_str());
		
		//initialize vector
		for(int i =0; i < row_size; i++){
			vector<double> tmpvec ; 
			for(int j=0;j< col_size; j++){
				tmpvec.push_back(0);
			}
			outvec.push_back(tmpvec);
		}
		//read values
		while(!filestr.eof()){
	 		getline(filestr,str);	
			if(str.length()>0){
				vector<string> tokvec = tok.Tokenize(str,"_");
				int i = atoi(tokvec[0].c_str());
				int j = atoi(tokvec[1].c_str());
				double val = atof(tokvec[2].c_str());
				outvec[i][j] = val;
			}
		}
	
	filestr.close();
	return outvec;
}

