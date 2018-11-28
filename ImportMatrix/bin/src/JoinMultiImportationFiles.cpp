#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>


using namespace std;


int main(int argc, char* argv[]){

    string outputfile(argv[1]);
    string outputfiletmp = outputfile+"_tmp";
    vector<string> inputfiles;
    
    //all the input files 
    for(int64_t kp = 2; kp<argc; kp++)
    {
      string paramc = argv[kp];
      inputfiles.push_back(paramc); //list of path of input files
    }
    string line;
    uint64_t Tfinal  = 0; //number of columns in final file

    for(uint64_t fi = 0; fi<inputfiles.size(); fi++){
      //loading the current input file in hash
      map<uint64_t, string> data;
      ifstream file_input_ptr(inputfiles[fi].c_str(), ios::in);
      if(!file_input_ptr){
        exit(0);
      }
      string header1, header2;
      getline(file_input_ptr, header1);
      getline(file_input_ptr, header2);
      
      vector<string> content00;
      boost::split(content00, header2, boost::is_any_of("\t"));
      uint64_t Tdata = content00.size()-1; //nombre de donnees nouvelles sans tenir compte de l'ID
      
      while (getline(file_input_ptr, line)){
        vector<string> content01;
        boost::split(content01, line, boost::is_any_of("\t"));
        string topush = "";
        for(uint64_t ic = 1; ic<=Tdata; ic++){
          topush+="\t"+content01[ic];
        }
        uint64_t indice_c = strtoull(content01[0].c_str(), NULL, 0);
        data[indice_c] = topush;
      }
      file_input_ptr.close();

      ofstream file_tmp_ptr(outputfiletmp.c_str(), ios::out);
      if(Tfinal>0){
        string lostone = ""; 
        for(uint64_t i=0;i<Tdata;i++){
          lostone+="\t0";
        }
        ifstream file_output_ptr(outputfile.c_str(), ios::in);
        //les headers
        string header1out, header2out;
        getline(file_output_ptr, header1out);
        getline(file_output_ptr, header2out);
        file_tmp_ptr << header1out << header1 << endl;
        file_tmp_ptr << header2out;
        for(uint64_t ic = 1; ic<=Tdata; ic++){
          file_tmp_ptr << "\t" << content00[ic];
        }
        file_tmp_ptr << endl;
        while(getline(file_output_ptr, line)){
          vector<string> content01;
          boost::split(content01, line, boost::is_any_of("\t"));
          uint64_t indice_c = strtoull(content01[0].c_str(), NULL, 0);
          if(data.find(indice_c)!=data.end()){
            file_tmp_ptr << line << data[indice_c] << endl;
            data.erase(indice_c);
          }
          else{
            file_tmp_ptr << line << lostone << endl;
          }
        }
      } //if(Tfinal>0)
      else{
        // les headers
        file_tmp_ptr << header1 << endl;
        file_tmp_ptr << header2 << endl;
      }
      //ecriture des nouveaux IDs
      string newone = ""; 
      for(uint64_t i=0;i<Tfinal;i++){
        newone+="\t0";
      }
      for(map<uint64_t, string>::iterator it = data.begin(); it!=data.end(); it++){
        uint64_t indice_c = it->first;
        file_tmp_ptr << indice_c << newone << data[indice_c] << endl;
      }
      file_tmp_ptr.close();
      Tfinal+=Tdata;
      int64_t resmv = rename(outputfiletmp.c_str(), outputfile.c_str()); //
    }
      





    return 0;
}

