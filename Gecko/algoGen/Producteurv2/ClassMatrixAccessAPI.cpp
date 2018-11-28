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
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/cxx11/is_permutation.hpp>
#include <stdlib.h>
#include <time.h>
#include <numeric>
#include <chrono>
#include <ctime>
#include "ClassMatrixAccess.cpp"


using namespace std;

/*!
* This script aims to handle binary matrix
* It relies on class ClassMatrixAccess, and gives the ability to filter binary matrix files on kmer sequence and groups
* the output are : another binary matrix by default, but can be in text 
**/


bool anotherParamAfter(int p_param, int p_argc){
    if((p_param+1)<p_argc){
        return true;
    }
    return false;
}

void help(string nameProg){
    cout << nameProg << " aims to handle binary matrix for GECKO" << endl;
    cout << "get the list of groups : " << nameProg << " -i <path to binary file> --listgroups" << endl;
    cout << "get the list of sequences : " << nameProg << " -i <path to binary file> --listkmers" << endl;
    cout << "get a submatrix : " << nameProg << " -i <path to binary file> -o <path to output file> [OPTIONS]" << endl;
    cout << "OPTIONS:" << endl;
    cout << "\t" <<  "-filterkmers <comma separated sequences> : only these sequences will be output" << endl;
    cout << "\t" <<  "-filtergroups <comma separated groups> : only these groups will be output" << endl;
    cout << "\t" <<  "-sizekmer <size of the kmers in nucleotides> : default 30" << endl;
    cout << "\t" <<  "--T : the output will be in ASCII" << endl;
    cout << "\t" <<  "--help or -h : this help" << endl;
    cout << endl << endl;
    cout << "Examples : " << endl;
    cout << nameProg << " -i /home/user/data/pathToMatrix.binary --listgroups" << endl;
    cout << nameProg << " -i /home/user/data/pathToMatrix.binary --listkmers" << endl;
    cout << nameProg << " -i /home/user/data/pathToMatrix.binary -o /home/user/data/newMatrix.binary -filterkmers AAAAAAAAAAAAAAAAAAAA,AAAAAAAAAAAAAAAAAAAC,AAAAAAAAAAAAAAAAAAAT,AAAAAAAAAAAAAAAAAAAG" << endl;
    cout << nameProg << " -i /home/user/data/pathToMatrix.binary -o /home/user/data/newMatrix.binary -filtergroups control,diagnostic1,diagnostic2" << endl;
    cout << nameProg << " -i /home/user/data/pathToMatrix.binary -o /home/user/data/newMatrix.binary -filterkmers AAAAAAAAAAAAAAAAAAAA,AAAAAAAAAAAAAAAAAAAC -filtergroups control,diagnostic1,diagnostic2" << endl;

}

int main(int argc, char *argv[]){

    
    /*************************
     * Parameters
     * To be completed
     * **********************/
    map<string, string> parameters;

    parameters["sizekmer"] = "30"; //default
    for(int param = 1; param < argc; param++) {
        string param_c(argv[param]);
        if(param_c.compare("-i")==0){
            if(anotherParamAfter(param,argc)){
                string param2_c(argv[param+1]);
                parameters["i"] = param2_c;
                param++;
            }
            else{
                cerr << "parameter -i must be followed by another parameter" << endl;
                return EXIT_FAILURE;
            }
        }

        if(param_c.compare("-o")==0){
            if(anotherParamAfter(param,argc)){
                string param2_c(argv[param+1]);
                parameters["o"] = param2_c;
                param++;
            }
            else{
                cerr << "parameter -o must be followed by another parameter" << endl;
                return EXIT_FAILURE;
            }
        }

        if(param_c.compare("-filterkmers")==0){
            if(anotherParamAfter(param,argc)){
                string param2_c(argv[param+1]);
                parameters["filterkmers"] = param2_c;
                param++;
            }
            else{
                cerr << "parameter -filterkmers must be followed by another parameter" << endl;
                return EXIT_FAILURE;
            }
        }

        if(param_c.compare("-filtergroups")==0){
            if(anotherParamAfter(param,argc)){
                string param2_c(argv[param+1]);
                parameters["filtergroups"] = "param2_c";
                param++;
            }
            else{
                cerr << "parameter -filtergroups must be followed by another parameter" << endl;
                return EXIT_FAILURE;
            }
        }

        if(param_c.compare("-sizekmer")==0){
            if(anotherParamAfter(param,argc)){
                string param2_c(argv[param+1]);
                parameters["sizekmer"] = param2_c;
                param++;
            }
            else{
                cerr << "parameter -sizekmer must be followed by another parameter" << endl;
                return EXIT_FAILURE;
            }
        }

        if(param_c.compare("--T")==0){
            parameters["T"] = "1";
        }

        if(param_c.compare("--listgroups")==0){
            parameters["listgroups"] = "1";
        }

        if(param_c.compare("--listkmers")==0){
            parameters["listkmers"] = "1";
        }

        if(param_c.compare("--fastakmers")==0){
            parameters["fastakmers"] = "1";
        }


        if( (param_c.compare("--help")==0) or (param_c.compare("-h")==0)){
            help(argv[0]);
            return EXIT_SUCCESS;
        }

    }
    
    /*************************
     * Test of parameters
     * To be completed
     * **********************/
    ClassMatrixAccess CMA_input;
    bool state_test = false; //indicator to stop or continue the execution. Will be false if a processus is initiated 

    if(parameters.find("i")!=parameters.end()){
        ClassMatrixAccess tmp = ClassMatrixAccess(parameters["i"], 1);
        CMA_input = tmp;
        state_test = CMA_input.isItABinaryFile();
    }



    if(!state_test){
        help(argv[0]);
        return EXIT_FAILURE;
    }
    

    /*************************
     * Execution
     * **********************/
    //get the list of groups to standard output
    if(state_test and parameters.find("fastakmers")!=parameters.end()){
        vector<string> listKmers = CMA_input.getSequenceIDs();
        for(uint64_t k = 0; k<listKmers.size(); k++){
            cout << "> kmer" << k+1 << endl;
            cout << CMA_input.idseq2fasta(stoull(listKmers[k]),30) << endl;
        }
        state_test = false;
    }

    //get the list of groups to standard output
    if(state_test and parameters.find("listgroups")!=parameters.end()){
        vector<string> listGroups = CMA_input.getGroups();
        for(uint64_t g = 0; g<listGroups.size(); g++){
            cout << listGroups[g] << endl;
        }
        state_test = false;
    }

    //get the list of kmers to standard output
    if(state_test and parameters.find("listkmers")!=parameters.end()){
        vector<string> listKmers = CMA_input.getSequenceIDs();
        for(uint64_t k = 0; k<listKmers.size(); k++){
            cout << CMA_input.idseq2fasta(stoull(listKmers[k]),30) << endl;
        }
        state_test = false;
    }

    //output in a file
    if(state_test and parameters.find("o")!=parameters.end()){
        state_test = false;
        bool stateGroup = false; //true if a filter on groups is required
        bool stateKmers = false; //true if a batch of sequence is given. indices are relative to getSequenceIDs result
        vector<uint64_t> indice_groups; // the indice in the matrix of required groups. indices are absolute to matrix ranks
        vector<uint64_t> indice_kmers; // the indice in the matrix of required kmers

        vector<string> listGroups  = CMA_input.getGroups();;
        vector<string> listKmers   = CMA_input.getSequenceIDs();;
        vector<string> headerNames = CMA_input.getNames();
        
        ///////////////////////////////////////////////////////////////////////
        //get indice according filters
        ///////////////////////////////////////////////////////////////////////
        if(parameters.find("filtergroups") !=parameters.end()){
            stateGroup = true;
            indice_groups.push_back(0); //the first column is required
            map<string, bool> requiredGroups;
            vector<string> contentgroups;
            boost::split(contentgroups, parameters["filtergroups"], boost::is_any_of(","));

            for(uint64_t i = 0; i<contentgroups.size(); i++){
                requiredGroups[contentgroups[i]] = true;
            }

            for(uint64_t i = 0; i<listGroups.size(); i++){
                if(requiredGroups.find(listGroups[i])!= requiredGroups.end()){
                    indice_groups.push_back(i+1);
                }
            }
        }

        if(parameters.find("filterkmers") !=parameters.end()){
            stateKmers = true;
            map<string, bool> requiredKmers;
            vector<string> contentkmers;
            boost::split(contentkmers, parameters["filterkmers"], boost::is_any_of(","));

            for(uint64_t i = 0; i<contentkmers.size(); i++){
                uint64_t indice_c = CMA_input.getIndice(contentkmers[i],0);
                requiredKmers[to_string(indice_c)] = true;
            }

            for(uint64_t i = 0; i<listKmers.size(); i++){
                if(requiredKmers.find(listKmers[i])!= requiredKmers.end()){
                    indice_kmers.push_back(i);
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////
        //opening the output file
        ///////////////////////////////////////////////////////////////////////
        ofstream ofs;
        if(parameters.find("T")!=parameters.end()){
            ofs.open(parameters["o"].c_str(), ios_base::out);
        }
        else{
            ofs.open(parameters["o"].c_str(), ios_base::binary);
            float ID = 421.21;
            uint64_t NumberOfLinesIncludingLabel = CMA_input._self_nbLines;
            uint64_t NumberOfColumnsIncludingLabel = CMA_input._self_nbColumns;
            if(stateGroup)
                NumberOfColumnsIncludingLabel = indice_groups.size();
            if(stateKmers)
                NumberOfLinesIncludingLabel = 2+indice_kmers.size();
            ofs.write(reinterpret_cast<const char*>(&ID), sizeof(float));
            ofs.write(reinterpret_cast<const char*>(&NumberOfLinesIncludingLabel), sizeof(uint64_t));
            ofs.write(reinterpret_cast<const char*>(&NumberOfColumnsIncludingLabel), sizeof(uint64_t));
        }

        ///////////////////////////////////////////////////////////////////////
        //Headers 
        ///////////////////////////////////////////////////////////////////////
        vector<string> realNames; //different of headerNames if stateGroup
        vector<string> realGroups; //different of listGroups if stateGroup
        if(stateGroup){
            for(uint64_t i = 0; i<indice_groups.size(); i++){
                realNames.push_back(headerNames[indice_groups[i]]);
                if(i==0){
                    realGroups.push_back("groups");
                }
                else{
                    realGroups.push_back(listGroups[indice_groups[i]-1]);
                }
            }
        }
        else{
            realNames = headerNames;
            realGroups = listGroups;
            realGroups.insert(realGroups.begin(),"groups");
        }


        if(parameters.find("T")!=parameters.end()){
            ofs << realNames[0];
            for(uint64_t i = 1 ; i<realNames.size() ; i++){
                ofs << "\t" << realNames[i];
            }
            ofs << endl;
            ofs << realGroups[0];
            for(uint64_t i = 1 ; i<realGroups.size() ; i++){
                ofs << "\t" << realGroups[i];
            }
            ofs << endl;
        }
        else{
            for(uint64_t i = 0 ; i<realNames.size() ; i++){
                char bufhead[64];
                strcpy(bufhead, realNames[i].c_str());
                ofs.write(reinterpret_cast<const char*>(bufhead), 64);
            }
            for(uint64_t i = 0 ; i<realGroups.size() ; i++){
                char bufhead[64];
                strcpy(bufhead, realGroups[i].c_str());
                ofs.write(reinterpret_cast<const char*>(bufhead), 64);
            }
        }

        ///////////////////////////////////////////////////////////////////////
        //Rest of the file
        ///////////////////////////////////////////////////////////////////////
        //open the binary input file
        ifstream file_in_ptr(parameters["i"].c_str(), ios_base::binary);
        if(!file_in_ptr){
            exit(EXIT_FAILURE);
        }
        uint64_t sizeOfHeader = sizeof(float)+2*sizeof(uint64_t)+(CMA_input._self_nbColumns)*64*2;
        uint64_t sizeOfLine = sizeof(uint64_t)+(CMA_input._self_nbColumns-1)*sizeof(float);

        if(stateKmers){
            for(uint64_t i = 0; i < indice_kmers.size(); i++){
                //move to the kmer
                uint64_t offset = sizeOfHeader+sizeOfLine*indice_kmers[i];
                file_in_ptr.seekg( offset,  ios_base::beg );

                //read the vector
                uint64_t seqID_tmp;
                vector<float> vect_tmp;
                file_in_ptr.read(reinterpret_cast<char*>(&seqID_tmp), sizeof(uint64_t));
                for(uint64_t j = 1; j < CMA_input._self_nbColumns; j++){
                    float value_tmp;
                    file_in_ptr.read(reinterpret_cast<char*>(&value_tmp), sizeof(float) );
                    vect_tmp.push_back(value_tmp);
                }

                //output
                if(parameters.find("T")==parameters.end()){
                    //write the ID and desired sample values
                    ofs.write(reinterpret_cast<const char*>(&seqID_tmp), sizeof(uint64_t));
                    if(stateGroup){
                        //get some values
                        for(uint64_t j = 1; j < indice_groups.size() ; j++){
                            float ldtmp = vect_tmp[indice_groups[j]-1];
                            ofs.write(reinterpret_cast<const char*>(&ldtmp), sizeof(float));
                        }
                    }
                    else{
                        //get all values
                        for(uint64_t j = 0; j < vect_tmp.size() ; j++){
                            float ldtmp = vect_tmp[j];
                            ofs.write(reinterpret_cast<const char*>(&ldtmp), sizeof(float));
                        }
                    }
                }
                else{
                    ofs << seqID_tmp;
                    if(stateGroup){
                        for(uint64_t j = 1; j < indice_groups.size() ; j++){
                            ofs << "\t" << vect_tmp[indice_groups[j]-1];
                        }
                    }
                    else{
                        for(uint64_t j = 0; j < vect_tmp.size() ; j++){
                            ofs << "\t" << vect_tmp[j];
                        }
                    }
                    ofs << endl;
                }

            }
        }
        else{
            //move to the 1st kmer
            uint64_t offset = sizeOfHeader;
            file_in_ptr.seekg( offset,  ios_base::beg );

            for(uint64_t i = 2; i < CMA_input._self_nbLines; i++){
                uint64_t seqID_tmp;
                vector<float> vect_tmp;

                file_in_ptr.read(reinterpret_cast<char*>(&seqID_tmp), sizeof(uint64_t));
                for(uint64_t j = 1; j < CMA_input._self_nbColumns; j++){
                    float value_tmp;
                    file_in_ptr.read(reinterpret_cast<char*>(&value_tmp), sizeof(float) );
                    vect_tmp.push_back(value_tmp);
                }

                //output
                if(parameters.find("T")==parameters.end()){
                    ofs.write(reinterpret_cast<const char*>(&seqID_tmp), sizeof(uint64_t));
                    if(stateGroup){
                        //get some values
                        for(uint64_t j = 1; j < indice_groups.size() ; j++){
                            float ldtmp = vect_tmp[indice_groups[j]-1];
                            ofs.write(reinterpret_cast<const char*>(&ldtmp), sizeof(float));
                        }
                    }
                    else{
                        //get all values
                        for(uint64_t j = 0; j < vect_tmp.size() ; j++){
                            float ldtmp = vect_tmp[j];
                            ofs.write(reinterpret_cast<const char*>(&ldtmp), sizeof(float));
                        }
                    }
                }
                else{
                    ofs << seqID_tmp;
                    if(stateGroup){
                        for(uint64_t j = 1; j < indice_groups.size() ; j++){
                            ofs << "\t" << vect_tmp[indice_groups[j]-1];
                        }
                    }
                    else{
                        for(uint64_t j = 0; j < vect_tmp.size() ; j++){
                            ofs << "\t" << vect_tmp[j];
                        }
                    }
                    ofs << endl;
                }
            }
        }
        
       
        ofs.close();
    }

    if(state_test){
        help(argv[0]);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;    
}
