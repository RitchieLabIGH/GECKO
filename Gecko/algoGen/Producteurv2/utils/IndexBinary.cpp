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

//! 
//! \brief  Cut a binary gecko matrix into a serie of matrix
//! the main with contain the general informations + 2 headers + seqID
//! then the values of corresponding kmers will be divided into several pieces
//!

int main(int argc, char* argv[]){

    string path2MatrixBinaryInput(argv[1]); 
    string path2MatrixOutput(argv[2]);
    uint64_t NKmerPerFile = stoull(argv[3]);

    /****************************************************
     * Verify that the file is a good matrix input
     * **************************************************/
    bool isBinaryInput = false;
    ifstream file_in_ptr(path2MatrixBinaryInput.c_str(), ios_base::binary);
    if(!file_in_ptr){
        cerr << "ClassMatrixAccess::isItABinaryFile Error the file " << path2MatrixBinaryInput << " does not exist" << endl;
        exit(EXIT_FAILURE);
    }
    float value;
    uint64_t numberOfLines, numberOfColumns;
    file_in_ptr.read(reinterpret_cast<char*>(&value), sizeof(float));
        
    if(round(value-421.21)==0){
        file_in_ptr.read(reinterpret_cast<char*>(&numberOfLines), sizeof(uint64_t));
        file_in_ptr.read(reinterpret_cast<char*>(&numberOfColumns), sizeof(uint64_t));
        isBinaryInput = true;
    }
    if(!isBinaryInput){
        return EXIT_FAILURE;
    }

    /****************************************************
     * Open the main output and begins to write the 
     * content
     * **************************************************/
    ofstream file_ao(path2MatrixOutput.c_str(), ios_base::binary);
    if (!file_ao) {
        cerr << "Error while opening " << path2MatrixOutput << " in write binary mode" << endl;
        exit(EXIT_FAILURE);
    }
    float ID = 233.33;
    file_ao.write(reinterpret_cast<const char*>(&ID), sizeof(float));
    file_ao.write(reinterpret_cast<const char*>(&numberOfLines), sizeof(uint64_t));
    file_ao.write(reinterpret_cast<const char*>(&numberOfColumns), sizeof(uint64_t));
    file_ao.write(reinterpret_cast<const char*>(&NKmerPerFile), sizeof(uint64_t));
    //the names
    vector<string> Vnames;
    for(uint64_t i = 0; i<numberOfColumns ; i++){
        char bufread[64];
        file_in_ptr.read(reinterpret_cast<char*>(bufread), 64);
        string group_c(bufread);
        Vnames.push_back(group_c);
    }
    //the the groups
    vector<string> Vgroups;
    for(uint64_t i = 0; i<numberOfColumns ; i++){
        char bufread[64];
        file_in_ptr.read(reinterpret_cast<char*>(bufread), 64);
        string group_c(bufread);
        Vgroups.push_back(group_c);
    }
    for (uint64_t i = 0; i < Vnames.size(); i++) {
        char bufhead[64];
        strcpy(bufhead, Vnames[i].c_str());
        file_ao.write(reinterpret_cast<const char*>(bufhead), 64);
    }
    for (uint64_t i = 0; i < Vgroups.size(); i++) {
        char bufhead[64];
        strcpy(bufhead, Vgroups[i].c_str());
        file_ao.write(reinterpret_cast<const char*>(bufhead), 64);
    }


    /****************************************************
     * Now for each kmer, write the ID into the main file
     * and values into pieces
     * **************************************************/
    ofstream file_data_ao;
    uint64_t sliceNumber = 0;
    uint64_t kmer_count=0;
    string current_path_output_datafile = path2MatrixOutput+"_"+to_string(sliceNumber);

    file_data_ao.open(current_path_output_datafile.c_str(), ios_base::binary);
    
    for(uint64_t kgeneral = 2; kgeneral<numberOfLines; kgeneral++){
        kmer_count+=1;
        uint64_t seqID;
        vector<float> valuesOfKmer;
        file_in_ptr.read(reinterpret_cast<char*>(&seqID), sizeof(uint64_t));
        for(uint64_t i = 1; i<numberOfColumns; i++){
            float value_tmp;
            file_in_ptr.read(reinterpret_cast<char*>(&value_tmp), sizeof(float) );
            valuesOfKmer.push_back(value_tmp);
        }
        file_ao.write(reinterpret_cast<const char*>(&seqID), sizeof(uint64_t));
        for (uint64_t i = 0; i < valuesOfKmer.size(); i++) {
            float ldtmp = valuesOfKmer[i];
            file_data_ao.write(reinterpret_cast<const char*>(&ldtmp), sizeof(float));
        }

        if(kmer_count==NKmerPerFile){
            file_data_ao.close();
            sliceNumber+=1;
            current_path_output_datafile = path2MatrixOutput+"_"+to_string(sliceNumber);
            file_data_ao.open(current_path_output_datafile.c_str(), ios_base::binary);
            kmer_count=0;
        }

    }
    file_data_ao.close();
    file_in_ptr.close();
    file_ao.close();






    return EXIT_SUCCESS;
}