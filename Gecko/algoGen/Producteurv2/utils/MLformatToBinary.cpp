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
//! \brief  Transforms the matrix (originally in ASCII )
//!

int main(int argc, char* argv[]){
    string path2MatrixText(argv[1]); //line x col : kmer x groups
    //string path2MatrixText = "/home/aubin/ritchie/gecko/peripheralBlood/PRJNA391912_norm_discret_Clean_realCounts.matrix";
    string path2MatrixBinary(argv[2]);
    //string path2MatrixBinary = "test.binary";

    /********************************************
    //Get the number of col plus number of lines
    ********************************************/
    ifstream file_ai(path2MatrixText.c_str(), ios::in);
    if (!file_ai) {
        cerr << "Error while opening " << path2MatrixText << " in read mode" << endl;
        exit(EXIT_FAILURE);
    }
    string line;
    vector<string> content;

    uint64_t NumberOfColumnsIncludingLabel = 0;
    uint64_t NumberOfLinesIncludingLabel   = 0;
    getline(file_ai, line);
    NumberOfLinesIncludingLabel++;
    getline(file_ai, line);
    NumberOfLinesIncludingLabel++;

    boost::split(content, line, boost::is_any_of("\t"));
    NumberOfColumnsIncludingLabel = content.size();
    while(getline(file_ai, line)){
        NumberOfLinesIncludingLabel++;
    }
    file_ai.close();
    cerr << "The matrix has " << NumberOfLinesIncludingLabel << " lines and " << NumberOfColumnsIncludingLabel << " columns" << endl;


    /*******************************************
    //read the file and binary transformation
    *******************************************/
    ifstream file_ai2(path2MatrixText.c_str(), ios::in);
    if (!file_ai2) {
        cerr << "Error while opening " << path2MatrixText << " in read mode" << endl;
        exit(EXIT_FAILURE);
    }
    ofstream file_ao(path2MatrixBinary.c_str(), ios_base::binary);
    if (!file_ao) {
        cerr << "Error while opening " << path2MatrixBinary << " in write binary mode" << endl;
        exit(EXIT_FAILURE);
    }


    /*********************************************
    //writing the number of lines, columns and ID 
    *********************************************/
    float ID = 421.21;
    file_ao.write(reinterpret_cast<const char*>(&ID), sizeof(float));
    file_ao.write(reinterpret_cast<const char*>(&NumberOfLinesIncludingLabel), sizeof(uint64_t));
    file_ao.write(reinterpret_cast<const char*>(&NumberOfColumnsIncludingLabel), sizeof(uint64_t));
    
    
    /*********************************************
    //writing the first lines 
    *********************************************/
    getline(file_ai2, line);
    boost::split(content, line, boost::is_any_of("\t"));
    content[0]="name";

    for (uint64_t i = 0; i < content.size(); i++) {
        char bufhead[64];
        strcpy(bufhead, content[i].c_str());
        file_ao.write(reinterpret_cast<const char*>(bufhead), 64);
    }
    getline(file_ai2, line);
    boost::split(content, line, boost::is_any_of("\t"));

    for (uint64_t i = 0; i < content.size(); i++) {
        char bufhead[64];
        strcpy(bufhead, content[i].c_str());
        file_ao.write(reinterpret_cast<const char*>(bufhead), 64);
    }


   /*********************************************
    //writing the other lines
    *********************************************/
    while(getline(file_ai2, line)){
        boost::split(content, line, boost::is_any_of("\t"));
        char buf[64];

        uint64_t inttmp = stoull(content[0].c_str(), NULL);
        file_ao.write(reinterpret_cast<const char*>(&inttmp), sizeof(uint64_t));
        
        for (uint64_t i = 1; i < content.size(); i++) {
            float ldtmp = stof(content[i].c_str(), NULL);
            file_ao.write(reinterpret_cast<const char*>(&ldtmp), sizeof(float));
        }
    }
    file_ai2.close();
    file_ao.close();

   
    return EXIT_SUCCESS;


}