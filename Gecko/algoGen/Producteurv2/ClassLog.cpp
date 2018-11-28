#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <string.h>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>
#include <stdlib.h>
#include <time.h>

using namespace std;


/*!
 * \file ClassLog.cpp
 * \brief class to instantiate the write into the logs 
 * \class Class ClassLog
 * 
 * The class stores a set of features that should be update at each step of the genetic algorithm
 * The iteration is considered as a token, the _self_data is the structure to store a value at each token for each feature
 * The list of features are stored into _self_elements
 * The class provides also methods to directly store data structures into files
 * 
*/
class ClassLog{
public:

    /// \param _self_data the data to be store
    /// \param _self_elements the features to be stored for yml result
    /// \param _self_self_token_data the number of elements to be stored
    map<uint64_t, map<string, string>> _self_data;
    map<string, bool> _self_elements;
    uint64_t _self_token;

    ClassLog(){
    }

    /// \brief this will initialize the lit of features to be stored
    void init(vector<string> listOfElements, uint64_t token_init){
        for(uint64_t i = 0; i< listOfElements.size(); i++ ){
            _self_elements[listOfElements[i]] = false;
        }
        _self_token = token_init;
    }

    /// \brief for the current token, store the value to the corresponding feature
    void addFeatureToToken(string feature, string value){
        if(_self_elements.find(feature)!=_self_elements.end()){
            _self_data[_self_token][feature] = value;
            _self_elements[feature] = true;
        }
        else{
            cerr << "Warning feature " << feature << " does not belong to log structure" << endl;
        }
    }

    /// \brief if the current token is passed, go to the next
    void updateToken(uint64_t newtoken){
        for(auto it = _self_elements.begin(); it!=_self_elements.end(); it++){
            if(it->second){
                _self_elements[it->first] = false;
            }
            else{
                _self_data[_self_token][it->first] = "NA";
            }
        }
        _self_token = newtoken;
    }

    /// \brief store a list of values for the series of tokens
    void addFeatureToTokenMultiGeneration(string feature, vector<string> values, uint64_t tokenstart){
        if(_self_elements.find(feature)!=_self_elements.end()){
            for(uint64_t t = 0; t<values.size(); t++){
                //cerr << tokenstart-t << "\t" << feature << "\t" << values[t] << endl;;
                _self_data[tokenstart-t][feature] = values[t];
                if((tokenstart-t)==_self_token){
                    _self_elements[feature] = true;
                }
            }
        }
        else{
            cerr << "Warning feature " << feature << " does not belong to log structure" << endl;
        }
    }

    /// \brief write all the content of _self_data to the file (yml format)
    void writeContent(string pathToDoc){
        for(auto it = _self_elements.begin(); it!=_self_elements.end(); it++){
            if(it->second){
                _self_elements[it->first] = false;
            }
            else{
                _self_data[_self_token][it->first] = "NA";
            }
        }
        ofstream file_ao(pathToDoc.c_str(), ios::out);
        if(file_ao){
            file_ao << "{" << endl;
            file_ao << "\t\"run\": {" << endl;

            vector<uint64_t> tmp_tokens;
            vector<string> tmp_elements;
            for(auto itd = _self_data.begin(); itd!=_self_data.end(); itd++){
                tmp_tokens.push_back(itd->first);
            }
            for(auto ite = _self_elements.begin(); ite!=_self_elements.end(); ite++){
                tmp_elements.push_back(ite->first);
            }

            for(uint64_t itd = 0; itd<tmp_tokens.size(); itd++){
                uint64_t token_c = tmp_tokens[itd];
                file_ao << "\t\t\"" << token_c << "\": {" << endl;
                for(uint64_t ite = 0; ite<tmp_elements.size(); ite++){
                    file_ao << "\t\t\t\"" << tmp_elements[ite] << "\": \"" << _self_data[token_c][tmp_elements[ite]];
                    if(ite==(tmp_elements.size()-1)){
                        file_ao << "\"" << endl;
                    }
                    else{
                        file_ao << "\"," << endl;
                    }
                    
                }
                if(itd==(tmp_tokens.size()-1)){
                    file_ao << "\t\t}" << endl;
                }
                else{
                    file_ao << "\t\t}," << endl;
                }  
            }
            file_ao << "\t}" << endl;
            file_ao << "}" << endl;
            file_ao.close();
        }

    }

    /// \brief write into a file the list of occurences for each kmers
    void writeOccKmer(string pathToDoc, vector<string> sequenceID, vector<uint64_t> sequenceCounts, vector<long double> initial_sequenceScores){
        string pathFinal = pathToDoc + "_OccKmers.csv";
        string pathFinalSort = pathToDoc+"_OccKmers.csv_sort";
        ofstream file_ao(pathFinal.c_str(), ios::out);
        if(file_ao){
            for(uint64_t i =0; i <sequenceID.size(); i++){
                if(sequenceCounts[i]>0)
                    file_ao << sequenceID[i] << "," << sequenceCounts[i] <<','<<initial_sequenceScores[i]<< endl;
            }
            file_ao.close();

            string command = "sort -t, -r -nk2 " + pathFinal + " > " + pathFinalSort;
            system(command.c_str());
            command = "mv "+pathFinalSort+" "+ pathFinal;
            system(command.c_str());
        }


    }

    /// \brief write a table of integers into file with add of _tableOfIndividuals.csv extension
    void writeTableOfIndividuals(string pathToTableFile, vector<vector<uint64_t>> table){
        string pathFinal = pathToTableFile + "_tableOfIndividuals.csv";
        ofstream file_ao(pathFinal.c_str(), ios::out);
        for(uint64_t i = 0; i<table.size(); i++){
            file_ao << table[i][0];
            for(uint64_t j=1; j<table[i].size(); j++){
                file_ao << "," << table[i][j];
            }
            file_ao << endl;
        }
        file_ao.close();
    }

    /// \brief write a table of integers into file with add of _tableOfIndividuals.csv extension
    void writeOuter(string pathToTableFile, vector<int> table){            
        string pathFinal = pathToTableFile + "_outerSelected.csv";
        ofstream file_ao(pathFinal.c_str(), ios::out);
        for(uint64_t i = 0; i<table.size(); i++){
            file_ao << table[i]<<endl;
        }
        file_ao.close();
    }

    void writeGoodIndividuals(string pathToFile, vector<string> sequenceID, vector <vector<uint64_t>> tableWinnersWithOuters_individuals, vector <uint64_t> tableWinnersWithOuters_iteration,
                                        vector <long double> tableWinnersWithOuters_test, vector <long double> tableWinnersWithOuters_outer){
        cerr << "ClassLog::writeGoodIndividuals" << endl;                                            
        ofstream file_ao(pathToFile.c_str(), ios::out | ios::app);
        if(file_ao){
            for(uint64_t p = 0; p<tableWinnersWithOuters_iteration.size(); p++){

                string listKmers = sequenceID[tableWinnersWithOuters_individuals[p][0]];
                for(uint64_t i = 1; i<tableWinnersWithOuters_individuals[p].size(); i++){
                    listKmers+=","+sequenceID[tableWinnersWithOuters_individuals[p][i]];
                }
                file_ao << tableWinnersWithOuters_iteration[p] << "\t" << listKmers << "\t" << tableWinnersWithOuters_test[p] << "\t" << tableWinnersWithOuters_outer[p] << endl;
            }
            file_ao.close();
        }
        else{
            cerr << "Error while opening " << pathToFile << " in write mode" << endl;
            exit(EXIT_FAILURE);
        }
    }


    void writeGoodIndividualsMultiPop(string pathToFile, vector<string> sequenceID, vector <vector<uint64_t>> tableWinnersWithOuters_individuals, vector <uint64_t> tableWinnersWithOuters_iteration,
                                        vector <long double> tableWinnersWithOuters_test, vector <long double> tableWinnersWithOuters_outer, uint64_t npop_current){
        cerr << "ClassLog::writeGoodIndividuals" << endl;                                            
        ofstream file_ao(pathToFile.c_str(), ios::out | ios::app);
        if(file_ao){
            for(uint64_t p = 0; p<tableWinnersWithOuters_iteration.size(); p++){

                string listKmers = sequenceID[tableWinnersWithOuters_individuals[p][0]];
                for(uint64_t i = 1; i<tableWinnersWithOuters_individuals[p].size(); i++){
                    listKmers+=","+sequenceID[tableWinnersWithOuters_individuals[p][i]];
                }
                file_ao << tableWinnersWithOuters_iteration[p] << "\t" << listKmers << "\t" << tableWinnersWithOuters_test[p] << "\t" << tableWinnersWithOuters_outer[p] << "\t" << npop_current << endl;
            }
            file_ao.close();
        }
        else{
            cerr << "Error while opening " << pathToFile << " in write mode" << endl;
            exit(EXIT_FAILURE);
        }
    }

};
