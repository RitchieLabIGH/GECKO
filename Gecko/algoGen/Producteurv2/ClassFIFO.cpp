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
//#include <boost/filesystem.hpp>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <thread>

using namespace std;

class ClassFIFO{
public:
    ClassFIFO(){
    }

    void writeML(string uid, map<uint64_t, vector< vector<long double>>> train, map<uint64_t, vector< vector<long double>>> test, uint64_t nbIndividuals, uint64_t nbsamplesPerIndividual, uint64_t nbkmers){

        //tailles
        //long double sldib = sizeof(long double); //size long double in bytes
        //cout<<"writeML test"<<endl;
        //ecriture au format text pour l'instant
        string filetmp = uid+"_tmp";
        string filefinal = uid+"_request";
        ofstream file_ao(filetmp.c_str(), ios::out);

        //mise en memoire
        for(auto it = train.begin(); it!=train.end(); it++){
            uint64_t guy = it->first;
            //le train
            for(uint64_t i = 0; i < it->second.size(); i++){
                file_ao << it->second[i][0];
                for(uint64_t j = 1; j<it->second[i].size(); j++){
                    file_ao << "," << it->second[i][j];
                }
                file_ao << endl;
            }
            //le test
            for(uint64_t i = 0; i < test[guy].size(); i++){
                file_ao << test[guy][i][0];
                for(uint64_t j = 1; j<test[guy][i].size(); j++){
                    file_ao << "," << test[guy][i][j];
                }
                file_ao << endl;
            }
        }
        file_ao.close();

        //move
        rename( filetmp.c_str() , filefinal.c_str() );
    }

    void writeMLbinary(string uid, 
                map<uint64_t, vector< vector<long double>>> train, 
                map<uint64_t, vector< vector<long double>>> test, 
                uint64_t nbIndividuals, 
                uint64_t nbsamplesPerIndividual, 
                uint64_t nbkmers){

        //tailles
        uint32_t sldib     = sizeof(long double); //size long double in bytes
        uint32_t nblines   = 0;
        uint32_t nbcolumns = 0;
        
        //cout<<"writeMLbinaries"<<endl;
        
        for(auto it = train.begin(); it!=train.end(); it++){
            uint64_t guy = it->first;
            nblines+=it->second.size();
            nblines+=test[guy].size();
            nbcolumns = test[guy][0].size();
        }
        
        //ecriture au format binary
        string filetmp = uid+"_tmp";
        string filefinal = uid+"_request";
        
        ofstream file_ao(filetmp.c_str(), ios::out | ios::binary);
        long double tmpval;
        //file_ao << sldib;
        file_ao.write ((char*)&sldib,sizeof(sldib));
        file_ao.write ((char*)&nblines,sizeof(sldib));
        file_ao.write ((char*)&nbcolumns,sizeof(sldib));

        //cout<< sldib<<" line "<< nblines <<" col "<< nbcolumns<<endl;
        //mise en memoire
        for(auto it = train.begin(); it!=train.end(); it++){
            uint64_t guy = it->first;
            //le train
            for(uint64_t i = 0; i < it->second.size(); i++){
                for(uint64_t j = 0; j<it->second[i].size(); j++){
                    tmpval = it->second[i][j];
                    file_ao.write ((char*)&tmpval,sizeof (tmpval));
                }
            }
            //le test
            for(uint64_t i = 0; i < test[guy].size(); i++){
                for(uint64_t j = 0; j<test[guy][i].size(); j++){
                    tmpval = test[guy][i][j];
                    file_ao.write ((char*)&tmpval,sizeof(tmpval));
                }
            }
        }
        file_ao.close();

        //move
        rename( filetmp.c_str() , filefinal.c_str() );
    }
    void writeLastMLbinary(string uid, map<uint64_t, vector< vector<long double>>> train, map<uint64_t, vector< vector<long double>>> test, 
        uint64_t nbsamplesPerIndividual, uint64_t nbkmers){
        //tailles
        uint32_t sldib     = sizeof(long double); //size long double in bytes
        uint32_t nblines   = 0;
        uint32_t nbcolumns = 0;
        //uint64_t kguy = 0;
    //    cout<<"nbIndividualsBegin "<<nbIndividualsBegin<<"   nbIndividualsEnd "<<nbIndividualsEnd<<endl;
        for(auto it = train.begin(); it!=train.end(); it++){
            //if((kguy>=nbIndividualsBegin) and (kguy<nbIndividualsEnd)){
                uint64_t guy = it->first;
                nblines+=it->second.size();
                nblines+=test[guy].size();
                nbcolumns = test[guy][0].size();
            //}
            //kguy++;
        }

        //ecriture au format text pour l'instant
        string filetmp = uid+"_tmp";
        string filefinal = uid+"_request";
        ofstream file_ao(filetmp.c_str(), ios::out| ios::binary);
        //cout<<"sizeof(sldib)"<<sizeof(sldib)<<endl;
        file_ao.write ((char*)&sldib,sizeof(sldib));
        file_ao.write ((char*)&nblines,sizeof(sldib));
        file_ao.write ((char*)&nbcolumns,sizeof(sldib));
        long double tmpval;
        //cout<< sldib<<" line "<< nblines <<" col "<< nbcolumns<<endl;
        //mise en memoire
        //kguy = 0;
        for(auto it = train.begin(); it!=train.end(); it++){
            //if((kguy>=nbIndividualsBegin) and (kguy<nbIndividualsEnd)){
                uint64_t guy = it->first;
                //le train
                for(uint64_t i = 0; i < it->second.size(); i++){
                    for(uint64_t j = 0; j<it->second[i].size(); j++){
                        tmpval = it->second[i][j];
                        file_ao.write ((char*)&tmpval,sizeof (tmpval));
                        //file_ao << "," << it->second[i][j];
                    }
                    //file_ao << endl;
                }
                //le test
                for(uint64_t i = 0; i < test[guy].size(); i++){
                    //file_ao << test[guy][i][0];
                    for(uint64_t j = 0; j<test[guy][i].size(); j++){
                        tmpval = test[guy][i][j];
                        file_ao.write ((char*)&tmpval,sizeof (tmpval));
                        //file_ao << "," << test[guy][i][j];
                    }
                    //file_ao << endl;
                }
           // }
           // kguy++;
        }
        file_ao.close();

        //move
        rename( filetmp.c_str() , filefinal.c_str() );
    }

    vector<long double> readML(string uid){
        string filefinal = uid+"_response";
        bool state = true;
        vector<long double> result;

        while(state){
            ifstream fai(filefinal.c_str(), ios::in);
            state = fai.fail();
            if(state){
                this_thread::sleep_for(std::chrono::milliseconds(100));
            }
            else{
                long double value;
                while(fai >> value){
                    result.push_back(value);
                }
                fai.close();
            }
        }
        remove(filefinal.c_str());

        return result;
    }

    void writeLastML(string uid, map<uint64_t, vector< vector<long double>>> train, map<uint64_t, vector< vector<long double>>> test, 
        uint64_t nbIndividualsBegin, uint64_t nbIndividualsEnd, uint64_t nbsamplesPerIndividual, uint64_t nbkmers){
               
        //ecriture au format text pour l'instant
        string filetmp = uid+"_tmp";
        string filefinal = uid+"_request";
        ofstream file_ao(filetmp.c_str(), ios::out);
        
        //mise en memoire
        uint64_t kguy = 0;
        for(auto it = train.begin(); it!=train.end(); it++){
            if((kguy>=nbIndividualsBegin) and (kguy<nbIndividualsEnd)){
                uint64_t guy = it->first;
                //le train
                for(uint64_t i = 0; i < it->second.size(); i++){
                    file_ao << it->second[i][0];
                    for(uint64_t j = 1; j<it->second[i].size(); j++){
                        file_ao << "," << it->second[i][j];
                    }
                    file_ao << endl;
                }
                //le test
                for(uint64_t i = 0; i < test[guy].size(); i++){
                    file_ao << test[guy][i][0];
                    for(uint64_t j = 1; j<test[guy][i].size(); j++){
                        file_ao << "," << test[guy][i][j];
                    }
                    file_ao << endl;
                }
                cout << endl;
            }
            kguy++;
        }
        file_ao.close();

        //move
        rename( filetmp.c_str() , filefinal.c_str() );
    }

    void killML(string uid){
        //ecriture au format text pour l'instant
        string filefinal = uid+"_stop";
        ofstream file_ao(filefinal.c_str(), ios::out);
        file_ao.close();

    }

    

};
