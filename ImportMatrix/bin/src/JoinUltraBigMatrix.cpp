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
    string file_output_final = argv[1];
    string file_input1 = argv[2];
    string file_input2 = argv[3];

    uint64_t size_sample  = 500000000; //200M

    /**
     * ouverture des fichiers de sortie
     * **/
    ofstream file_output_final_ptr(file_output_final.c_str(), ios::out);

    /**
     * recuperation des infos sur les fichiers input
     * **/
    uint64_t nbline1 = 0;
    uint64_t nbline2 = 0;
    uint64_t nbcol1 = 0;
    uint64_t nbcol2 = 0;

    ifstream file_input_ptr1(file_input1.c_str(), ios::in);
    ifstream file_input_ptr2(file_input2.c_str(), ios::in);

    if(!file_input_ptr1 or !file_input_ptr2){
        return EXIT_FAILURE;
    }

    string file1_header1, file1_header2, file2_header1, file2_header2;
    string line;
    getline(file_input_ptr1, file1_header1); getline(file_input_ptr1, file1_header2);
    getline(file_input_ptr2, file2_header1); getline(file_input_ptr2, file2_header2);

    vector<string> content00;
    boost::split(content00, file1_header2, boost::is_any_of("\t"));
    nbcol1 = content00.size()-1; //nombre de donnees nouvelles sans tenir compte de l'ID
    boost::split(content00, file2_header2, boost::is_any_of("\t"));
    nbcol2 = content00.size()-1; //nombre de donnees nouvelles sans tenir compte de l'ID
    cerr << nbcol1 << " samples for file1 and " << nbcol2 << " samples for file2" << endl;

    while(getline(file_input_ptr1, line)){
        nbline1++;
    }
    while(getline(file_input_ptr2, line)){
        nbline2++;
    }
    cerr << nbline1 << " kmers for file1 and " << nbline2 << " kmers for file2" << endl;
    file_input_ptr1.close();
    file_input_ptr2.close();

    string h1 = file1_header1;
    string h2 = file1_header2;
    boost::split(content00, file2_header1, boost::is_any_of("\t"));
    for(uint64_t i = 1; i < content00.size(); i++){
      h1+="\t"+content00[i];
    }
    boost::split(content00, file2_header2, boost::is_any_of("\t"));
    for(uint64_t i = 1; i < content00.size(); i++){
      h2+="\t"+content00[i];
    }

    file_output_final_ptr << h1 << endl;
    file_output_final_ptr << h2 << endl;



    /**
     * traitement par paquet de file1
     * **/
    uint64_t nbpaquet1 = nbline1/size_sample;
    uint64_t nbpaquet2 = nbline2/size_sample;
    cerr << "nbpaquet1 = " << nbpaquet1 << endl;
    cerr << "nbpaquet2 = " << nbpaquet2 << endl;

    uint64_t log_cpt_total1 = 0;
    uint64_t log_cpt_total2 = 0;

    for(uint64_t paquet = 0; paquet<nbpaquet1; paquet++){
        cerr << "paquet " << paquet << endl;
        uint64_t log_cpt_tmp = 0;
        //pointeur à la bonne ligne
        ifstream file_ptr1(file_input1.c_str(), ios::in);
        for(uint64_t cpt = 0; cpt<(2+paquet*size_sample); cpt++){
            getline(file_ptr1, line);
        }
        cerr << 2+paquet*size_sample << " lines ignored" << endl;

        //creation dictionnaire
        map<uint64_t, string> dictionnaire;
        for(uint64_t cpt = 0; cpt<size_sample; cpt++){
            getline(file_ptr1, line);
            vector<string> content;
            boost::split(content, line, boost::is_any_of("\t"));
            dictionnaire[strtoull(content[0].c_str(), NULL, 0)] = line;
        }
        file_ptr1.close();
        cerr << size_sample << " kmers read" << endl;

        //lecture file2
        ifstream file_ptr2(file_input2.c_str(), ios::in);
        while(getline(file_ptr2, line)){
            vector<string> content;
            boost::split(content, line, boost::is_any_of("\t"));
            uint64_t idc = strtoull(content[0].c_str(), NULL, 0);
            map<uint64_t, string>::iterator it = dictionnaire.find(idc);
            if(it!=dictionnaire.end()){
                log_cpt_tmp++;
                //ecriture dans output
                string ltw = it->second;
                for(uint64_t i = 1; i < content.size(); i++){
                  ltw+="\t"+content[i];
                }
                file_output_final_ptr << ltw << endl;

                //suppression de l'entree  dictionnaire
                dictionnaire.erase (it);
            }
        }
        file_ptr2.close();
        cerr << log_cpt_tmp << " found in 2nd file" << endl;

        //traitement des fails
        log_cpt_tmp = 0;
        string basis = "";
        for(uint64_t i = 0; i<nbcol2; i++){
            basis=basis+"\t0";
        }
        for(map<uint64_t, string>::iterator it = dictionnaire.begin(); it!=dictionnaire.end(); it++){
            file_output_final_ptr << it->second << basis << endl;
            log_cpt_tmp++;
        }
        cerr << log_cpt_tmp << " kmers written" << endl << endl;

    }
    if(true){
        ifstream file_ptr1(file_input1.c_str(), ios::in);
        //pointeur à la bonne ligne
        for(uint64_t cpt = 0; cpt<(2+nbpaquet1*size_sample); cpt++){
            getline(file_ptr1, line);
        }
        //creation dictionnaire
        map<uint64_t, string> dictionnairefinal;
        while(getline(file_ptr1, line)){
            vector<string> content;
            boost::split(content, line, boost::is_any_of("\t"));
            dictionnairefinal[strtoull(content[0].c_str(), NULL, 0)] = line;
        }
        file_ptr1.close();
        //lecture file2
        ifstream file_ptr2(file_input2.c_str(), ios::in);
        while(getline(file_ptr2, line)){
            vector<string> content;
            boost::split(content, line, boost::is_any_of("\t"));
            uint64_t idc = strtoull(content[0].c_str(), NULL, 0);
            map<uint64_t, string>::iterator it = dictionnairefinal.find(idc);
            if(it!=dictionnairefinal.end()){
                //ecriture dans output
                string ltw = it->second;
                for(uint64_t i = 1; i < content.size(); i++){
                  ltw+="\t"+content[i];
                }
                file_output_final_ptr << ltw << endl;

                //suppression de l'entree  dictionnaire
                dictionnairefinal.erase (it);
            }
        }
        file_ptr2.close();
        //traitement des fails
        string basis = "";
        for(uint64_t i = 0; i<nbcol2; i++){
            basis=basis+"\t0";
        }
        for(map<uint64_t, string>::iterator it = dictionnairefinal.begin(); it!=dictionnairefinal.end(); it++){
            file_output_final_ptr << it->second << basis << endl;
        }
    }


    /**
     * traitement par paquet de file2
     * **/

    for(uint64_t paquet = 0; paquet<nbpaquet2; paquet++){
        //pointeur à la bonne ligne
        ifstream file_ptr2(file_input2.c_str(), ios::in);
        for(uint64_t cpt = 0; cpt<(2+paquet*size_sample); cpt++){
            getline(file_ptr2, line);
        }

        //creation dictionnaire
        map<uint64_t, string> dictionnaire;
        for(uint64_t cpt = 0; cpt<size_sample; cpt++){
            getline(file_ptr2, line);
            vector<string> content;
            boost::split(content, line, boost::is_any_of("\t"));
            dictionnaire[strtoull(content[0].c_str(), NULL, 0)] = line;
        }
        file_ptr2.close();

        //lecture file1
        ifstream file_ptr1(file_input1.c_str(), ios::in);
        while(getline(file_ptr1, line)){
            vector<string> content;
            boost::split(content, line, boost::is_any_of("\t"));
            uint64_t idc = strtoull(content[0].c_str(), NULL, 0);
            map<uint64_t, string>::iterator it = dictionnaire.find(idc);
            if(it!=dictionnaire.end()){
                //ecriture dans output
                //file_output_final_ptr << it->second << "\t" << line << endl;
                //suppression de l'entree  dictionnaire
                dictionnaire.erase (it);
            }
        }
        file_ptr1.close();

        //traitement des fails
        string basis = "";
        for(uint64_t i = 0; i<nbcol2; i++){
            basis=basis+"\t0";
        }
        for(map<uint64_t, string>::iterator it = dictionnaire.begin(); it!=dictionnaire.end(); it++){
            vector<string> content;
            boost::split(content, it->second, boost::is_any_of("\t"));
            string interstice = "";
            for(uint64_t i =1; i<content.size(); i++){
                interstice+="\t"+content[i];
            }
            file_output_final_ptr << content[0] << basis << interstice << endl;
        }

    }
    if(true){
        //pointeur à la bonne ligne
        ifstream file_ptr2(file_input2.c_str(), ios::in);
        for(uint64_t cpt = 0; cpt<(2+nbpaquet2*size_sample); cpt++){
            getline(file_ptr2, line);
        }

        //creation dictionnaire
        map<uint64_t, string> dictionnairefinal2;
        while(getline(file_ptr2, line)){
            vector<string> content;
            boost::split(content, line, boost::is_any_of("\t"));
            dictionnairefinal2[strtoull(content[0].c_str(), NULL, 0)] = line;
        }
        file_ptr2.close();

        //lecture file1
        ifstream file_ptr1(file_input1.c_str(), ios::in);
        while(getline(file_ptr1, line)){
            vector<string> content;
            boost::split(content, line, boost::is_any_of("\t"));
            uint64_t idc = strtoull(content[0].c_str(), NULL, 0);
            map<uint64_t, string>::iterator it = dictionnairefinal2.find(idc);
            if(it!=dictionnairefinal2.end()){
                //ecriture dans output
                //file_output_final_ptr << it->second << "\t" << line << endl;
                //suppression de l'entree  dictionnaire
                dictionnairefinal2.erase (it);
            }
        }
        file_ptr1.close();

        //traitement des fails
        string basis = "";
        for(uint64_t i = 0; i<nbcol2; i++){
            basis=basis+"\t0";
        }
        for(map<uint64_t, string>::iterator it = dictionnairefinal2.begin(); it!=dictionnairefinal2.end(); it++){
            vector<string> content;
            boost::split(content, it->second, boost::is_any_of("\t"));
            string interstice = "";
            for(uint64_t i =1; i<content.size(); i++){
                interstice+="\t"+content[i];
            }
            file_output_final_ptr << content[0] << basis << interstice << endl;
        }

    }

    file_output_final_ptr.close();
    return EXIT_SUCCESS;
}
