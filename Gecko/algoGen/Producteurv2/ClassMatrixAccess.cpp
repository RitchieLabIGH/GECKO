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


using namespace std;

/*!
 * \file ClassMatrixAccess.cpp
 * \brief class to instantiate access to matrix file
 * \class Class ClassMatrixAccess
 * In case of accessing a binary file, it handles the access through direct access plus buffer.
 *
 * Currently the class handles only binary file but not ASCII version.
 *
*/



class ClassMatrixAccess{
public:
    string _self_file; //the file to access
    bool _self_isInBinary;
    int _self_modebinary; //1 is the normal binary format, 2 for splitted binary
    uint64_t _self_nKmerPerSlice; //in case of _self_modebinary==2, the number of kmers per slice ==> (_self_nbLines-2)/_self_nKmerPerSlice is the number of files
    uint64_t _self_nbColumns;
    uint64_t _self_nbLines;
    uint64_t _self_maxbuffersize;
    //vector <uint64_t> _self_realposition; //id position in the file give position in cache _self_data_x
    //vector <vector<float>> _self_data_x; // the big matrix (l)x(c) = (sample)x(kmer)
    map<uint64_t, int64_t> _self_mapScore;
    map<uint64_t, vector<float>> _self_mapValue;
    int64_t _self_map_scoreMin,_self_map_scoreMax;

    /// \brief the constructor takes as parameter a path to a file, look if the format is binary, and if so reads the number of lines and columns
    ClassMatrixAccess(string p_pathtofile, uint64_t p_maxBufferSize = 100000){
        _self_file = p_pathtofile;
        _self_isInBinary = isItABinaryFile();
        _self_maxbuffersize = p_maxBufferSize;
        _self_map_scoreMin = 1;
        _self_map_scoreMax = 0;
    }
    ClassMatrixAccess(){
    }

    /// \brief look if the file is in binary format by looking the flag (the first float that should be 412.12 (mode1) or 233.33(mode2))
    bool isItABinaryFile(){
        bool result = false;
        ifstream file_in_ptr(_self_file.c_str(), ios_base::binary);
        if(!file_in_ptr){
            cerr << "ClassMatrixAccess::isItABinaryFile Error the file " << _self_file << " does not exist" << endl;
            exit(EXIT_FAILURE);
        }
        float value;
        uint64_t nbl, nbc, nbkps;
        file_in_ptr.read(reinterpret_cast<char*>(&value), sizeof(float));

        if(round(value-421.21)==0){
            _self_modebinary = 1;
            file_in_ptr.read(reinterpret_cast<char*>(&nbl), sizeof(uint64_t));
            file_in_ptr.read(reinterpret_cast<char*>(&nbc), sizeof(uint64_t));
            _self_nbColumns = nbc;
            _self_nbLines   = nbl;
            result = true;
        }
        if(round(value-233.33)==0){
            _self_modebinary = 2;
            file_in_ptr.read(reinterpret_cast<char*>(&nbl), sizeof(uint64_t));
            file_in_ptr.read(reinterpret_cast<char*>(&nbc), sizeof(uint64_t));
            file_in_ptr.read(reinterpret_cast<char*>(&nbkps), sizeof(uint64_t));
            _self_nKmerPerSlice = nbkps;
            _self_nbColumns = nbc;
            _self_nbLines   = nbl;
            result = true;
        }
        file_in_ptr.close();
        return result;
    }

    /// \brief get back all the group ids
    vector<string> getGroups(){
        vector<string> result;
        ifstream file_in_ptr(_self_file.c_str(), ios_base::binary);
        if(!file_in_ptr){
            cerr << "ClassMatrixAccess::getGroups Error the file " << _self_file << " does not exist" << endl;
            exit(EXIT_FAILURE);
        }

        uint64_t sizeOfHeader;
        if(_self_modebinary==1){
            sizeOfHeader = sizeof(float)+2*sizeof(uint64_t)+(_self_nbColumns)*64+64;
        }
        else{
          if(_self_modebinary==2){
            sizeOfHeader = sizeof(float)+3*sizeof(uint64_t)+(_self_nbColumns)*64+64;
          }
          else{
            exit(EXIT_FAILURE);
          }
        }


        file_in_ptr.seekg( sizeOfHeader,  ios_base::beg );
        for(uint64_t i = 1; i<_self_nbColumns ; i++){
            char bufread[64];
            file_in_ptr.read(reinterpret_cast<char*>(bufread), 64);
            string group_c(bufread);
            result.push_back(group_c);
        }
        file_in_ptr.close();
        return result;
    }

    /// \brief get back all the name ids
    vector<string> getNames(){
        vector<string> result;

        ifstream file_in_ptr(_self_file.c_str(), ios_base::binary);
        if(!file_in_ptr){
            cerr << "ClassMatrixAccess::getNames Error the file " << _self_file << " does not exist" << endl;
            exit(EXIT_FAILURE);
        }

        uint64_t sizeOfHeader;
        if(_self_modebinary==1){
            sizeOfHeader = sizeof(float)+2*sizeof(uint64_t);
        }
        if(_self_modebinary==2){
            sizeOfHeader = sizeof(float)+3*sizeof(uint64_t);
        }
        file_in_ptr.seekg( sizeOfHeader,  ios_base::beg );

        for(uint64_t i = 0; i<_self_nbColumns ; i++){
            char bufread[64];
            file_in_ptr.read(reinterpret_cast<char*>(bufread), 64);
            string group_c(bufread);
            result.push_back(group_c);
        }
        file_in_ptr.close();
        return result;
    }

    /// \brief get back all the groups
    /// http://manpagesfr.free.fr/man/man2/lseek.2.html for whence (3e param of gzseek)
    vector<string> getSequenceIDs(){
        vector<string> result;
        ifstream file_in_ptr(_self_file.c_str(), ios_base::binary);
        if(!file_in_ptr){
            cerr << "ClassMatrixAccess::getSequenceIDs Error the file " << _self_file << " does not exist" << endl;
            exit(EXIT_FAILURE);
        }
        if(_self_modebinary==1){
            uint64_t sizeOfHeader = sizeof(float)+2*sizeof(uint64_t)+(_self_nbColumns)*64*2;
            uint64_t sizeOfLineWithoutKmerID = (_self_nbColumns-1)*sizeof(float);
            uint64_t value_tmp;

            file_in_ptr.seekg( sizeOfHeader,  ios_base::beg );
            for(uint64_t i = 2; i<_self_nbLines; i++){
                file_in_ptr.read(reinterpret_cast<char*>(&value_tmp), sizeof(uint64_t));
                string value_tmp_string = to_string(value_tmp);
                result.push_back(value_tmp_string);
                file_in_ptr.seekg( sizeOfLineWithoutKmerID,  ios_base::cur);
            }
        }
        if(_self_modebinary==2){
            uint64_t value_tmp;
            uint64_t sizeOfHeader = sizeof(float)+3*sizeof(uint64_t)+(_self_nbColumns)*64*2;
            file_in_ptr.seekg( sizeOfHeader,  ios_base::beg );
            for(uint64_t i = 2; i<_self_nbLines; i++){
                file_in_ptr.read(reinterpret_cast<char*>(&value_tmp), sizeof(uint64_t));
                string value_tmp_string = to_string(value_tmp);
                result.push_back(value_tmp_string);
            }
        }

        file_in_ptr.close();
        return result;
    }

    /// \brief get back a single value (indices are 0-based)
    vector<float> getValuesFromSingleKmer(uint64_t p_kmerPosition, bool updateBuf = true){
        vector<float> result;

        if(_self_mapScore.find(p_kmerPosition)!=_self_mapScore.end()){
            result = _self_mapValue[p_kmerPosition];
        }
        else{
            if(_self_modebinary==1){
                ifstream file_in_ptr(_self_file.c_str(), ios_base::binary);
                if(!file_in_ptr){
                    cerr << "ClassMatrixAccess::getValues Error the file " << _self_file << " does not exist" << endl;
                    exit(EXIT_FAILURE);
                }
                uint64_t sizeOfHeader = sizeof(float)+2*sizeof(uint64_t)+(_self_nbColumns)*64*2;
                uint64_t sizeOfLine = sizeof(uint64_t)+(_self_nbColumns-1)*sizeof(float);
                uint64_t offset = sizeOfHeader+sizeOfLine*p_kmerPosition+sizeof(uint64_t);

                file_in_ptr.seekg( offset,  ios_base::beg ); //place cursor to first line (move from file begin)
                for(uint64_t i =0; i<(_self_nbColumns-1); i++){
                    float value_tmp;
                    file_in_ptr.read(reinterpret_cast<char*>(&value_tmp), sizeof(float) );
                    result.push_back(value_tmp);
                }
                file_in_ptr.close();

                if(updateBuf){
                  _self_mapScore[p_kmerPosition] = 1;
                  _self_mapValue[p_kmerPosition] = result;
                }
            }
            if(_self_modebinary==2){
                uint64_t SliceNumber = floor((long double) p_kmerPosition/_self_nKmerPerSlice);
                string realfilePath = _self_file+"_"+to_string(SliceNumber);
                ifstream file_in_ptr(realfilePath.c_str(), ios_base::binary);
                if(!file_in_ptr){
                    cerr << "ClassMatrixAccess::getValues Error the file " << realfilePath << " does not exist" << endl;
                    exit(EXIT_FAILURE);
                }

                uint64_t sizeOfLine = (_self_nbColumns-1)*sizeof(float);
                uint64_t offset = sizeOfLine*(p_kmerPosition-SliceNumber*_self_nKmerPerSlice);

                file_in_ptr.seekg( offset,  ios_base::beg ); //place cursor to first line (move from file begin)
                for(uint64_t i =0; i<(_self_nbColumns-1); i++){
                    float value_tmp;
                    file_in_ptr.read(reinterpret_cast<char*>(&value_tmp), sizeof(float) );
                    result.push_back(value_tmp);
                }
                file_in_ptr.close();

                if(updateBuf){
                  _self_mapScore[p_kmerPosition] = 1;
                  _self_mapValue[p_kmerPosition] = result;
                }
            }

        }
        return result;
    }

    ///\brief load values of the kmers that are stored in table format (see Population) into buffer
    // beginIndiv and endInviv indice begin and end of subtable to load ( use for computeHDinScore)
    // defaults values = 0 or non affected will load the entire table
    void getValuesFromTable(vector <vector<uint64_t>> p_table,uint64_t beginIndiv=0,uint64_t endInviv=0){


        ////////////////////////////////////////////
        //load the kmers from the table in parameter
        uint64_t nbIndivTable = p_table.size();
        uint64_t nbKmerPerIndivTable = p_table[0].size();
        //For default value endInviv load all individuals from the table
        if (endInviv == 0){
            endInviv= nbIndivTable;
        }
        map<uint64_t, bool> mapKmerSeen; //if the kmer has been found this time


        if(_self_modebinary==1){
            ifstream file_in_ptr(_self_file.c_str(), ios_base::binary);
            if(!file_in_ptr){
                cerr << "ClassMatrixAccess::getValuesFromTable Error the file " << _self_file << " does not exist" << endl;
                exit(EXIT_FAILURE);
            }
            uint64_t sizeOfHeader = sizeof(float)+2*sizeof(uint64_t)+(_self_nbColumns)*64*2;
            uint64_t sizeOfLine = sizeof(uint64_t)+(_self_nbColumns-1)*sizeof(float);
            vector<uint64_t> kmerToImport;

            for(uint64_t i_indiv = beginIndiv; i_indiv<endInviv; i_indiv++){
                for(uint64_t i_kmer = 0; i_kmer<nbKmerPerIndivTable; i_kmer++){

                    uint64_t Kkmer = p_table[i_indiv][i_kmer];
                    if(_self_mapScore.find(Kkmer)!=_self_mapScore.end()){
                        //the kmer is already in memory
                        _self_mapScore[Kkmer]++;
                        mapKmerSeen[Kkmer]=true;
                    }
                    else{
                        //the kmer is new, get the data
                        kmerToImport.push_back(Kkmer);
                    }
                }
            }
            sort(kmerToImport.begin(),kmerToImport.end());

            for(uint64_t ik = 0; ik<kmerToImport.size(); ik++){
                uint64_t Kkmer= kmerToImport[ik];
                uint64_t offset = sizeOfHeader+sizeOfLine*Kkmer+sizeof(uint64_t);
                file_in_ptr.seekg( offset,  ios_base::beg ); //place cursor to first line (move from file begin)
                vector<float> vect_tmp;
                for(uint64_t i =0; i<(_self_nbColumns-1); i++){
                    float value_tmp;
                    file_in_ptr.read(reinterpret_cast<char*>(&value_tmp), sizeof(float) );
                    vect_tmp.push_back(value_tmp);
                }
                _self_mapValue[Kkmer] = vect_tmp;
                _self_mapScore[Kkmer]=1;
                mapKmerSeen[Kkmer]=true;
            }
            file_in_ptr.close();
        }

        if(_self_modebinary==2){
            uint64_t sizeOfLine = (_self_nbColumns-1)*sizeof(float);
            map<uint64_t, vector<uint64_t>> kmerPerSlice;

            for(uint64_t i_indiv = beginIndiv; i_indiv<endInviv; i_indiv++){
                for(uint64_t i_kmer = 0; i_kmer<nbKmerPerIndivTable; i_kmer++){

                    uint64_t Kkmer = p_table[i_indiv][i_kmer];
                    if(_self_mapScore.find(Kkmer)!=_self_mapScore.end()){
                        //the kmer is already in memory
                        _self_mapScore[Kkmer]++;
                        mapKmerSeen[Kkmer]=true;
                    }
                    else{
                        //the kmer is new, get the data
                        uint64_t SliceNumber = floor((long double) Kkmer/_self_nKmerPerSlice);
                        if(kmerPerSlice.find(SliceNumber)==kmerPerSlice.end()){
                            vector<uint64_t> tmp;
                            kmerPerSlice[SliceNumber] = tmp;
                        }
                        kmerPerSlice[SliceNumber].push_back(Kkmer);

                        vector<float> tmpvectvalues((_self_nbColumns-1),0);
                        _self_mapValue[Kkmer] = tmpvectvalues;
                        _self_mapScore[Kkmer] = 1;
                        mapKmerSeen[Kkmer]=true;
                    }
                }
            }

            vector<uint64_t> vectorOfSliceNumber;
            for(auto itslice = kmerPerSlice.begin(); itslice!=kmerPerSlice.end(); itslice++){
                vectorOfSliceNumber.push_back(itslice->first);
            }

            uint64_t vsn = vectorOfSliceNumber.size();
            #pragma omp parallel for
            for(uint64_t islice = 0; islice<vsn; islice++){
                uint64_t SliceNumber = vectorOfSliceNumber[islice];
                sort(kmerPerSlice[SliceNumber].begin(), kmerPerSlice[SliceNumber].end());
                string realfilePath = _self_file+"_"+to_string(SliceNumber);

                ifstream file_in_ptr(realfilePath.c_str(), ios_base::binary);
                if(!file_in_ptr){
                    cerr << "ClassMatrixAccess::getValues Error the file " << realfilePath << " does not exist" << endl;
                    exit(EXIT_FAILURE);
                }
                for(uint64_t ik = 0; ik<(kmerPerSlice[SliceNumber].size()); ik++){
                    uint64_t Kkmer = kmerPerSlice[SliceNumber][ik];
                    if(ik==0){
                        uint64_t offset = sizeOfLine*(Kkmer-(SliceNumber*_self_nKmerPerSlice));
                        file_in_ptr.seekg( offset,  ios_base::beg ); //place cursor to first line (move from file begin)
                    }
                    else{
                        uint64_t KkmerPrev = (kmerPerSlice[SliceNumber][ik-1])+1;
                        uint64_t offset = sizeOfLine*(Kkmer-KkmerPrev);
                        file_in_ptr.seekg( offset,  ios_base::cur ); //place cursor to first line (move from file begin)
                    }
                    vector<float> vect_tmp;
                    for(uint64_t i =0; i<(_self_nbColumns-1); i++){
                        float value_tmp;
                        file_in_ptr.read(reinterpret_cast<char*>(&value_tmp), sizeof(float) );
                        vect_tmp.push_back(value_tmp);
                    }
                    _self_mapValue[Kkmer] = vect_tmp;

                }
                file_in_ptr.close();
            }

        }

        /////////////////////////////////////////
        //penality to kmers that have disappeared
        for(auto it_kmerref = _self_mapScore.begin(); it_kmerref!=_self_mapScore.end(); it_kmerref++){
            uint64_t Kkmer = it_kmerref->first;
            if(mapKmerSeen.find(Kkmer)==mapKmerSeen.end())
                _self_mapScore[Kkmer]=0;
        }

        /////////////////////////////////////////
        //remove kmers if memory is exceeded
        cerr << _self_mapScore.size() << " elements in buffer" << endl;
        if(_self_mapScore.size()>_self_maxbuffersize){
            int64_t diff = _self_mapScore.size() - _self_maxbuffersize;
            cerr << diff << " to be removed" << endl;

            while(diff>0){
                //cerr << diff << " to be removed" << endl;

                auto it_ms_init = _self_mapScore.begin();
                _self_map_scoreMax = it_ms_init->second;
                _self_map_scoreMin = it_ms_init->second;
                for(auto it_ms = _self_mapScore.begin(); it_ms!=_self_mapScore.end(); it_ms++){
                    if(it_ms->second > _self_map_scoreMax)
                        _self_map_scoreMax = it_ms->second;
                    if(it_ms->second < _self_map_scoreMin)
                        _self_map_scoreMin = it_ms->second;
                }
                /*cerr << _self_map_scoreMin << " as value for _self_map_scoreMin" << endl;
                cerr << _self_map_scoreMax << " as value for _self_map_scoreMax" << endl;*/

                vector<uint64_t> KmerToRemove;

                for(auto it_kmerref = _self_mapScore.begin(); it_kmerref!=_self_mapScore.end(); it_kmerref++){
                    uint64_t Kkmer = it_kmerref->first;
                    int64_t Skmer = it_kmerref->second;

                    if(diff>0){
                        if(Skmer==_self_map_scoreMin){
                            if(mapKmerSeen.find(Kkmer)==mapKmerSeen.end()){
                            //_self_mapScore.erase(Kkmer);
                            //_self_mapValue.erase(Kkmer);
                                KmerToRemove.push_back(Kkmer);
                                diff--;
                            }
                        }
                    }
                }

                //cerr << "delete batch of " << KmerToRemove.size() << " elements" <<endl;
                for(uint64_t ki = 0; ki<KmerToRemove.size(); ki++){
                    _self_mapScore.erase(KmerToRemove[ki]);
                    _self_mapValue.erase(KmerToRemove[ki]);
                }
                KmerToRemove.erase (KmerToRemove.begin(),KmerToRemove.end());

            }
            //cerr << "end remove kmers" << endl;
        }
    }



    float getValue(uint64_t p_samplePosition, uint64_t p_kmerPosition, string debug=""){
        if(_self_mapScore.find(p_kmerPosition)!=_self_mapScore.end()){
            //cout << "In memory" << endl;
            return _self_mapValue[p_kmerPosition][p_samplePosition];
        }
        else{
            cout << debug << ": need to preload" << endl;
            vector<float> result = getValuesFromSingleKmer(p_kmerPosition);
            return result[p_samplePosition];
        }
    }

    //! \fn uint64_t getIndice(string kmer, int64_t rank=0)
    //! \brief transform a sequence into an integer value that corresponds to {A,C,T,G} in base-4 algebra
    //! \return (uint64_t) indice that corresponds to sequence
    //!
    uint64_t getIndice(string kmer, int64_t rank=0){
        if(kmer.size()==rank){
            return 0;
        }
        else{
            char character_c = kmer[kmer.size()-rank-1];
            if(character_c=='A'){
                return getIndice(kmer, rank+1)+0;
            }
            else{
                if(character_c=='C'){
                    return getIndice(kmer, rank+1)+1*powl(4,rank);
                }
                else{
                    if(character_c=='T'){
                        return getIndice(kmer, rank+1)+2*powl(4,rank);
                    }
                    else{
                        return getIndice(kmer, rank+1)+3*powl(4,rank);
                    }
                }
            }
        }
    }

    //! transform sequence ID to real nucleotidic sequence
    string idseq2fasta(uint64_t realid, uint64_t nbBasePerKmer){
        string result = "";
        int64_t k = nbBasePerKmer - 1;
        vector<string> letters;
        letters.push_back("A");
        letters.push_back("C");
        letters.push_back("T");
        letters.push_back("G");

        while (k >= 0){
            uint64_t kmin = pow(4, k);
            if(kmin <= realid){
                uint64_t alpha = 3;
                bool statealpha = true;

                while (statealpha){
                    if (realid >= (alpha * kmin)){
                        realid -= (alpha * kmin);
                        statealpha = false;
                        result += letters[alpha];
                    }
                    else{
                        alpha -= 1;
                    }
                }
            }
            else{
                result += letters[0];
            }
            k -= 1;
        }
        return result;
    }


};


/*
int main(){
    chrono::time_point<std::chrono::system_clock> start1, end1, start2, end2, start3, end3, start4, end4;
    //string file = "/home/aubin/ritchie/gecko/TCGA/TCGA_filtered_realCounts.bin";
    string file = "/home/aubin/ritchie/gecko/peripheralBlood/PRJNA391912_norm_discret_Clean_realCounts.bin";

    string file2 = "/home/aubin/ritchie/gecko/LargeScaleCLLGenome/05-ML/dev/importGECKO/bin/src/test.binary.old";
    ClassMatrixAccess testNew = ClassMatrixAccess(file,150000);
    //ClassMatrixAccess testOld = ClassMatrixAccess(file2);

    if(true){
        cout << testNew._self_isInBinary << "\t" << testNew._self_nbLines << "\t" << testNew._self_nbColumns << endl;
        //cout << testOld._self_isInBinary << "\t" << testOld._self_nbLines << "\t" << testOld._self_nbColumns << endl;
    }

    // groups
    if(true){
        vector<string> testgroups = testNew.getGroups();
        cout << "found " << testgroups.size() << " groups" << endl;
        for(int i =0; i<testgroups.size(); i++){
            cout << testgroups[i] << endl;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sequences ID
    if(true){
        vector<string> testseq = testNew.getSequenceIDs();
        cout << "found " << testseq.size() << " sequences" << endl;

        for(int i = 0; i<8 ;i++){
            cout << testseq[i] << endl;
        }
    }

    return 0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // getValues
    if(false){
        vector<float> v1 = testNew.getValuesFromSingleKmer(0);
        vector<float> v2 = testNew.getValuesFromSingleKmer(1);

        for(uint64_t k =0; k<5; k++){
            cout << v1[k] << " , " << v2[k] << endl;
        }
    }
    if(true){

        cout << "GO !" << endl;
        for(uint64_t i = 80000; i>0; i--){
            vector<float> v3 = testNew.getValuesFromSingleKmer(i);
        }

        cout << "GO !" << endl;
        for(uint64_t i = 80000; i>0; i--){
            vector<float> v3 = testNew.getValuesFromSingleKmer(i);
        }


        vector<vector<uint64_t>> table1;
        uint64_t k=100000;
        for(uint64_t i =0; i<1000; i++){
            vector<uint64_t> tmp;
            for(uint64_t j =0; j<50; j++){
                tmp.push_back(k);
                k+=5;
            }
            table1.push_back(tmp);
        }


        vector<vector<uint64_t>> table2;
        k=10000;
        for(uint64_t i =0; i<1000; i++){
            table2.push_back(table1[i]);
            i++;
            vector<uint64_t> tmp;
            for(uint64_t j =0; j<50; j++){
                tmp.push_back(k);
                k+=2;
            }
            table2.push_back(tmp);
        }

        vector<vector<uint64_t>> table3;
        k=20000;
        for(uint64_t i =0; i<1000; i++){
            vector<uint64_t> tmp;
            for(uint64_t j =0; j<50; j++){
                tmp.push_back(k);
                k+=20;
            }
            table3.push_back(tmp);
            table3.push_back(table1[i]);
            table3.push_back(table2[i]);
            table3.push_back(table1[i+1]);
            table3.push_back(table2[i+1]);
            i+=4;
        }

        vector<vector<uint64_t>> table4;
        k=100000;
        for(uint64_t i =0; i<1000; i++){
            vector<uint64_t> tmp;
            for(uint64_t j =0; j<50; j++){
                tmp.push_back(k);
                k+=1;
            }
            table4.push_back(tmp);
        }

        start3 = chrono::system_clock::now();
        cout << "GO !" << endl;
        testNew.getValuesFromTable(table1);
        end3 = chrono::system_clock::now();

        start4 = chrono::system_clock::now();
        cout << "GO !" << endl;
        testNew.getValuesFromTable(table2);
        end4 = chrono::system_clock::now();

        start1 = chrono::system_clock::now();
        cout << "GO !" << endl;
        testNew.getValuesFromTable(table3);
        end1 = chrono::system_clock::now();

        start2 = chrono::system_clock::now();
        cout << "GO !" << endl;
        testNew.getValuesFromTable(table4);
        end2 = chrono::system_clock::now();
    }

    int elapsed_seconds1 = chrono::duration_cast<std::chrono::milliseconds> (end1-start1).count();
    int elapsed_seconds2 = chrono::duration_cast<std::chrono::milliseconds> (end2-start2).count();
    cout << "time1 = " << elapsed_seconds1 << endl;
    cout << "time2 = " << elapsed_seconds2 << endl;

    int elapsed_seconds3 = chrono::duration_cast<std::chrono::milliseconds> (end3-start3).count();
    int elapsed_seconds4 = chrono::duration_cast<std::chrono::milliseconds> (end4-start4).count();
    cout << "time3 = " << elapsed_seconds1 << endl;
    cout << "time4 = " << elapsed_seconds2 << endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // getValue
    if(false){
        cout << testNew.getValue(0,1) << endl;
        cout << testNew.getValue(1,1) << endl;
        cout << testNew.getValue(2,1) << endl;
        cout << testNew.getValue(3,1) << endl;
        cout << testNew.getValue(4,1) << endl;
    }



    return EXIT_SUCCESS;
}*/
