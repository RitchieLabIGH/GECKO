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
#include <boost/algorithm/cxx11/is_permutation.hpp>
#include <stdlib.h>
#include <time.h>
#include <numeric>


using namespace std;

class ClassGoodIndividuals{
public:

    vector <vector<long double>> _self_Network; // the big square matrix
    vector <string> _self_SequenceID; //the sequence ID that refers to matrix
    map<string, uint64_t> _self_kmer_position; //the position of a kmer in _self_SequenceID

    vector<string> _self_indiv_seq; //the lists of kmers 
    vector<float> _self_indiv_scores; //the corresponding score
    map<string, uint64_t> _self_indiv_seq_occ; //the number of occurences of the list of kmers
    map<string, bool> _self_indiv_seq_state; //if the sequence can be taken in account

    vector<uint64_t> _self_indiv_iter; //for new kmer per iteration
    vector<uint64_t> _self_newKmer; //for new kmer per iteration

    string _self_file; //the csv file aka goodCandidates.csv


    ClassGoodIndividuals(string p_file){
        _self_file = p_file;
    }
    
    /*!
    * \brief import the csv file 
    * p_score can be :
    *           - test
    *           - outer
    *           - prod : test*outer
    *           - adjust : test*outer-overfit
    * p_filterocc : the individuals are considered only if they appears at least filterOcc times
    * */
    void importCSVfile(string p_score, uint64_t p_filterOcc = 1, float p_valtestmin = 0.8, bool getBesOcc=false){
        ifstream file_ai(_self_file.c_str(), ios::in);
        if (!file_ai) {
            cerr << "Error while opening " << _self_file << " in read mode" << endl;
            exit(EXIT_FAILURE);
        }
        string line;
        vector<string> content;

        string iter="-1";//pour ne garder que la meilleure occurence
        string iter_c; //pour ne garder que la meilleure occurence
        float valIterMax; //pour ne garder que la meilleure occurence
        string combIterMax; //pour ne garder que la meilleure occurence

        //first we store all individuals with scores and ranked IDs
        //cerr << "ClassGoodIndividuals::importCSVfile:: loading data " << endl;
        while(getline(file_ai, line)){
            boost::split(content, line, boost::is_any_of("\t"));
            string seq_reformat = reformatSequences(content[1]); //sort kmers in the list in ascending order
            _self_indiv_seq.push_back(seq_reformat);
            _self_indiv_iter.push_back(stoull(content[0]));

            float testValue = stof(content[2]);
            float outerValue = stof(content[3]);
            _self_indiv_scores.push_back( getScores(testValue, outerValue, p_score) );
            
            if(_self_indiv_seq_occ.find(seq_reformat)==_self_indiv_seq_occ.end()){
                _self_indiv_seq_occ[seq_reformat] = 1;
                _self_indiv_seq_state[seq_reformat] = true;
            }
            else{
                _self_indiv_seq_occ[seq_reformat]+=1;
            }
            if(testValue<p_valtestmin){
                _self_indiv_seq_state[seq_reformat] = false;
            }
            if(getBesOcc){
                iter_c = content[0];
                if(iter_c.compare(iter)!=0){
                    iter = iter_c;
                    valIterMax = testValue;
                    combIterMax = seq_reformat;
                }
                else{
                    if(valIterMax<testValue){
                        valIterMax = testValue;
                        _self_indiv_seq_state[combIterMax] = false;
                        combIterMax = seq_reformat;
                    }
                }
            }
            
        }
        file_ai.close();

        //look how many times each individual appears, and mark if under filter value
        //cerr << "ClassGoodIndividuals::importCSVfile:: filtering on occurences " << endl;
        if (p_filterOcc>1){
            for(uint64_t i = 0; i<_self_indiv_seq.size(); i++){
                if(_self_indiv_seq_occ[_self_indiv_seq[i]]<p_filterOcc){
                    _self_indiv_seq_state[_self_indiv_seq[i]] = false;
                }
            }
        }



        //create the data structure
        //cerr << "ClassGoodIndividuals::computeNetwork:: creating data structure " << endl;
        uint64_t k = 0;
        //store the kmers that compose the good individuals
        for(uint64_t i=0; i<_self_indiv_seq.size(); i++){
            if(_self_indiv_seq_state[_self_indiv_seq[i]]){
                
                boost::split(content, _self_indiv_seq[i], boost::is_any_of(","));
                for(uint64_t j = 0; j<content.size(); j++){
                    
                    if(_self_kmer_position.find(content[j])==_self_kmer_position.end()){
                        _self_SequenceID.push_back(content[j]);
                        _self_kmer_position[content[j]] = k;
                        k++;
                    }
                }
            }
        }

        
    }



    /*!
    * \brief compute the network
    * the method can be :
    *           - max : for each relationship only the maxima is taken
    *           - average : for each relationship the mean is computed
    *           - sum : for each relationship the sum is computed
    *           - singleAverage : each individuals are unique, the average is then computed
    *           - singleSum : each individuals are unique, the cumulation is then computed
    * 
    * */
    void computeNetwork(string p_method, vector<string> SequenceID){

        cerr <<  SequenceID.size() << " kmers to be loaded" << endl;
        map<string, uint64_t> kmer_position;
        _self_Network.erase (_self_Network.begin(),_self_Network.end());

        for(uint64_t i=0; i<SequenceID.size(); i++){
            kmer_position[SequenceID[i]] = i;
            vector<long double> tmp(SequenceID.size(), 0.0);
            _self_Network.push_back(tmp);
        }
        for(uint64_t i=0; i<SequenceID.size(); i++){
            for(uint64_t j=0; j<SequenceID.size(); j++){
                _self_Network[i][j] = 0.0;
            }
        }


        if(p_method.compare("max")==0){
            for(uint64_t i=0; i<_self_indiv_seq.size(); i++){
                if(_self_indiv_seq_state[_self_indiv_seq[i]]){
                    
                    vector<vector<string>> couples = getCouplesOfIndex(_self_indiv_seq[i]);
                    for(uint64_t c = 0 ; c<couples.size(); c++){
                        
                        if((kmer_position.find(couples[c][0])!=kmer_position.end()) and (kmer_position.find(couples[c][1])!=kmer_position.end())){
                            uint64_t pos1 = kmer_position[couples[c][0]];
                            uint64_t pos2 = kmer_position[couples[c][1]];
                            if(_self_Network[pos1][pos2] < _self_indiv_scores[i]){
                                _self_Network[pos1][pos2] = _self_indiv_scores[i];
                                _self_Network[pos2][pos1] = _self_indiv_scores[i];
                            }
                        }
                    }

                }
            }
        
        }

        if(p_method.compare("average")==0){
            vector<vector<uint64_t>> Network_counts;
            for(uint64_t i=0; i<SequenceID.size(); i++){
                vector<uint64_t> tmp(SequenceID.size(), 0);
                Network_counts.push_back(tmp);
            }

            for(uint64_t i=0; i<_self_indiv_seq.size(); i++){
                if(_self_indiv_seq_state[_self_indiv_seq[i]]){
                    vector<vector<string>> couples = getCouplesOfIndex(_self_indiv_seq[i]);
                    
                    for(uint64_t c = 0 ; c<couples.size(); c++){
                        if((kmer_position.find(couples[c][0])!=kmer_position.end()) and (kmer_position.find(couples[c][1])!=kmer_position.end())){
                            uint64_t pos1 = kmer_position[couples[c][0]];
                            uint64_t pos2 = kmer_position[couples[c][1]];
                            
                            _self_Network[pos1][pos2]+=_self_indiv_scores[i];
                            _self_Network[pos2][pos1]+=_self_indiv_scores[i];
                            Network_counts[pos1][pos2]+=1;
                            Network_counts[pos2][pos1]+=1;
                        }
                    }
                }
            }
            for(uint64_t i=0; i<SequenceID.size(); i++){
                for(uint64_t j=i+1; j<SequenceID.size(); j++){
                    if(Network_counts[i][j]!=0){
                        _self_Network[i][j]/=Network_counts[i][j];
                        _self_Network[j][i]/=Network_counts[j][i];
                    }
                }
            }
        
        }

        if(p_method.compare("sum")==0){
            for(uint64_t i=0; i<_self_indiv_seq.size(); i++){
                if(_self_indiv_seq_state[_self_indiv_seq[i]]){
                    vector<vector<string>> couples = getCouplesOfIndex(_self_indiv_seq[i]);
                    for(uint64_t c = 0 ; c<couples.size(); c++){
                        if((kmer_position.find(couples[c][0])!=kmer_position.end()) and (kmer_position.find(couples[c][1])!=kmer_position.end())){
                            uint64_t pos1 = kmer_position[couples[c][0]];
                            uint64_t pos2 = kmer_position[couples[c][1]];
                            
                            _self_Network[pos1][pos2]+=_self_indiv_scores[i];
                            _self_Network[pos2][pos1]+=_self_indiv_scores[i];
                        }
                    }
                }
            }
            
        }

        if(p_method.compare("singleAverage")==0){
            map<string,bool> indiv_seq_state (_self_indiv_seq_state);

            vector<vector<uint64_t>> Network_counts;
            for(uint64_t i=0; i<SequenceID.size(); i++){
                vector<uint64_t> tmp(SequenceID.size(), 0);
                Network_counts.push_back(tmp);
            }

            for(uint64_t i=0; i<_self_indiv_seq.size(); i++){
                if(indiv_seq_state[_self_indiv_seq[i]]){
                    indiv_seq_state[_self_indiv_seq[i]] = false;
                    vector<vector<string>> couples = getCouplesOfIndex(_self_indiv_seq[i]);
                    
                    for(uint64_t c = 0 ; c<couples.size(); c++){
                        if((kmer_position.find(couples[c][0])!=kmer_position.end()) and (kmer_position.find(couples[c][1])!=kmer_position.end())){
                            uint64_t pos1 = kmer_position[couples[c][0]];
                            uint64_t pos2 = kmer_position[couples[c][1]];
                            
                            _self_Network[pos1][pos2]+=_self_indiv_scores[i];
                            _self_Network[pos2][pos1]+=_self_indiv_scores[i];
                            Network_counts[pos1][pos2]+=1;
                            Network_counts[pos2][pos1]+=1;
                        }
                    }
                }
            }
            for(uint64_t i=0; i<SequenceID.size(); i++){
                for(uint64_t j=i+1; j<SequenceID.size(); j++){
                    if(Network_counts[i][j]!=0){
                        _self_Network[i][j]/=Network_counts[i][j];
                        _self_Network[j][i]/=Network_counts[j][i];
                    }
                }
            }
        }

        if(p_method.compare("singleSum")==0){
            map<string,bool> indiv_seq_state (_self_indiv_seq_state);
            for(uint64_t i=0; i<_self_indiv_seq.size(); i++){
                if(indiv_seq_state[_self_indiv_seq[i]]){
                    indiv_seq_state[_self_indiv_seq[i]] = false;
                    vector<vector<string>> couples = getCouplesOfIndex(_self_indiv_seq[i]);
                    
                    for(uint64_t c = 0 ; c<couples.size(); c++){
                        if((kmer_position.find(couples[c][0])!=kmer_position.end()) and (kmer_position.find(couples[c][1])!=kmer_position.end())){
                            uint64_t pos1 = kmer_position[couples[c][0]];
                            uint64_t pos2 = kmer_position[couples[c][1]];
                            
                            _self_Network[pos1][pos2]+=_self_indiv_scores[i];
                            _self_Network[pos2][pos1]+=_self_indiv_scores[i];
                        }
                    }
                }
            }
        }

    }

    /*!
    * \brief get all the scores from test and outer scores (see description of importCSVfile function)
    * used in importCSVfile
    * */
    float getScores(float p_test, float p_outer, string p_score){
        float overfit = p_test-p_outer;
        if(overfit<0){
            overfit = 0;
        }
        if(p_score.compare("test")==0)
            return p_test;
        if(p_score.compare("outer")==0)
            return p_outer;
        if(p_score.compare("pro")==0)
            return p_test*p_outer;
        if(p_score.compare("adjust")==0)
        return p_test*p_outer-overfit;
        return 0.0;
    }

    /*!
    * \brief reformat a list of IDs in ascending order
    * used in importCSVfile
    * */
    string reformatSequences(string p_seq){
        vector<string> contentSeq;
        boost::split(contentSeq, p_seq, boost::is_any_of(","));
        sort(contentSeq.begin(), contentSeq.end());
        string result = contentSeq[0];
        for(uint i = 1; i<contentSeq.size(); i++){
            result+=","+contentSeq[i];
        }
        return result;
    }

    /*!
    * \brief get all couples of sequences 
    * used in computeNetwork
    * */
    vector<vector<string>> getCouplesOfIndex(string p_seq){
        //ptrdiff_t pos = distance(Names.begin(), find(Names.begin(), Names.end(), old_name_));

        vector<vector<string>> result;
        vector<string> contentSeq;
        boost::split(contentSeq, p_seq, boost::is_any_of(","));
        
        for(uint i = 0; i<contentSeq.size(); i++){
            string id1 = contentSeq[i];
            for(uint j = i+1; j<contentSeq.size(); j++){

                string id2 = contentSeq[j];
                vector<string> tmp;
                tmp.push_back(id1);
                tmp.push_back(id2);
                result.push_back(tmp);
            }
        }
        return result;
    }

    /*!
    * \brief output matrix
    * */
    void outputMatrix(string p_file, vector<string> SequenceID){
        
        ofstream file_ao(p_file.c_str(), ios::out);
        if(file_ao){
            file_ao << SequenceID[0];
            for(uint64_t i=1; i<SequenceID.size(); i++){
                file_ao << "," << SequenceID[i];
            }
            file_ao << endl;
            for(uint64_t i=0; i<SequenceID.size(); i++){
                file_ao << SequenceID[i];
                for(uint64_t j = 0; j<SequenceID.size(); j++){
                    file_ao << "," << _self_Network[i][j];
                }
                file_ao << endl;
            }
            file_ao.close();
        }
        else{
            cerr << "ClassGoodIndividuals::outputMatrix Error while opening" << p_file << " in write mode" << endl;
        }

    }

    /*!
    * \brief cross _self_SequenceID with other _SequenceID
    * method can be : merge/intersect
    * */
    vector<string> crossSequencesIDs(vector<string> otherSequenceID, string p_method="merge"){
        vector<string> result;
        if(p_method.compare("merge") == 0 ){
            map<string, bool> tmpdict;
            for(uint64_t ipos = 0; ipos<_self_SequenceID.size(); ipos++){
                if(tmpdict.find(_self_SequenceID[ipos])==tmpdict.end()){
                   tmpdict[_self_SequenceID[ipos]] = true;
                   result.push_back(_self_SequenceID[ipos]);
                }
            }
            for(uint64_t ipos = 0; ipos<otherSequenceID.size(); ipos++){
                if(tmpdict.find(otherSequenceID[ipos])==tmpdict.end()){
                   tmpdict[otherSequenceID[ipos]] = true;
                   result.push_back(otherSequenceID[ipos]);
                }
            }
        }
        if(p_method.compare("intersect") == 0 ){
            map<string, bool> tmpdict;
            for(uint64_t ipos = 0; ipos<_self_SequenceID.size(); ipos++){
                tmpdict[_self_SequenceID[ipos]] = true;
            }
            for(uint64_t ipos = 0; ipos<otherSequenceID.size(); ipos++){
                if(tmpdict.find(otherSequenceID[ipos])!=tmpdict.end()){
                   result.push_back(otherSequenceID[ipos]);
                }
            }
        }
        return result;
    }

    /*!
    * \brief count the overlap between _self_SequenceID and other SequenceID
    * */
    vector<uint64_t> overlapSequencesIDs(vector<string> otherSequenceID){
        vector<uint64_t> result(3,0);
        result[0] = _self_SequenceID.size();
        result[2] = otherSequenceID.size();
                
        map<string, bool> tmpdict;
        for(uint64_t ipos = 0; ipos<_self_SequenceID.size(); ipos++){
            tmpdict[_self_SequenceID[ipos]] = true;
        }
        for(uint64_t ipos = 0; ipos<otherSequenceID.size(); ipos++){
            if(tmpdict.find(otherSequenceID[ipos])!=tmpdict.end()){
                result[1]+=1;
            }
        }
        
        return result;
    }

    /*!
    * \brief get histogram of counts of kmers
    * */
    map<string, uint64_t> getHistogram(){
        map<string, uint64_t> result;
        for(uint64_t i=0; i<_self_indiv_seq.size(); i++){
            if(_self_indiv_seq_state[_self_indiv_seq[i]]){
                vector<string> content;
                boost::split(content, _self_indiv_seq[i], boost::is_any_of(","));
                for(uint64_t j = 0; j<content.size(); j++){
                    
                    if(result.find(content[j])==result.end()){
                        result[content[j]] = 0;
                    }
                    
                    result[content[j]]+= 1;
                    
                }
            }
        }
        /*for(auto it = result.begin(); it!=result.end(); it++){
            cerr << it->first << "\t" << it->second << endl;
        }*/
        return result;
    }

    /*!
    * \brief get the sum of the lines in the network for each kmer
    * !!! work only if the network has been computed for the _self_SequenceID
    * */
    map<string, long double> getHistogramOfContacts(){
        map<string, long double> result;
        cerr << _self_SequenceID.size() << " elements " << endl;
        for(uint64_t iseq = 0; iseq<_self_SequenceID.size(); iseq++){
            result[_self_SequenceID[iseq]]=0;

            for(uint64_t ibinome = 0; ibinome < _self_SequenceID.size(); ibinome++){
                result[_self_SequenceID[iseq]]+=_self_Network[iseq][ibinome];
            }
            //cerr << _self_SequenceID[iseq] << " : " << result[_self_SequenceID[iseq]] << endl;
        }
        return result;
    }

    /*!
    * \brief get curve number of new kmers per iteration in a file
    * */
    void getNewKmerPerIteration(string p_file){
        map<string,bool> tmpmap;
        map<uint64_t, uint64_t> result;
        uint64_t itermin = _self_indiv_iter[0];
        uint64_t itermax=0;
        
        for(uint64_t i = 0; i<_self_indiv_iter.size(); i++){
            
            if(result.find(_self_indiv_iter[i])==result.end()){
                result[_self_indiv_iter[i]] = 0;
            }

            if(_self_indiv_seq_state[_self_indiv_seq[i]]){
                vector<string> content;
                boost::split(content, _self_indiv_seq[i], boost::is_any_of(","));
                for(uint64_t j = 0; j<content.size(); j++){
                    if(tmpmap.find(content[j])==tmpmap.end()){
                       tmpmap[content[j]]=true;
                       result[_self_indiv_iter[i]]+=1;
                    }
                }
            }

            if(_self_indiv_iter[i]>itermax){
                itermax = _self_indiv_iter[i];
            }
            if(_self_indiv_iter[i]<itermin){
                itermin = _self_indiv_iter[i];
            }
        }

        ofstream file_ao(p_file.c_str(), ios::out);
        if(file_ao){
            file_ao << "iteration\tNewKmers" << endl;
            for(uint64_t p = itermin; p<=itermax; p++){
                if(result.find(p)!=result.end()){
                    file_ao << p << "\t" << result[p] << endl;
                }
            }
            file_ao.close();
        }
        else{
            cerr << "ClassGoodIndividuals::getNewKmerPerIteration Error while opening" << p_file << " in write mode" << endl;
        }
    }


    /*!
    * \brief compute if plateau is obtained
    * first computes the curve of new kmers per iteration
    * then interpolation
    * */
    bool getConvergence(uint64_t NBINS, long double threshold){
        
        map<string,bool> tmpmap;
        map<uint64_t, uint64_t> result;
        uint64_t itermin = _self_indiv_iter[0];
        uint64_t itermax=0;
        
        for(uint64_t i = 0; i<_self_indiv_iter.size(); i++){
            
            if(result.find(_self_indiv_iter[i])==result.end()){
                result[_self_indiv_iter[i]] = 0;
            }

            if(_self_indiv_seq_state[_self_indiv_seq[i]]){
                vector<string> content;
                boost::split(content, _self_indiv_seq[i], boost::is_any_of(","));
                for(uint64_t j = 0; j<content.size(); j++){
                    if(tmpmap.find(content[j])==tmpmap.end()){
                       tmpmap[content[j]]=true;
                       result[_self_indiv_iter[i]]+=1;
                    }
                }
            }

            if(_self_indiv_iter[i]>itermax){
                itermax = _self_indiv_iter[i];
            }
            if(_self_indiv_iter[i]<itermin){
                itermin = _self_indiv_iter[i];
            }
        }

        /*for(uint64_t p = itermin; p<=itermax; p++){
            if(result.find(p)!=result.end()){
                file_ao << p << "\t" << result[p] << endl;
            }
        }*/
        
        //interpolation in NBINS (interpolation by the maximum value)
        uint64_t iterstep = (long double) (itermax/NBINS);
        vector<long double> IX (NBINS,0);
        vector<long double> IY (NBINS,0);

        for(uint64_t i = itermin; i<=itermax; i++){
            if(result.find(i)!=result.end()){
                uint64_t iter_current = i;
                uint64_t value_current =  result[i];
                uint64_t bin_c = floor((long double)iter_current/(long double)iterstep);
                if(value_current>IY[bin_c])
                    IY[bin_c]=(long double)value_current;
            }
        }

        for(uint64_t i_x = 0; i_x<NBINS; i_x++){
            IX[i_x] = i_x*iterstep+iterstep/2;
            //cerr << IX[i_x] << "\t" << IY[i_x] << endl;
        }

        //see if convergence is obtained
        for(uint64_t i_x = 0; i_x<NBINS; i_x++){
            if(IY[i_x]<=threshold){
                return true;
            }
        }
        return false;
        
    }

    
};

/*int main(){*/

    /****************************************
     * IMPORT
     * *************************************/

    /*string file_in1 = "/home/aubin/ritchie/gecko/LargeScaleCLLGenome/05-ML/dev/GECKO_source/algoGen/testMic/subpop1_9_Dir/0_2/goodIndividuals.csv";
    string file_in2 = "/home/aubin/ritchie/gecko/LargeScaleCLLGenome/05-ML/dev/GECKO_source/algoGen/testMic/subpop2_9_Dir/0_2/goodIndividuals.csv";
    string file_in3 = "/home/aubin/ritchie/gecko/LargeScaleCLLGenome/05-ML/dev/GECKO_source/algoGen/testMic/subpop3_9_Dir/0_2/goodIndividuals.csv";*/

    //les 3 replicat de microRNA
/*    string file_in1 = "/home/aubin/ritchie/gecko/LargeScaleCLLGenome/05-ML/Logs_microRNA_rep1_Dir/0_15/goodIndividuals.csv";
    ClassGoodIndividuals inst1 = ClassGoodIndividuals(file_in1);
    inst1.importCSVfile("adjust",1); cerr << inst1._self_SequenceID.size() << " kmers for sample1" << endl;
    cout << inst1.getConvergence(20,4) << endl;*/

    /*string file_in2 = "/home/aubin/ritchie/gecko/LargeScaleCLLGenome/05-ML/Logs_microRNA_rep2_Dir/0_20/goodIndividuals.csv";
    string file_in3 = "/home/aubin/ritchie/gecko/LargeScaleCLLGenome/05-ML/Logs_microRNA_rep3_Dir/0_20/goodIndividuals.csv";
    
    ClassGoodIndividuals inst1 = ClassGoodIndividuals(file_in1);
    ClassGoodIndividuals inst2 = ClassGoodIndividuals(file_in2);
    ClassGoodIndividuals inst3 = ClassGoodIndividuals(file_in3);

    ClassGoodIndividuals finst1 = ClassGoodIndividuals(file_in1);
    ClassGoodIndividuals finst2 = ClassGoodIndividuals(file_in2);
    ClassGoodIndividuals finst3 = ClassGoodIndividuals(file_in3);

    ClassGoodIndividuals ffinst1 = ClassGoodIndividuals(file_in1);
    ClassGoodIndividuals ffinst2 = ClassGoodIndividuals(file_in2);
    ClassGoodIndividuals ffinst3 = ClassGoodIndividuals(file_in3);

    //importation standard
    inst1.importCSVfile("adjust",1); cerr << inst1._self_SequenceID.size() << " kmers for sample1" << endl;
    inst2.importCSVfile("adjust",1); cerr << inst2._self_SequenceID.size() << " kmers for sample2" << endl;
    inst3.importCSVfile("adjust",1); cerr << inst3._self_SequenceID.size() << " kmers for sample3" << endl;

    //importation standard avec filtrage  0.9
    finst1.importCSVfile("adjust",1,0.9); cerr << finst1._self_SequenceID.size() << " kmers for sample1" << endl;
    finst2.importCSVfile("adjust",1,0.9); cerr << finst2._self_SequenceID.size() << " kmers for sample2" << endl;
    finst3.importCSVfile("adjust",1,0.9); cerr << finst3._self_SequenceID.size() << " kmers for sample3" << endl;

    //importation standard avec filtrage à 0.9 + juste le meilleur
    ffinst1.importCSVfile("adjust",1,0.9, true); cerr << ffinst1._self_SequenceID.size() << " kmers for sample1" << endl;
    ffinst2.importCSVfile("adjust",1,0.9, true); cerr << ffinst2._self_SequenceID.size() << " kmers for sample2" << endl;
    ffinst3.importCSVfile("adjust",1,0.9, true); cerr << ffinst3._self_SequenceID.size() << " kmers for sample3" << endl;

    //inst1.computeNetwork("max", inst1._self_SequenceID);
    //inst2.computeNetwork("max", inst2._self_SequenceID);
    //inst3.computeNetwork("max", inst3._self_SequenceID);

    finst1.computeNetwork("max", finst1._self_SequenceID);
    finst2.computeNetwork("max", finst2._self_SequenceID);
    finst3.computeNetwork("max", finst3._self_SequenceID);

    ffinst1.computeNetwork("max", ffinst1._self_SequenceID);
    ffinst2.computeNetwork("max", ffinst2._self_SequenceID);
    ffinst3.computeNetwork("max", ffinst3._self_SequenceID);*/

    /****************************************
     * CONVERGENCE
     * *************************************/

    //output of convergence
    /*inst1.getNewKmerPerIteration("testMic/rep1_std.csv");
    inst2.getNewKmerPerIteration("testMic/rep2_std.csv");
    inst3.getNewKmerPerIteration("testMic/rep3_std.csv");
    
    finst1.getNewKmerPerIteration("testMic/rep1_filt1.csv");
    finst2.getNewKmerPerIteration("testMic/rep2_filt1.csv");
    finst3.getNewKmerPerIteration("testMic/rep3_filt1.csv");

    ffinst1.getNewKmerPerIteration("testMic/rep1_filt2.csv");
    ffinst2.getNewKmerPerIteration("testMic/rep2_filt2.csv");
    ffinst3.getNewKmerPerIteration("testMic/rep3_filt2.csv");*/
    


    /****************************************
     * OVERLAPS
     * *************************************/

    /*vector<uint64_t> line1 = inst1.overlapSequencesIDs(inst2._self_SequenceID);
    vector<uint64_t> line2 = inst1.overlapSequencesIDs(inst3._self_SequenceID);
    vector<uint64_t> line3 = inst2.overlapSequencesIDs(inst3._self_SequenceID);
    vector<uint64_t> line4 = inst1.overlapSequencesIDs(inst2.crossSequencesIDs(inst3._self_SequenceID,"intersect"));
    cerr << line1[0] << " , " << line1[1] << " , " << line1[2] << endl;
    cerr << line2[0] << " , " << line2[1] << " , " << line2[2] << endl;
    cerr << line3[0] << " , " << line3[1] << " , " << line3[2] << endl;
    cerr << line4[1] << endl;
    cerr << endl << endl ;

    vector<uint64_t> fline1 = finst1.overlapSequencesIDs(finst2._self_SequenceID);
    vector<uint64_t> fline2 = finst1.overlapSequencesIDs(finst3._self_SequenceID);
    vector<uint64_t> fline3 = finst2.overlapSequencesIDs(finst3._self_SequenceID);
    vector<uint64_t> fline4 = finst1.overlapSequencesIDs(finst2.crossSequencesIDs(finst3._self_SequenceID,"intersect"));
    cerr << fline1[0] << " , " << fline1[1] << " , " << fline1[2] << endl;
    cerr << fline2[0] << " , " << fline2[1] << " , " << fline2[2] << endl;
    cerr << fline3[0] << " , " << fline3[1] << " , " << fline3[2] << endl;
    cerr << fline4[1] << endl;
    cerr << endl << endl ;

    vector<uint64_t> ffline1 = ffinst1.overlapSequencesIDs(ffinst2._self_SequenceID);
    vector<uint64_t> ffline2 = ffinst1.overlapSequencesIDs(ffinst3._self_SequenceID);
    vector<uint64_t> ffline3 = ffinst2.overlapSequencesIDs(ffinst3._self_SequenceID);
    vector<uint64_t> ffline4 = ffinst1.overlapSequencesIDs(ffinst2.crossSequencesIDs(ffinst3._self_SequenceID,"intersect"));
    cerr << ffline1[0] << " , " << ffline1[1] << " , " << ffline1[2] << endl;
    cerr << ffline2[0] << " , " << ffline2[1] << " , " << ffline2[2] << endl;
    cerr << ffline3[0] << " , " << ffline3[1] << " , " << ffline3[2] << endl;
    cerr << ffline4[1] << endl;
    cerr << endl << endl ;*/


    /****************************************
     * HISTOGRAMS
     * *************************************/

    /*map<string, uint64_t> result1 = inst1.getHistogram();
    map<string, uint64_t> result2 = inst2.getHistogram();
    map<string, uint64_t> result3 = inst3.getHistogram();

    map<string, long double> result10 = finst1.getHistogramOfContacts();
    map<string, long double> result20 = finst2.getHistogramOfContacts();
    map<string, long double> result30 = finst3.getHistogramOfContacts();*/

    /*for(auto it = result1.begin(); it!=result1.end(); it++){
            cerr << it->first << "\t" << it->second << endl;
    }*/

    /*vector<string> allkmers = inst1.crossSequencesIDs(inst2.crossSequencesIDs(inst3._self_SequenceID,"merge"), "merge");
    cout << "ID\treplicate1\treplicate2\treplicate3\tcontact1\tcontact2\tcontact3" << endl;
    for(uint64_t k = 0; k<allkmers.size(); k++){
        uint64_t ck1=0;
        uint64_t ck2=0;
        uint64_t ck3=0;

        long double ck10=0;
        long double ck20=0;
        long double ck30=0;

        if(result1.find(allkmers[k])!=result1.end()){
            ck1 = result1[allkmers[k]];
        }
        if(result2.find(allkmers[k])!=result2.end()){
            ck2 = result2[allkmers[k]];
        }
        if(result3.find(allkmers[k])!=result3.end()){
            ck3 = result3[allkmers[k]];
        }

        if(result10.find(allkmers[k])!=result10.end()){
            ck10 = result10[allkmers[k]];
        }
        if(result20.find(allkmers[k])!=result20.end()){
            ck20 = result20[allkmers[k]];
        }
        if(result30.find(allkmers[k])!=result30.end()){
            ck30 = result30[allkmers[k]];
        }


        cout << allkmers[k] << "\t" << ck1 << "\t" << ck2 << "\t" << ck3 << "\t" << ck10 << "\t" << ck20 << "\t" << ck30 << endl;
    }*/


    /****************************************
     * MATRIX
     * *************************************/
    
    
    /*inst1.outputMatrix("inst1_std_max.matrix", inst1._self_SequenceID);
    inst2.outputMatrix("inst2_std_max.matrix", inst2._self_SequenceID);
    inst3.outputMatrix("inst3_std_max.matrix", inst3._self_SequenceID);*/

//    finst1.outputMatrix("inst1_filt1_max.matrix", finst1._self_SequenceID);
//    finst2.outputMatrix("inst2_filt1_max.matrix", finst2._self_SequenceID);
//    finst3.outputMatrix("inst3_filt1_max.matrix", finst3._self_SequenceID);

//    ffinst1.outputMatrix("inst1_filt2_max.matrix", ffinst1._self_SequenceID);
//    ffinst2.outputMatrix("inst2_filt2_max.matrix", ffinst2._self_SequenceID);
//    ffinst3.outputMatrix("inst3_filt2_max.matrix", ffinst3._self_SequenceID);
    



    /*test.computeNetwork("max", intersectKmers);
    test2.computeNetwork("max", intersectKmers);
    test3.computeNetwork("max", intersectKmers);
    test4.computeNetwork("max", intersectKmers);*/
    
    /*test.outputMatrix("sample1.tsv", intersectKmers);
    test2.outputMatrix("sample3.tsv", intersectKmers);
    test3.outputMatrix("sample4.tsv", intersectKmers);
    test4.outputMatrix("sample5.tsv", intersectKmers);*/
    
    
    //void computeNetwork(string p_method, vector<string> indiv_seq, vector<float> indiv_scores, map<string, map<string, bool> indiv_seq_state, map<string, uint64_t> kmer_position){
    /*test.computeNetwork("max", test._self_SequenceID);
    test2.computeNetwork("average", test2._self_SequenceID);
    test3.computeNetwork("sum", test3._self_SequenceID);
    test4.computeNetwork("singleAverage", test4._self_SequenceID);
    test5.computeNetwork("singleSum", test5._self_SequenceID);
        

    test.outputMatrix("methodmax.tsv");
    test2.outputMatrix("methodaverage.tsv");
    test3.outputMatrix("methodsum.tsv");
    test4.outputMatrix("methodsingleaverage.tsv");
    test5.outputMatrix("methodsinglesum.tsv");*/
/*    return EXIT_SUCCESS;
}*/
