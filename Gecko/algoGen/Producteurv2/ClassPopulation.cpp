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
#include "ClassMatrixAccess.cpp"
#include "ExtractNetwork.cpp"


using namespace std;

/*!
 * \file ClassPopulation.cpp
 * \brief class to instantiate population model for genetic algorithm.
 * \class Class ClassPopulation
 * The genetic algorithm is applied on a population of individuals that are made of k-mers.
 * Each individual is composed with fix number of k-mers.
 * The k-mers are imported from a matrix file in csv format. columns are k-mers and lines are samples. The values are k-mer enrichments,
 * except for the 1st column that refers to the group of the samples.
 * The role of the class ClassPopulation is to import the matrix, create population, divide for this population the samples.
 * The samples are divided into 3 parts : training, test and outer (outter, sorry).
 *
 * The class can store for each k-mer the number of time it appears through the population.
 * The class can store the resulys of several methods at a time
 *
*/

struct RankAndValue {
  uint64_t rank;
  float value;
} ;


bool compareByValue(const RankAndValue &a, const RankAndValue &b){
    return a.value < b.value;
}
bool compareByValueInv(const RankAndValue &a, const RankAndValue &b){
    return a.value > b.value;
}


class ClassPopulation {
public:
    /// \param initial_sequenceID the vector of sequence ids in the matrix (order is conserved)
    /// \param initial_sequenceCounts the vector of counts of k-mers (number of times it appears)
    /// \param initial_sequenceCountsScore the vector of cumulative test scores for each k-mer

    /// \param initial_data_x the big matrix (l)x(c) = (sample)x(kmer)
    /// \param initial_data_y the vector of groups in string ex : Yes No NA Yes No NA No Yes NA
    /// \param initial_data_group data_y put in integer values ex : 1 2 3 1 2 3 2 1 3
    /// \param initial_data_replicate the cumulative vector of apparition of the group through samples ex : 1 1 1 2 2 2 3 3 3
    /// \param initial_data_state : vector that indicates if the samples is considered as training, test or outer (respectively 0, SUM(2^crossvalidation), -1)

    /// \param dictocc hash group->occurence ex : Yes->3
    /// \param dim_groups number of groups
    /// \param dim_Kmers number of k-mers
    /// \param _self_sizeKmer size of kmers in nucleotides
    /// \param dim_outter number of outers

    /// \param table 2D that refers to list of individuals
    /// \param tableWinners table of best individuals at each iteration, key related to method
    vector <string> initial_sequenceID; // the kmer ids
    vector <uint64_t> initial_sequenceCounts; //the counts of occurences for each kmer
    vector <long double> initial_sequenceCountsScore;//the sum of testscores for each occurences of each kmer

    vector <vector<long double>> initial_data_x; // the big matrix (l)x(c) = (sample)x(kmer)
    vector <string> initial_data_y; // groups in string
    vector <string> initial_data_indiv_name; // groups in string
    vector <uint64_t> initial_data_group; // groups in uint >0
    vector <uint64_t> initial_data_replicate; // replicates of groups
    vector <int> initial_data_state; //-1 : outter ; 1 : test ; 0 : train
    map< string, uint64_t> dictocc; //occurences of sample per group

    uint64_t dim_groups;
    uint64_t dim_Kmers;
    uint64_t dim_outter;
    uint64_t _self_sizeKmer;

    vector <vector<uint64_t>> table; //table of individuals composed of Kmers
    map<uint64_t, vector <vector<uint64_t>>> tableWinners; //table of best individuals at each iteration, key related to method

    map<uint64_t, vector <vector<uint64_t>>> tableWinnersWithOuters_individuals; //table of best individuals at each iteration, key related to method
    map<uint64_t, vector <uint64_t>> tableWinnersWithOuters_iteration; //iteration respectively to tableWinnersWithOuters_individuals , key related to method
    map<uint64_t, vector <long double>> tableWinnersWithOuters_test; //test score respectively to tableWinnersWithOuters_individuals , key related to method
    map<uint64_t, vector <long double>> tableWinnersWithOuters_outer; //outer score respectively to tableWinnersWithOuters_individuals , key related to method

    ClassMatrixAccess initial_binary_file; //object that will handle binary file if the matrix is proposed as binary format

    /*!
     *  \brief constructor
     *
     *  take as input the number of methods of scores (default 1)
     */
    ClassPopulation(uint64_t pNumOfMethods = 1) {
        for (uint64_t i = 0; i < pNumOfMethods; i++) {
            vector <vector < uint64_t>> tmp;
            vector < uint64_t> tmp2;
            vector < long double> tmp3;
            tableWinners[i] = tmp;
            tableWinnersWithOuters_individuals[i] = tmp;
            tableWinnersWithOuters_iteration[i] = tmp2;
            tableWinnersWithOuters_test[i] = tmp3;
            tableWinnersWithOuters_outer[i] = tmp3;
        }
    }


    /*!
    * \brief return true if convergence is obtained
    * the method use the implementation of the class ClassGoodIndividuals:
    * */
    bool convergence(string path_goodCandidates, string path_outputfile){
        ClassGoodIndividuals instGC = ClassGoodIndividuals(path_goodCandidates, _self_sizeKmer);
        return instGC.getConvergence(path_outputfile);
    }


    /*!
     *  \brief return the number of sample analysis except outers
     */
    uint64_t getDimSampleAnalysis() {
        return initial_data_y.size() - 2 * dim_outter;
    }

    /*!
     *  \brief return if the candidate as individual is already in the table
     *   To avoid having twice the same individual (rare event)
     */
    bool checkifunique(vector<uint64_t> src,uint64_t NKmer){
        //vector<uint64_t> srcsort=std::sort(table[isrc]);
        bool bEqual=true;
        for(uint64_t i=0;i<table.size();i++){
            bEqual = boost::algorithm::is_permutation(src.begin(),src.begin()+NKmer,table[i].begin());

            if (bEqual==true){
                break;
            }

        }
        return !bEqual;
    }


    /*!
    * \ brief string to int for groups + compute group numbers
    * transform the groups in string as group in int value
    * and counts the occurences of each group
    */
    void computeGroups(bool debug) {
        map< string, uint64_t> dict;
        dim_groups = 0;
        for (uint64_t i = 0; i < initial_data_y.size(); i++) {
            if (dict.find(initial_data_y[i]) == dict.end()) {
                dict[initial_data_y[i]] = dim_groups + 1;
                dictocc[initial_data_y[i]] = 0;
                if (debug) {
                    cerr << "group = " << initial_data_y[i] << ",id = " << dim_groups + 1 << endl;
                }
                dim_groups++;
            }
            dictocc[initial_data_y[i]] += 1;
            initial_data_group.push_back(dict[initial_data_y[i]]);
            initial_data_replicate.push_back(dictocc[initial_data_y[i]]);
            initial_data_state.push_back(0);
        }
        if (debug) {
            cerr << "dim_groups = " << dim_groups << endl;
            cerr << "group\treplicate\tstate\tdictocc" << endl;
            for (uint64_t d = 0; d < initial_data_y.size(); d++) {
                cerr << initial_data_group[d] << "\t" <<
                        initial_data_replicate[d] << "\t" <<
                        initial_data_state[d] << "\t" <<
                        dictocc[initial_data_y[d]] << endl;
            }
        }

    }

    //! \brief transform sequence ID to real nucleotidic sequence
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


    /*!
     *  \brief return the size of the kmers
     * if all the kmers of the sample begins with A, then decrease
     */
    void getSizeKmers(float percentage_dict=0.3, uint64_t size_init = 31, uint64_t result_percentage_limit = 5000){
        srand (time(NULL));
        uint64_t numberOfKmers = floor(percentage_dict*dim_Kmers);
        if(result_percentage_limit<numberOfKmers){
            numberOfKmers = result_percentage_limit;
        }
        bool state = true;
        while(state){
            uint64_t nAstart = 0;
            for(uint64_t k=0; k<numberOfKmers; k++){
                string tmp = initial_binary_file.idseq2fasta(strtoull(initial_sequenceID[rand() % dim_Kmers].c_str(),NULL,0), size_init);
                if(tmp[0]=='A'){
                    nAstart++;
                }
            }
            if(nAstart==numberOfKmers){
                size_init--;
            }
            else{
                _self_sizeKmer = size_init;
                state = false;
                cerr << "Kmer size detected is " << _self_sizeKmer << endl;
            }
        }
    }


    /*!
    * \brief get the number of samples for the most little group
    *
    * */
    uint64_t getMinGroup(bool debug = false) {
        uint64_t result = initial_data_group.size();
        for (map< string, uint64_t>::iterator it = dictocc.begin(); it != dictocc.end(); it++) {
            if (debug) {
                cerr << it->first << "\t" << it->second << endl;
            }
            if (it->second < result) {
                result = it->second;
            }
        }
        return result;
    }


    /*!
    * \brief the loadCSV method is used to integrate a matrix. This matrix can be in text or binary format
    *
    * */
    void loadCSV(string param_pathToFile, bool debug = false) {
        ClassMatrixAccess tmp     = ClassMatrixAccess(param_pathToFile,100000);
        initial_binary_file       = tmp;

        if(initial_binary_file._self_isInBinary){
            loadCSV_binary(param_pathToFile, debug);
        }
        else{
            loadCSV_text(param_pathToFile, debug);
            //loadCSV_text_sauv(param_pathToFile, debug);
        }
        getSizeKmers();//autodetect the size of the kmers
        loadIndivName(param_pathToFile);
    }
    /*!
    * \brief the loadCSV method is used to load the individuals name from a matrix. This matrix can be in text or binary format
    *
    * */
    void loadIndivName(string param_pathToFile) {

        if(initial_binary_file._self_isInBinary){
            initial_data_indiv_name=initial_binary_file.getNames();
            initial_data_indiv_name.erase(initial_data_indiv_name.begin());

        }
        else{
            ifstream file_ai(param_pathToFile.c_str(), ios::in);
            if (!file_ai) {
                cerr << "Error while opening " << param_pathToFile << " in read mode" << endl;
                exit(EXIT_FAILURE);
            }

            string header;
            vector<string> contentheader;
            getline(file_ai, header);
            boost::split(contentheader, header, boost::is_any_of("\t"));
            for (uint64_t i = 1; i < contentheader.size(); i++) {
                initial_data_indiv_name.push_back(contentheader[i]);
            }
            file_ai.close();

        }
    }



    /*!
    * \brief load the csv file and computes the groups
    */
    void loadCSV_text(string param_pathToFile, bool debug = false) {
        //boost::split(content, line, boost::is_any_of("\t"));
        ifstream file_ai(param_pathToFile.c_str(), ios::in);
        if (!file_ai) {
            cerr << "Error while opening " << param_pathToFile << " in read mode" << endl;
            exit(EXIT_FAILURE);
        }

        string header, line;
        vector<string> contentheader;
        getline(file_ai, header);
        getline(file_ai, header);
        boost::split(contentheader, header, boost::is_any_of("\t"));
        for (uint64_t i = 1; i < contentheader.size(); i++) {
            initial_data_y.push_back(contentheader[i]);
            vector<long double> ldvtmp;
            initial_data_x.push_back(ldvtmp);
        }
        while (getline(file_ai, line)) {
            vector<string> content;
            boost::split(content, line, boost::is_any_of("\t"));
            initial_sequenceID.push_back(content[0]);
            initial_sequenceCounts.push_back(0);
            initial_sequenceCountsScore.push_back(0);
            for (uint64_t i = 1; i < content.size(); i++)
                initial_data_x[i-1].push_back(strtold(content[i].c_str(), NULL));

        }

        file_ai.close();
        //computing groups
        computeGroups(debug);
        dim_Kmers = initial_data_x[0].size();
        dim_outter = 0;
    }
    void loadCSV_text_sauv(string param_pathToFile, bool debug = false) {
        //boost::split(content, line, boost::is_any_of("\t"));
        ifstream file_ai(param_pathToFile.c_str(), ios::in);
        if (!file_ai) {
            cerr << "Error while opening " << param_pathToFile << " in read mode" << endl;
            exit(EXIT_FAILURE);
        }

        string header, line;
        vector<string> contentheader;
        getline(file_ai, header);
        boost::split(contentheader, header, boost::is_any_of(","));
        for (uint64_t i = 1; i < contentheader.size(); i++) {
            initial_sequenceID.push_back(contentheader[i]);
            initial_sequenceCounts.push_back(0);
            initial_sequenceCountsScore.push_back(0);
        }

        dim_Kmers = initial_sequenceID.size();
        if (debug) {
            cerr << "dim_Kmers = " << dim_Kmers << endl;
        }
        while (getline(file_ai, line)) {
            vector<string> content;
            boost::split(content, line, boost::is_any_of(","));
            initial_data_y.push_back(content[0]);
            vector<long double> vtmp(content.size() - 1, 0.0);
            for (uint64_t i = 1; i < content.size(); i++)
                vtmp[i - 1] = strtold(content[i].c_str(), NULL);
            initial_data_x.push_back(vtmp);

        }
        file_ai.close();
        //computing groups
        computeGroups(debug);
        dim_outter = 0;

    }

    /*!
    * \brief load the binary file and computes the groups
    */
    void loadCSV_binary(string param_pathToFile, bool debug = false) {

        initial_sequenceID = initial_binary_file.getSequenceIDs();
        dim_Kmers = initial_sequenceID.size();
        //dim_Kmers = initial_binary_file._self_nbLines-2; //2 columns of headers

        vector<uint64_t> tmpi (dim_Kmers,0);
        vector<long double> tmpld (dim_Kmers,0.0);
        initial_sequenceCounts = tmpi;
        initial_sequenceCountsScore = tmpld;

        if (debug) {
            cerr << "dim_Kmers = " << dim_Kmers << endl;
        }

        //computing groups
        initial_data_y = initial_binary_file.getGroups();
        computeGroups(debug);

        dim_outter = 0;
    }

    /*!
    * \brief get the param_maxNumberOfKmers most variant kmers and output into a text matrix.
    * The CSV must have been loaded before the launch
    * for each kmer the ratio variance/mean is computed
    * then they are ranked, divided into 20 bins.
    * then a z-score is applied inside each bin, and a z-score-cutoff is chosen to find the param_maxNumberOfKmers most variant genes
    * based on "Spatial reconstruction of single-cell gene expression data, Rahul Satija, Jeffrey A Farrell, David Gennert, Alexander F Schier & Aviv Regev"
    */
    void MakeMatrixFromMostVariantKmers(string param_pathToOutputMatrixMostVariant, string param_pathToOutputMatrixLessVariant, uint64_t param_maxNumberOfKmers, bool debug = false){
        RankAndValue tmp;
        tmp.rank = 0;
        tmp.value = 0;
        vector<RankAndValue> ValueKmers(dim_Kmers,tmp);

        //compute ratio variance/mean for every kmers
        if(debug){
            cerr << "Computing ratios with dim_Kmers = " << dim_Kmers << endl;
        }
        #pragma omp parallel for
        for(uint64_t k = 0; k<dim_Kmers; k++){
            float mean = 0.0;
            float variance = 0.0;
            if(debug){
                cerr << "k = "<< k << endl;
            }

            //get the value of the kmer
            uint64_t n_samples = initial_data_y.size();
            vector<float> result;
            if(initial_binary_file._self_isInBinary){
                result = initial_binary_file.getValuesFromSingleKmer(k, false);
            }
            else{
                for (uint64_t i = 0; i < n_samples; i++) {
                    result.push_back(initial_data_x[i][k]);
                }
            }
            if(debug){
                cerr << "\tinitial values : "  << result[0];
                for (uint64_t i = 1; i < n_samples; i++){
                    cerr << "," << result[i];
                }
                cerr << endl;
            }

            //compute mean
            for (uint64_t i = 0; i < n_samples; i++){
                mean+=result[i];
            }
            mean/=n_samples;
            if(debug){
                cerr << "\tmean  = " << mean << endl;
            }

            //compute variance
            for (uint64_t i = 0; i < n_samples; i++){
                variance+=((result[i]-mean)*(result[i]-mean));
            }
            variance/=n_samples;
            if(debug){
                cerr << "\tvariance  = " << variance << endl;
            }

            ValueKmers[k].rank = k;
            ValueKmers[k].value = mean/variance;
            if(debug){
                cerr << "\tratio  = " << ValueKmers[k].value << endl;
            }
        }

        //sort the values
        sort(ValueKmers.begin(), ValueKmers.end(), compareByValue);
        if(debug){
            cerr << endl << endl;
            for(uint64_t idebug = 0; idebug<dim_Kmers; idebug++){
                cerr <<  ValueKmers[idebug].rank << "\t" << ValueKmers[idebug].value << endl;
            }
        }


        //cut ValueKmers into 20 bins, and computes the z-score of each kmer inside each bin
        uint64_t sizeBin = floor(dim_Kmers/20);
        #pragma omp parallel for
        for(uint64_t bc = 0; bc<20; bc++){
            float mean = 0.0;
            float variance = 0.0;
            uint64_t kstart = sizeBin*bc;
            uint64_t kstop = sizeBin*(bc+1);
            uint64_t ksize = sizeBin;
            if(bc==19){
                kstop = dim_Kmers;
                ksize = kstop-kstart;
            }

            for(uint64_t k = kstart; k<kstop; k++){
                mean+=ValueKmers[k].value;
            }
            mean/=ksize;

            for(uint64_t k = kstart; k<kstop; k++){
                variance+= ((ValueKmers[k].value-mean)*(ValueKmers[k].value-mean));
            }
            variance=sqrt(variance/(ksize-1));
            //compute z-score
            for(uint64_t k = kstart; k<kstop; k++){
                ValueKmers[k].value = (ValueKmers[k].value-mean)/variance;
            }
        }

        //sort the value
        sort(ValueKmers.begin(), ValueKmers.end(), compareByValueInv);


        //create the output matrix
        vector<string> sample_names;
        if(initial_binary_file._self_isInBinary){
            sample_names = initial_binary_file.getNames();
        }
        else{

        }
        ofstream file_ao(param_pathToOutputMatrixMostVariant.c_str(), ios::out);
        ofstream file_ao2(param_pathToOutputMatrixLessVariant.c_str(), ios::out);

        for(uint64_t iy=0; iy<initial_data_y.size(); iy++){
            file_ao << "," << initial_data_indiv_name[iy];
            file_ao2 << "," << initial_data_indiv_name[iy];
        }
        file_ao << endl;
        file_ao2 << endl;
        for(uint64_t iy=0; iy<initial_data_y.size(); iy++){
            file_ao << "," << initial_data_y[iy];
            file_ao2 << "," << initial_data_y[iy];
        }
        file_ao << endl;
        file_ao2 << endl;

        for(uint64_t i = 0; i<param_maxNumberOfKmers ; i++){
            file_ao << idseq2fasta(strtoull(initial_sequenceID[ValueKmers[i].rank].c_str(),NULL,0), _self_sizeKmer);
            for(uint64_t iy=0; iy<initial_data_y.size(); iy++){
                if(initial_binary_file._self_isInBinary){
                    file_ao << "," << initial_binary_file.getValue(iy, ValueKmers[i].rank, "Population::MakeMatrixFromMostVariantKmers");
                }
                else{
                    file_ao << "," << initial_data_x[iy][ValueKmers[i].rank];
                }
                
            }
            file_ao << endl;
        }

        for(uint64_t i = ValueKmers.size()-1; i>ValueKmers.size()-1-param_maxNumberOfKmers ; i--){
            file_ao2 << idseq2fasta(strtoull(initial_sequenceID[ValueKmers[i].rank].c_str(),NULL,0), _self_sizeKmer);
            for(uint64_t iy=0; iy<initial_data_y.size(); iy++){
                if(initial_binary_file._self_isInBinary){
                    file_ao2 << "," << initial_binary_file.getValue(iy, ValueKmers[i].rank, "Population::MakeMatrixFromMostVariantKmers");
                }
                else{
                    file_ao2 << "," << initial_data_x[iy][ValueKmers[i].rank];
                }
                
            }
            file_ao2 << endl;
        }

        file_ao.close();
        file_ao2.close();

    }

    /*void loadKmerValues(){
        vector<uint64_t> kmerListGeneration;

        if(initial_binary_file._self_isInBinary){

            for(uint64_t indiv=0; indiv<table.size(); indiv++){
                for(uint64_t kmer=0; kmer<table[indiv].size(); kmer++){
                    if (find(kmerListGeneration.begin(),kmerListGeneration.end(),table[indiv][kmer])==kmerListGeneration.end()){
                        kmerListGeneration.push_back(table[indiv][kmer]);
                        //cout<<kmerListGeneration.size()<<"   kmer:"<<table[indiv][kmer]<<endl;
                    }


                }
            }
          //cout<<"number of individual in loadKmerValues = "<<table.size()<<endl;
            initial_binary_file.loadValue(kmerListGeneration,table.size(),table[0].size());
        }


    }*/


    /*void loadKmerWinnersValues(uint64_t kellogs,uint64_t begini,uint64_t endi){
        vector<uint64_t> kmerListGeneration;

        if(initial_binary_file._self_isInBinary){

            for (uint64_t guy = begini; guy < endi; guy++){
                for(uint64_t kmer=0; kmer<tableWinners[kellogs][guy].size(); kmer++){
                    if (find(kmerListGeneration.begin(),kmerListGeneration.end(),tableWinners[kellogs][guy][kmer])==kmerListGeneration.end()){
                        kmerListGeneration.push_back(tableWinners[kellogs][guy][kmer]);
                        //cout<<kmerListGeneration.size()<<"   kmer:"<<table[indiv][kmer]<<endl;
                    }


                }
            }
            cout<<"number of individual in loadKmerWinnersValues = "<<endi-begini<<endl;
            initial_binary_file.loadValue(kmerListGeneration,endi-begini,tableWinners[kellogs][begini].size());
        }


    }*/

    /*!
    * \brief print a submatrix to output file
    * the output file is param_pathToFile
    * the input file is param_path2occs
    * threshold for the minimum kmer count to keep the kmer
    *
    * NOT BINARY CAPABLE
    */
    void printDatas(string param_path2occs,string param_pathToFile,u_int64_t threshold) {

        ofstream file_ao(param_pathToFile.c_str(), ios::out);
        ifstream file_afi(param_path2occs.c_str(), ios::in);
        u_int64_t countkm;

        if (!file_afi) {
            cerr << "Error while opening " << param_path2occs << " in read mode" << endl;
            exit(EXIT_FAILURE);
        }
        //map<string, long double> catalog;
        string line;
        vector<string> sortkmerid;
        cerr << "reading " << param_path2occs.c_str() << "...";
        while (getline(file_afi, line)) {
            if(!line.empty()){
                //cout<<"line "<<line<<"fin"<<endl;
                vector<string> content;
                boost::split(content, line, boost::is_any_of(","));
                boost::trim(content[0]);
                boost::trim(content[1]);
                countkm=strtold(content[1].c_str(),NULL);
                if (countkm>=threshold)
                    sortkmerid.push_back(content[0].c_str());
            }
            //boost::trim(content[1]);
            //cout<<content[0]<<","<<content[1]<<endl;

        }
        cerr << "Done" << endl;
        file_afi.close();


        //is the kmer inside the matrix
        vector <uint64_t> resortid;
        for(uint64_t sortkmer_indice = 0; sortkmer_indice < sortkmerid.size(); sortkmer_indice++) {
            for (uint64_t kmer_indice = 0; kmer_indice < dim_Kmers; kmer_indice++) {
                if (initial_sequenceID[kmer_indice]==sortkmerid[sortkmer_indice]){
                    resortid.push_back(kmer_indice);
                    break;
                }
            }
        }


        if (file_ao) {
            for (uint64_t grp_indice = 0; grp_indice < initial_data_y.size(); grp_indice++) {
               file_ao <<"\t"<<initial_data_indiv_name[grp_indice];
            }
            file_ao << endl;
            for (uint64_t grp_indice = 0; grp_indice < initial_data_y.size(); grp_indice++) {
               file_ao <<"\t"<<initial_data_y[grp_indice];
            }
            file_ao << endl;
            for (uint64_t kmer_indice = 0; kmer_indice < resortid.size(); kmer_indice++) {
                file_ao <<initial_sequenceID[resortid[kmer_indice]];
                for (uint64_t grp_indice = 0; grp_indice < initial_data_y.size(); grp_indice++) {

                    if(initial_binary_file._self_isInBinary){
                        file_ao <<"\t"<<initial_binary_file.getValue(grp_indice,resortid[kmer_indice],"printDatas");
                    }
                    else{
                        file_ao <<"\t"<<initial_data_x[grp_indice][resortid[kmer_indice]];
                    }
                 }
                 file_ao << endl;
            }

            file_ao.close();
        }
    }
     void printDatas_old_beforeDetranslate(string param_path2occs,string param_pathToFile,u_int64_t threshold) {

        ofstream file_ao(param_pathToFile.c_str(), ios::out);
        ifstream file_afi(param_path2occs.c_str(), ios::in);
        u_int64_t countkm;

        if (!file_afi) {
            cerr << "Error while opening " << param_path2occs << " in read mode" << endl;
            exit(EXIT_FAILURE);
        }
        //map<string, long double> catalog;
        string line;
        vector<string> sortkmerid;
        cerr << "reading " << param_path2occs.c_str() << "...";
        while (getline(file_afi, line)) {
            if(!line.empty()){
                //cout<<"line "<<line<<"fin"<<endl;
                vector<string> content;
                boost::split(content, line, boost::is_any_of(","));
                boost::trim(content[0]);
                boost::trim(content[1]);
                countkm=strtold(content[1].c_str(),NULL);
                if (countkm>=threshold)
                    sortkmerid.push_back(content[0].c_str());
            }
            //boost::trim(content[1]);
            //cout<<content[0]<<","<<content[1]<<endl;

        }
        cerr << "Done" << endl;
        file_afi.close();


        //is the kmer inside the matrix
        vector <uint64_t> resortid;
        for(uint64_t sortkmer_indice = 0; sortkmer_indice < sortkmerid.size(); sortkmer_indice++) {
            for (uint64_t kmer_indice = 0; kmer_indice < dim_Kmers; kmer_indice++) {
                if (initial_sequenceID[kmer_indice]==sortkmerid[sortkmer_indice]){
                    resortid.push_back(kmer_indice);
                    break;
                }
            }
        }


        if (file_ao) {
            file_ao << "groups" ;
            for (uint64_t kmer_indice = 0; kmer_indice < resortid.size(); kmer_indice++) {
                file_ao << ","<< initial_sequenceID[resortid[kmer_indice]]  ;

            }
            file_ao <<endl;
            for (uint64_t grp_indice = 0; grp_indice < initial_data_y.size(); grp_indice++) {

                file_ao <<initial_data_y[grp_indice];

                for (uint64_t kmer_indice = 0; kmer_indice < resortid.size(); kmer_indice++) {
                    if(initial_binary_file._self_isInBinary){
                        file_ao <<","<<initial_binary_file.getValue(grp_indice,resortid[kmer_indice],"printDatas");
                    }
                    else{
                        file_ao <<","<<initial_data_x[grp_indice][resortid[kmer_indice]];
                    }
                }
                file_ao <<endl;
            }
            file_ao.close();
        }
    }
    /*!
    * \brief load the csv file with constraints and computes the groups
    *
    * path2occs is a csv file that contains k-mer and its scores
    * then import the raw matrix in csv
    *
    * NOT BINARY CAPABLE
    * */
    void loadCSVwithConstraint(string param_pathToFile, string param_path2occs, long double threshold, bool debug = false) {
        string header, line;
        ifstream file_afi(param_path2occs.c_str(), ios::in);
        ifstream file_ai(param_pathToFile.c_str(), ios::in);
        if (!file_ai) {
            cerr << "Error while opening " << param_pathToFile << " in read mode" << endl;
            exit(EXIT_FAILURE);
        }
        if (!file_afi) {
            cerr << "Error while opening " << param_path2occs << " in read mode" << endl;
            exit(EXIT_FAILURE);
        }

        map<string, long double> catalog;
        while (getline(file_afi, line)) {
            vector<string> content;
            boost::split(content, line, boost::is_any_of(","));
            boost::trim(content[0]);
            boost::trim(content[1]);
            //cout<<content[0]<<","<<content[1]<<endl;
            if (catalog.find(content[0]) == catalog.end())
                catalog[content[0]] = 0;
            catalog[content[0]] += strtold(content[1].c_str(), NULL);
        }
        file_afi.close();

        vector<string> toremove;
        if (threshold!=0){
            for (auto it = catalog.begin(); it != catalog.end(); it++) {
                if (it-> second < threshold) {
                    toremove.push_back(it->first);
                }
            }
            for (uint64_t i = 0; i < toremove.size(); i++) {
                catalog.erase(toremove[i]);
            }
        }

        //now importing file
        vector<string> contentheader;
        vector<bool> contentheaderIndic;
        getline(file_ai, header);
        boost::split(contentheader, header, boost::is_any_of(","));

        for (uint64_t i = 1; i < contentheader.size(); i++) {
            if (catalog.find(contentheader[i]) != catalog.end()) {

                contentheaderIndic.push_back(true);
                initial_sequenceID.push_back(contentheader[i]);
                initial_sequenceCounts.push_back(0);
                initial_sequenceCountsScore.push_back(0);
            } else {
                contentheaderIndic.push_back(false);
            }
        }

        dim_Kmers = initial_sequenceID.size();
        if (debug) {
            cerr << "dim_Kmers = " << dim_Kmers << endl;
        }
        while (getline(file_ai, line)) {
            vector<string> content;
            boost::split(content, line, boost::is_any_of(","));
            initial_data_y.push_back(content[0]);

            vector<long double> vtmp;
            for (uint64_t i = 1; i < content.size(); i++) {
                if (contentheaderIndic[i-1]==true){
                    vtmp.push_back(strtold(content[i].c_str(), NULL));
                }
            }
            initial_data_x.push_back(vtmp);

        }
        file_ai.close();
        //computing groups
        computeGroups(debug);
        dim_outter = 0;

    }

    /*!
    * \brief create the table of individuals from existing pre-made composition
    *
    * Each individual is made of param_NKmer k-mers
    * the file is a csv file, each line is an individual with param_NKmer (comma separated)
    *
    * */
    void readTable(uint64_t param_NKmer = 50, uint64_t paramNIndividuals = 10, string param_pathToFile = "bestindiv.csv", bool debug = false) {
        ifstream file_afi(param_pathToFile.c_str(), ios::in);
        string line;
        if (!file_afi) {
            cerr << "Error while opening " << param_pathToFile << " in read mode" << endl;
            exit(EXIT_FAILURE);
        }
        for (uint64_t ipop = 0; ipop < paramNIndividuals; ipop++) {
            vector<uint64_t> popk;
            vector<string> content;
            if (!getline(file_afi, line)) {
                if (!file_afi) {
                    cerr << "Error while opening " << param_pathToFile << " in read mode : not enough inviduals in the file : " << ipop << "/" << paramNIndividuals << endl;
                    exit(EXIT_FAILURE);
                }
            }

            boost::split(content, line, boost::is_any_of(","));
            for (unsigned int i = 0; i < content.size(); i++) {
                popk.push_back(strtold(content[i].c_str(), NULL));
            }
            if (debug) {
                for (uint64_t d = 0; d < popk.size(); d++) {
                    cerr << "\t" << popk[d];
                }
                cerr << endl;
            }
            table.push_back(popk);

        }
        file_afi.close();

    }

    /*!
    * \brief create the table of individuals
    *
    * Each individual is made of param_NKmer k-mers
    * there are paramNIndividuals in the table
    *
    * */
    void createTable(uint64_t param_NKmer = 50, uint64_t paramNIndividuals = 10, bool debug = false) {
    	if(debug){
    		cerr << "createTable( " << param_NKmer << " , " << paramNIndividuals << " )" << endl;
    	}
        for (uint64_t ipop = 0; ipop < paramNIndividuals; ipop++) {
            if(debug)
            	cerr << "pop " << ipop << endl;
            vector<uint64_t> popk;
            while (popk.size() < param_NKmer) {
                int pos_rand = rand() % dim_Kmers;
                bool state = false;
                //do not want twice the same kmer
                for (uint64_t isec = 0; isec < popk.size(); isec++) {
                    if (popk[isec] == (uint64_t) pos_rand) {
                        state = true;
                        isec = popk.size();
                    }
                } //for isec
                if (!state) {
                    popk.push_back((uint64_t) pos_rand);
                }
            } //while popk.size()
            if (debug) {
                for (uint64_t d = 0; d < popk.size(); d++) {
                    cerr << "\t" << popk[d];
                }
                cerr << endl;
            }
            table.push_back(popk);
        } //for ipop

        if(initial_binary_file._self_isInBinary){
            initial_binary_file.getValuesFromTable(table); //chargement tampon
        }
    }

    /*!
    * \brief update the counts of the k-mers and their scores
    * the scores are given in parameters
    *
    * */
    void updateCounts(vector <long double> resRandTree) {
        for (uint64_t i = 0; i < table.size(); i++) {
            for (uint64_t j = 0; j < table[i].size(); j++) {
                initial_sequenceCounts[table[i][j]] += 1;
                initial_sequenceCountsScore[table[i][j]] += resRandTree[i];

            }
        }
    }

    /*!
    * \brief update the individuals that have a sufficient score in outer and test
    * the scores are given in parameters
    *
    * */
    void updateOuterTest(uint64_t kellogs, vector <long double> resRandTree, vector <long double> resRandTreeOuter, long double testmin, long double outermin, uint64_t iteration) {
        for (uint64_t i = 0; i < table.size(); i++) {
            if(resRandTree[i]>=testmin and resRandTreeOuter[i]>=outermin){
                vector<uint64_t> tmp(table[i]);
                tableWinnersWithOuters_individuals[kellogs].push_back(tmp);
                tableWinnersWithOuters_iteration[kellogs].push_back(iteration);
                tableWinnersWithOuters_test[kellogs].push_back(resRandTree[i]);
                tableWinnersWithOuters_outer[kellogs].push_back(resRandTreeOuter[i]);
            }
        }
        cerr << endl << "ClassPopulation::updateOuterTest : " << tableWinnersWithOuters_outer[kellogs].size() << " good individuals" << endl;
    }

    /**
     * \brief clear the content of containers of good individuals
     * */
    void clearGoodIndividuals(uint64_t kellogs){
        tableWinnersWithOuters_individuals[kellogs].erase(tableWinnersWithOuters_individuals[kellogs].begin(), tableWinnersWithOuters_individuals[kellogs].end());
        tableWinnersWithOuters_iteration[kellogs].erase(tableWinnersWithOuters_iteration[kellogs].begin(), tableWinnersWithOuters_iteration[kellogs].end());
        tableWinnersWithOuters_test[kellogs].erase(tableWinnersWithOuters_test[kellogs].begin(), tableWinnersWithOuters_test[kellogs].end());
        tableWinnersWithOuters_outer[kellogs].erase(tableWinnersWithOuters_outer[kellogs].begin(), tableWinnersWithOuters_outer[kellogs].end());
    }

    /*!
    * \brief add the best individual to tableWinners
    *
    * */
    //
    void addWinner(uint64_t indiceTable, uint64_t kellogs) {
        vector<uint64_t> tmp(table[indiceTable]);
        tableWinners[kellogs].push_back(tmp);
    }

    /*!
    * \brief delete individual to be replaced by another one
    * the indice of the individual to be deleted is given as indiceTable
    * */
    void IndividualEraser(uint64_t indiceTable, bool debug = false) {
        uint64_t param_NKmer = table[0].size();
        if (debug) {
            cerr << "param_NKmer = " << param_NKmer << endl;
        }
        table.erase(table.begin() + indiceTable);

        //then add a new individual
        vector<uint64_t> popk;
        do{
            popk.clear();
            while (popk.size() < param_NKmer) {
                int pos_rand = rand() % dim_Kmers;
                bool state = false;
                for (uint64_t isec = 0; isec < popk.size(); isec++) {
                    if (popk[isec] == (uint64_t) pos_rand) {
                        state = true;
                        isec = popk.size();
                    }
                } //for isec
                if (!state) {
                    popk.push_back((uint64_t) pos_rand);
                }
            } //while popk.size()
           // cout<<"IndividualEraser-----"<<endl;
        }while(!checkifunique(popk,param_NKmer));
        //cout<<"eraserdebug :"<<popk.size()<<endl;

        table.push_back(popk);

        if (debug) {
            for (uint64_t d = 0; d < table.size(); d++) {
                for (uint64_t f = 0; f < param_NKmer; f++) {
                    cerr << "\t" << table[d][f];
                }
                cerr << endl;
            }
            cerr << endl;
        }
    }

    /*!
    * \brief apply a mutation of 1 kmer on 1 individual
    * the indice of the individual to be mutated is given as indiceTable
    * */
    void applyMutation_only1kmer(uint64_t indiceTable, bool debug = false) {
        bool stateLoop = true;
        while (stateLoop) {
            uint64_t newKmer = rand() % dim_Kmers;
            bool state = true;
            for (uint64_t isec = 0; isec < table[indiceTable].size(); isec++) {
                if (table[indiceTable][isec] == newKmer) {
                    state = false;
                    isec = table[indiceTable].size();
                }
            }

            if (state) {
                uint64_t pos_mutation = rand() % (table[indiceTable].size());
                vector <uint64_t>popk=table[indiceTable];
                popk[pos_mutation]=newKmer;
               // cout<<"applyMutation-----"<<endl;
                if (checkifunique(popk,table[indiceTable].size())==true){
                    table[indiceTable][pos_mutation] = newKmer;

                    //cout<<"mutdebug :"<<table[indiceTable].size()<<endl;

                    stateLoop = false;
                }
            }
        }
        if (debug) {
            for (uint64_t d = 0; d < table[indiceTable].size(); d++) {
                cerr << "\t" << table[indiceTable][d];
            }
            cerr << endl;
        }
        //std::sort(table[indiceTable].begin(), table[indiceTable].end());
    }

    /*!
    * \brief apply a mutation on 1 individual
    * the indice of the individual to be mutated is given as indiceTable
    * */
    void applyMutation(uint64_t indiceTable, float tauxMutation =0.5, bool debug = false) {
        //for each kmer of the individual we decide if we mutate it or not
        vector <uint64_t>popk = table[indiceTable];

        for(uint64_t kmer=0; kmer<table[indiceTable].size(); kmer++){
            float russianRoulette = (float) (rand() % 10000) / 10000;
            if (russianRoulette < tauxMutation) {
                //if the kmer has to be mutated, we pick another one randomly
                uint64_t newKmer = rand() % dim_Kmers;
                popk[kmer]=newKmer;

                if (debug) {
                    for (uint64_t d = 0; d < table[indiceTable].size(); d++) {
                        cerr << "\t" << table[indiceTable][d];
                    }
                    cerr << endl;
                }
            }
        }
        //if the individual is not unique we randomly change on of its kmer
        while (!checkifunique(popk,table[indiceTable].size())){
            uint64_t newKmer = rand() % dim_Kmers;
            uint64_t pos_mutation = rand() % (table[indiceTable].size());
            popk[pos_mutation] = newKmer;
        }
        table[indiceTable] = popk;
        //std::sort(table[indiceTable].begin(), table[indiceTable].end());
    }

   /*!
    * \brief apply a translocation between 2 individuals
    * the indice of the individuals to be translocated are given as indiceTable1&2
    * */
    bool applyTranslocation(uint64_t indiceTable1, uint64_t indiceTable2, bool debug = false) {


        if (indiceTable1 != indiceTable2) {
            vector <uint64_t> tmporga1;

                tmporga1 =table[indiceTable1];
                uint64_t originlength1 = table[indiceTable1].size();
                uint64_t originlength2 = table[indiceTable2].size();
                uint64_t length1 = rand() % (originlength1);
                uint64_t length2 = length1;
                uint64_t pos_mutation1 = rand() % (originlength1 - length1);
                uint64_t pos_mutation2 = rand() % (originlength2 - length2);

                if (debug) {
                    cerr << "with " << indiceTable1 << ", " << indiceTable2 << endl;
                    cerr << "\t" << originlength1 << "\t" << originlength2 << endl << "\t\t"
                            << length1 << "\t" << length2 << endl << "\t\t"
                            << pos_mutation1 << "\t" << pos_mutation2 << endl;
                }

                //vector<uint64_t> mem1;
                vector<uint64_t> mem2;

                //sauvegarde des indices kmers à etre transposés
                if (debug)
                    cerr << "\tbackups" << endl;
                //for (uint64_t i = pos_mutation1; i < (pos_mutation1 + length1); i++) {
                //    mem1.push_back(table[indiceTable1][i]);
                //}
                for (uint64_t i = pos_mutation2; i < (pos_mutation2 + length2); i++) {
                    mem2.push_back(table[indiceTable2][i]);
                }

                //suppression des kmers
                if (debug) {
                    cerr << "\terase" << endl;
                    for (uint64_t ideb = 0; ideb < originlength1; ideb++) {
                        cerr << "\t" << table[indiceTable1][ideb];
                    }
                    cerr << endl;
                    for (uint64_t ideb = 0; ideb < originlength1; ideb++) {
                        cerr << "\t" << table[indiceTable2][ideb];
                    }
                    cerr << endl;
                }
                for (uint64_t cpt = 0; cpt < length1; cpt++) {
                    tmporga1.erase(tmporga1.begin() + pos_mutation1);
                    //table[indiceTable2].erase(table[indiceTable2].begin() + pos_mutation2);
                }

                //ajout des kmers
                if (debug)
                    cerr << "\tswitching" << endl;
                /*for (uint64_t i = 0; i < mem1.size(); i++) {
                    bool state = true;
                    for (uint64_t pos = 0; pos < table[indiceTable2].size(); pos++) {
                        if (table[indiceTable2][pos] == mem1[i]) {
                            state = false;
                            pos = table[indiceTable2].size();
                        }
                    }
                    if (state) {
                        table[indiceTable2].push_back(mem1[i]);
                    }
                }*/
                for (uint64_t i = 0; i < mem2.size(); i++) {
                    bool state = true;
                    for (uint64_t pos = 0; pos < tmporga1.size(); pos++) {
                        if (tmporga1[pos] == mem2[i]) {
                            state = false;
                            pos = tmporga1.size();
                        }
                    }
                    if (state) {
                        tmporga1.push_back(mem2[i]);
                    }
                }

                //ajouts de nouveaux kmers pour compléter si besoin
                if (debug)
                    cerr << "\tadding" << endl;
                while (tmporga1.size() < originlength1) {
                    int pos_rand = rand() % dim_Kmers;
                    bool state = false;
                    for (uint64_t isec = 0; isec < tmporga1.size(); isec++) {
                        if (tmporga1[isec] == (uint64_t) pos_rand) {
                            state = true;
                            isec = tmporga1.size();
                        }
                    } //for isec
                    if (!state) {
                        tmporga1.push_back((uint64_t) pos_rand);
                    }
                }


                    if (!checkifunique(tmporga1,tmporga1.size())){
                        return false;
                     //   cout<<"translocation fail identic organism, tmporga1.size()="<<tmporga1.size()<<"originlength1="<<originlength1<<endl;

                    }else{
                        table[indiceTable1] = tmporga1;
                        return true;
                    }

            //std::sort(table[indiceTable1].begin(), table[indiceTable1].end());

            /*while (table[indiceTable2].size() < originlength2) {
                int pos_rand = rand() % dim_Kmers;
                bool state = false;
                for (uint64_t isec = 0; isec < table[indiceTable2].size(); isec++) {
                    if (table[indiceTable2][isec] == (uint64_t) pos_rand) {
                        state = true;
                        isec = table[indiceTable2].size();
                    }
                } //for isec
                if (!state) {
                    table[indiceTable2].push_back((uint64_t) pos_rand);
                }
            }*/


            /*if(debug){
                for(uint64_t d = 0; d<table[indiceTable1].size(); d++){
                    cerr << "\t" << table[indiceTable1][d];
                }
                cerr << endl;
                for(uint64_t d = 0; d<table[indiceTable2].size(); d++){
                    cerr << "\t" << table[indiceTable2][d];
                }
                cerr << endl;
            }*/
        }
        return false;
    }

    ///definition des partitions inner:modèle et outter:test du modèle à la fin uniquement.
    ///outter sera de percentage*plus bas effectif
    void createPartitionOutter(float percentage, bool debug = false) {
        uint64_t NPerGroup = uint64_t(percentage * (float) getMinGroup());
        if (debug)
            cerr << "createPartitionOutter : NPerGroup outter = " << NPerGroup << endl;
        dim_outter = NPerGroup;
        //recuperation d'une liste aleatoire des replicats marqués comme outter
        //pour chaque groupe
        for (map< string, uint64_t>::iterator it = dictocc.begin(); it != dictocc.end(); it++) {
            //selection des replicats cibles de façon aleatoire
            map<uint64_t, int> selection;
            uint64_t ks = 0;
            while (ks < NPerGroup) {
                uint64_t krandom = rand() % it->second + 1;
                if (selection.find(krandom) == selection.end()) {
                    selection[krandom] = 1;
                    ks++;
                }
            }
            //tag des replicats choisis
            for (uint64_t i = 0; i < initial_data_y.size(); i++) {
                if ((initial_data_y[i].compare(it->first) == 0) and (selection.find(initial_data_replicate[i]) != selection.end())) {
                    initial_data_state[i] = -1;
                }
            }
        }

        if(debug){
            cerr << "group\treplicate\tstate\tdictocc" << endl;
            for(uint64_t d = 0; d<initial_data_y.size() ; d++){
                cerr << initial_data_group[d] << "\t" <<
                    initial_data_replicate[d] << "\t" <<
                    initial_data_state[d] << "\t" <<
                    dictocc[initial_data_y[d]] << endl;
            }
        }
    }

    ///definition des partitions inner:modèle et outter:test du modèle à la fin uniquement.
    ///outter sera de percentage*plus bas effectif
    void readPartitionOutter(string param_pathToFile, bool debug = false) {

        ifstream file_afi(param_pathToFile.c_str(), ios::in);
        string content;
        uint64_t cptout=0;
        if (!file_afi) {
            cerr << "Error while opening " << param_pathToFile << " in read mode" << endl;
            exit(EXIT_FAILURE);
        }
        cout<<"initial_data_state.size "<<initial_data_state.size()<<endl;
        for (uint64_t i = 0; i < initial_data_state.size(); i++) {

            if (!getline(file_afi, content)) {
                if (!file_afi) {
                    cerr << "Error while opening " << param_pathToFile << " in read mode : not enough sample in the file : " << i << "/" << initial_data_state.size() << endl;
                    exit(EXIT_FAILURE);
                }
            }
            initial_data_state[i]=(stoi(content.c_str(), NULL));
            if (initial_data_state[i]==-1)
                cptout++;

            if (debug) {
                cerr <<i<<"\t"<< initial_data_state[i]<<endl;
            }
        }
        file_afi.close();


        vector<u_int64_t> uniquegrp(initial_data_group);
        auto it = unique (uniquegrp.begin(), uniquegrp.end());   // 10 20 30 20 10 ?  ?  ?  ?
        //uniquegrp.resize( distance(uniquegrp.begin(),it) );
        dim_outter = cptout/distance(uniquegrp.begin(),it);
        cout << "readPartitionOutter : NPerGroup outter = " << dim_outter << endl;
        //tag des replicats choisis

    }

    ///definition des partitions train test, basé sur le pourcentage de test
    ///outter sera de percentage*plus bas effectif
    /// nbTest : nombre de cross validation
    // then load in buffer for binary file
    void createPartitionTrainTest(float percentageTest, uint64_t nbTest, bool debug = false) {
        uint64_t NPerGroup = uint64_t(percentageTest * ((float) getMinGroup() - dim_outter));
        if (debug) {
            cerr << "createPartitionTrainTest : dim_outter = " << dim_outter << endl;
            cerr << "createPartitionTrainTest : diff = " << ((float) getMinGroup() - dim_outter) << endl;
            cerr << "createPartitionTrainTest : NPerGroup test = " << NPerGroup << endl;
        }

        //remise à 0, sauf si initialement destiné à la condition outer (-1)
        for (uint64_t i = 0; i < initial_data_y.size(); i++) {
            if (initial_data_state[i] != -1) {
                initial_data_state[i] = 0;
            }
        }
        if (debug) {
            cerr << "\tgroup\treplicate\tstate\tdictocc" << endl;
            for (uint64_t d = 0; d < initial_data_y.size(); d++) {
                cerr << "d = " << d << "\t" << initial_data_group[d] << "\t" <<
                        initial_data_replicate[d] << "\t" <<
                        initial_data_state[d] << "\t" <<
                        dictocc[initial_data_y[d]] << endl;
            }
        }
        //pour chaque test
        for (unsigned int ntest=0;ntest<nbTest ; ntest++){
            //pour chaque groupe
            for (map< string, uint64_t>::iterator it = dictocc.begin(); it != dictocc.end(); it++) {
                //selection des replicats cibles de façon aleatoire
                if (debug)
                    cerr << "\t" << it->first << "\t" << it->second << endl;
                map<uint64_t, int> selection;

                //on ajoute les outers
                for (uint64_t i = 0; i < initial_data_state.size(); i++) {
                    if ((initial_data_state[i] == -1) and (initial_data_y[i].compare(it->first) == 0)) {
                        selection[initial_data_replicate[i]] = -1;
                        if (debug)
                            cerr << "\t\t add " << i << " in selection as outer as replicate " << initial_data_replicate[i] << endl;
                    }
                }

                //les tests
                uint64_t ks = 0;
                while (ks < NPerGroup) {
                    uint64_t krandom = rand() % it->second + 1;
                    if (selection.find(krandom) == selection.end()) {
                        if (debug)
                            cerr << "\t\t add " << krandom << " in selection as candidate as trainer " << endl;
                        selection[krandom] = 1;
                        //cout<<"add test "<<krandom<<"\tlabel ="<<pow(2,ntest)<<endl;
                        ks++;
                    }
                }

                //tag des replicats choisis
                for (uint64_t i = 0; i < initial_data_y.size(); i++) {
                    if ((initial_data_y[i].compare(it->first) == 0) and (selection.find(initial_data_replicate[i]) != selection.end()) and (initial_data_state[i] != -1)) {
                        initial_data_state[i]  += pow(2,ntest);
                    }
                }
            }
        }
        if (debug) {
            cerr << "group\treplicate\tstate\tdictocc" << endl;
            for (uint64_t d = 0; d < initial_data_y.size(); d++) {
                cerr << initial_data_group[d] << "\t" <<
                        initial_data_replicate[d] << "\t" <<
                        initial_data_state[d] << "\t" <<
                        dictocc[initial_data_y[d]] << endl;
            }
        }
        if(initial_binary_file._self_isInBinary){
            initial_binary_file.getValuesFromTable(table); //chargement tampon
        }
    }

    ///prepare the data for FIFO layer
    map<uint64_t, vector< vector<long double>>> getTestFifo(bool debug = false) {
        if (debug) {
            cerr << "ClassPopulation::getTestFifo table.size =  " << table.size() << endl;
            cerr << "ClassPopulation::getTestFifo initial_data_y.size =  " << initial_data_y.size() << endl;
        }

        map<uint64_t, vector< vector<long double>>> result;
        for (uint64_t guy = 0; guy < table.size(); guy++) {
            vector< vector<long double>> vvtmp;
            if (debug)
                cerr << "guy = " << guy << endl;

            for (uint64_t i = 0; i < initial_data_y.size(); i++) {
                if (initial_data_state[i] >= 1) {
                    if (debug)
                        cerr << i<<"\t is a test " << endl;
                    //i is a test sample
                    vector<long double> vtmp;
                    for (uint64_t kmer_c = 0; kmer_c < table[0].size(); kmer_c++) {
                        if(initial_binary_file._self_isInBinary){
                            vtmp.push_back(initial_binary_file.getValue(i,table[guy][kmer_c],"getTestFifo"));
                        }
                        else{
                            vtmp.push_back(initial_data_x[i][ table[guy][kmer_c] ]);
                        }
                    }
                    vtmp.push_back(initial_data_group[i]);
                    vtmp.push_back(initial_data_state[i]);
                    vtmp.push_back(guy);
                    vvtmp.push_back(vtmp);
                }
            }
            result[guy] = vvtmp;
        }
        if (debug) {
            cerr << "ClassPopulation::getTestFifo table.size =  " << table.size() << endl;
            cerr << "ClassPopulation::getTestFifo initial_data_y.size =  " << initial_data_y.size() << endl;
            cerr << "ClassPopulation::getTestFifo map.size =  " << result.size() << endl;
            cerr << "ClassPopulation::getTestFifo map.vector.size =  " << result[0].size() << endl;
            cerr << "ClassPopulation::getTestFifo map.vector.vector.size =  " << result[0][0].size() << endl;
        }
        return result;
    }

    ///prepare the data for FIFO layer
    map<uint64_t, vector< vector<long double>>> getTrainFifo(bool debug = false) {
        if (debug) {
            cerr << "ClassPopulation::getTestFifo table.size =  " << table.size() << endl;
            cerr << "ClassPopulation::getTestFifo initial_data_y.size =  " << initial_data_y.size() << endl;
        }

        map<uint64_t, vector< vector<long double>>> result;
        for (uint64_t guy = 0; guy < table.size(); guy++) {
            vector< vector<long double>> vvtmp;
            for (uint64_t i = 0; i < initial_data_y.size(); i++) {
                if (initial_data_state[i] == 0) {
                    //i is a train sample
                    vector<long double> vtmp;
                    for (uint64_t kmer_c = 0; kmer_c < table[0].size(); kmer_c++) {
                        if(initial_binary_file._self_isInBinary){
                            vtmp.push_back(initial_binary_file.getValue(i,table[guy][kmer_c],"getTrainFifo"));
                        }
                        else{
                            vtmp.push_back(initial_data_x[i][ table[guy][kmer_c] ]);
                        }
                    }
                    vtmp.push_back(initial_data_group[i]);
                    vtmp.push_back(0);
                    vtmp.push_back(guy);
                    vvtmp.push_back(vtmp);
                }
            }
            result[guy] = vvtmp;
        }
        if (debug) {
            cerr << "ClassPopulation::getTrainFifo table.size =  " << table.size() << endl;
            cerr << "ClassPopulation::getTrainFifo initial_data_y.size =  " << initial_data_y.size() << endl;
            cerr << "ClassPopulation::getTrainFifo map.size =  " << result.size() << endl;
            cerr << "ClassPopulation::getTrainFifo map.vector.size =  " << result[0].size() << endl;
            cerr << "ClassPopulation::getTrainFifo map.vector.vector.size =  " << result[0][0].size() << endl;
        }
        return result;
    }

    ///prepare the data for FIFO layer
    map<uint64_t, vector< vector<long double>>> getTestOutterFifo(uint64_t kellogs,uint64_t begini,uint64_t endi, bool debug = false) {

        //preload kmer into buffer for binary file
        if(initial_binary_file._self_isInBinary){
            initial_binary_file.getValuesFromTable(tableWinners[kellogs],begini,endi);
        }

        map<uint64_t, vector< vector<long double>>> result;
        for (uint64_t guy = begini; guy < endi; guy++) {
            vector< vector<long double>> vvtmp;
            if (debug)
                cerr << "guy = " << guy << endl;

            for (uint64_t i = 0; i < initial_data_y.size(); i++) {
                if (initial_data_state[i] == -1) {
                    if (debug)
                        cerr << "\ti is a good one " << endl;
                    //i is a outer sample
                    vector<long double> vtmp;
                    for (uint64_t kmer_c = 0; kmer_c < tableWinners[kellogs][0].size(); kmer_c++) {
                        if(initial_binary_file._self_isInBinary){
                            vtmp.push_back(initial_binary_file.getValue(i,tableWinners[kellogs][guy][kmer_c],"getTestOutterFifo" ));
                        }
                        else{
                            vtmp.push_back(initial_data_x[i][tableWinners[kellogs][guy][kmer_c]]);
                        }
                    }
                    vtmp.push_back(initial_data_group[i]);
                    vtmp.push_back(-1);
                    vtmp.push_back(guy);
                    vvtmp.push_back(vtmp);
                }
            }
            result[guy] = vvtmp;
        }
        if (debug) {
            cerr << "ClassPopulation::getTestFifo tableWinners.size =  " << tableWinners[kellogs].size() << endl;
            cerr << "ClassPopulation::getTestFifo initial_data_y.size =  " << initial_data_y.size() << endl;
            cerr << "ClassPopulation::getTestFifo map.size =  " << result.size() << endl;
            cerr << "ClassPopulation::getTestFifo map.vector.size =  " << result[0].size() << endl;
            cerr << "ClassPopulation::getTestFifo map.vector.vector.size =  " << result[0][0].size() << endl;
        }
        return result;
    }

    ///prepare the data for FIFO layer
    map<uint64_t, vector< vector<long double>>> getTrainOutterFifo(uint64_t kellogs,uint64_t begini,uint64_t endi, bool debug = false) {
        //remise à 0
        for (uint64_t i = 0; i < initial_data_y.size(); i++) {
            if (initial_data_state[i] != -1) {
                initial_data_state[i] = 0;
            }
        }

        map<uint64_t, vector< vector<long double>>> result;
        for (uint64_t guy = begini; guy < endi; guy++) {
            vector< vector<long double>> vvtmp;
            for (uint64_t i = 0; i < initial_data_y.size(); i++) {
                if (initial_data_state[i] == 0) {
                    //i is a train sample
                    vector<long double> vtmp;
                    for (uint64_t kmer_c = 0; kmer_c < tableWinners[kellogs][0].size(); kmer_c++) {
                        if(initial_binary_file._self_isInBinary){
                            vtmp.push_back(initial_binary_file.getValue(i,tableWinners[kellogs][guy][kmer_c],"getTrainOutterFifo" ));
                        }
                        else{
                            vtmp.push_back(initial_data_x[i][tableWinners[kellogs][guy][kmer_c]]);
                        }
                    }
                    vtmp.push_back(initial_data_group[i]);
                    vtmp.push_back(0);
                    vtmp.push_back(guy);
                    vvtmp.push_back(vtmp);
                }
            }
            result[guy] = vvtmp;
        }
        if (debug) {
            cerr << "ClassPopulation::getTrainFifo tableWinners.size =  " << tableWinners[kellogs].size() << endl;
            cerr << "ClassPopulation::getTrainFifo initial_data_y.size =  " << initial_data_y.size() << endl;
            cerr << "ClassPopulation::getTrainFifo map.size =  " << result.size() << endl;
            cerr << "ClassPopulation::getTrainFifo map.vector.size =  " << result[0].size() << endl;
            cerr << "ClassPopulation::getTrainFifo map.vector.vector.size =  " << result[0][0].size() << endl;
        }
        return result;
    }

    ///prepare the data for FIFO layer
    map<uint64_t, vector< vector<long double>>> getTestOutterFifoCurrentG(vector<uint64_t>  resRandTreeRank, bool debug = false) {

        map<uint64_t, vector< vector<long double>>> result;
        uint64_t guy;
         for (uint64_t rankguy = 0; rankguy < table.size(); rankguy++) {

            guy = resRandTreeRank[rankguy];
            vector< vector<long double>> vvtmp;
            if (debug)
                cerr << "guy = " << guy << endl;

            for (uint64_t i = 0; i < initial_data_y.size(); i++) {
                if (initial_data_state[i] == -1) {
                    if (debug)
                        cerr << i<<"\t is a good one " << endl;
                    //i is a test sample
                    vector<long double> vtmp;
                    for (uint64_t kmer_c = 0; kmer_c < table[0].size(); kmer_c++) {
                        if(initial_binary_file._self_isInBinary){
                            vtmp.push_back(initial_binary_file.getValue(i,table[guy][kmer_c],"GetTestOutterFifoCurrentG"));
                        }
                        else{
                            vtmp.push_back(initial_data_x[i][table[guy][kmer_c]]);
                        }

                    }
                    vtmp.push_back(initial_data_group[i]);
                    vtmp.push_back(-2);
                    vtmp.push_back(rankguy);
                    vvtmp.push_back(vtmp);
                }
            }
            result[rankguy] = vvtmp;
        }


        return result;
    }

    ///prepare the data for FIFO layer
    map<uint64_t, vector< vector<long double>>> getTrainOutterFifoCurrentG(vector<uint64_t>  resRandTreeRank, bool debug = false) {
        //remise à 0
        for (uint64_t i = 0; i < initial_data_y.size(); i++) {
            if (initial_data_state[i] != -1) {
                initial_data_state[i] = 0;
            }
        }

        map<uint64_t, vector< vector<long double>>> result;
        uint64_t guy;
        for (uint64_t rankguy = 0; rankguy < table.size(); rankguy++) {

            guy = resRandTreeRank[rankguy];
            //cout<<"train rankguy "<<guy<<endl;
            vector< vector<long double>> vvtmp;
            for (uint64_t i = 0; i < initial_data_y.size(); i++) {
                if (initial_data_state[i] == 0) {
                    //i is a train sample
                    vector<long double> vtmp;
                    for (uint64_t kmer_c = 0; kmer_c < table[0].size(); kmer_c++) {
                        if(initial_binary_file._self_isInBinary){
                            vtmp.push_back(initial_binary_file.getValue(i,table[guy][kmer_c],"getTrainOutterFifoCurrentG"));
                        }
                        else{
                            vtmp.push_back(initial_data_x[i][table[guy][kmer_c]]);
                        }
                    }
                    vtmp.push_back(initial_data_group[i]);
                    vtmp.push_back(0);
                    vtmp.push_back(rankguy);
                    vvtmp.push_back(vtmp);
                }
            }
            result[rankguy] = vvtmp;
        }


        return result;
    }

    ///prepare the data for FIFO layer in case of computation for all individuals
    map<uint64_t, vector< vector<long double>>> getTrainOutterFifo_allPop(bool debug = false) {

        map<uint64_t, vector< vector<long double>>> result;
        for (uint64_t guy = 0; guy < table.size(); guy++) {
            vector< vector<long double>> vvtmp;
            for (uint64_t i = 0; i < initial_data_y.size(); i++) {
                if (initial_data_state[i] != -1) {
                    //i is a train sample
                    vector<long double> vtmp;
                    for (uint64_t kmer_c = 0; kmer_c < table[0].size(); kmer_c++) {
                        if(initial_binary_file._self_isInBinary){
                            vtmp.push_back(initial_binary_file.getValue(i,table[guy][kmer_c],"getTrainOutterFifo_allPop"));
                        }
                        else{
                            vtmp.push_back(initial_data_x[i][table[guy][kmer_c]]);
                        }
                    }
                    vtmp.push_back(initial_data_group[i]);
                    vtmp.push_back(0);
                    vtmp.push_back(guy);
                    vvtmp.push_back(vtmp);
                }
            }
            result[guy] = vvtmp;
        }
        if (debug) {
            cerr << "ClassPopulation::getTrainOutterFifo_allPop table.size =  " << table.size() << endl;
            cerr << "ClassPopulation::getTrainOutterFifo_allPop initial_data_y.size =  " << initial_data_y.size() << endl;
            cerr << "ClassPopulation::getTrainOutterFifo_allPop map.size =  " << result.size() << endl;
            cerr << "ClassPopulation::getTrainOutterFifo_allPop map.vector.size =  " << result[0].size() << endl;
            cerr << "ClassPopulation::getTrainOutterFifo_allPop map.vector.vector.size =  " << result[0][0].size() << endl;
        }
        return result;
    }

    ///prepare the data for FIFO layer in case of computation for all individuals
    map<uint64_t, vector< vector<long double>>> getTestOutterFifo_allPop(bool debug = false) {

        map<uint64_t, vector< vector<long double>>> result;
        for (uint64_t guy = 0; guy < table.size(); guy++) {
            vector< vector<long double>> vvtmp;
            for (uint64_t i = 0; i < initial_data_y.size(); i++) {
                if (initial_data_state[i] == -1) {
                    //i is a train sample
                    vector<long double> vtmp;
                    for (uint64_t kmer_c = 0; kmer_c < table[0].size(); kmer_c++) {
                        if(initial_binary_file._self_isInBinary){
                            vtmp.push_back(initial_binary_file.getValue(i,table[guy][kmer_c],"getTestOutterFifo_allPop"));
                        }
                        else{
                            vtmp.push_back(initial_data_x[i][table[guy][kmer_c]]);
                        }
                    }
                    vtmp.push_back(initial_data_group[i]);
                    vtmp.push_back(1);
                    vtmp.push_back(guy);
                    vvtmp.push_back(vtmp);
                }
            }
            result[guy] = vvtmp;
        }
        if (debug) {
            cerr << "ClassPopulation::getTestOutterFifo_allPop table.size =  " << table.size() << endl;
            cerr << "ClassPopulation::getTestOutterFifo_allPop initial_data_y.size =  " << initial_data_y.size() << endl;
            cerr << "ClassPopulation::getTestOutterFifo_allPop map.size =  " << result.size() << endl;
            cerr << "ClassPopulation::getTrainOuttergetTestOutterFifo_allPopFifo_allPop map.vector.size =  " << result[0].size() << endl;
            cerr << "ClassPopulation::getTestOutterFifo_allPop map.vector.vector.size =  " << result[0][0].size() << endl;
        }
        return result;
    }


};





/*******************************************************************************/
/*int main(){

    ClassPopulation popA = ClassPopulation();
    string file_csv = "../InOUtMatrix_onlyInformativeKMers_ML.txt";
    popA.loadCSV(file_csv);
    popA.getMinGroup();
    popA.createTable(5, 2, true);
    popA.createPartitionOutter(0.17, false);
    popA.createPartitionTrainTest(0.3, true);
    popA.getTestFifo(true);
    popA.computeVariability();
    popA.printVariability("variabilityTest.txt");



    return EXIT_SUCCESS;
}*/
/*******************************************************************************/
