#include <cstdio>
#include <iostream>
#include <iomanip>
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
#include <unistd.h>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>

//#include <boost/filesystem.hpp>
#include <sys/stat.h>

#include <stdlib.h>
#include <time.h>
#include <chrono>
#include "ClassPopulation.cpp"
#include "ClassFIFO.cpp"
#include "ClassLog.cpp"

using namespace std;

/*!
 * \file ClassGeneticAlg.cpp
 * \brief class to instantiate genetic algorithm.
 * \class Class ClassGeneticAlg
 * The genetic algorithm is applied on a population of individuals that is instantiated in ClassPopulation
 *
 *
 *
*/


/*!
* \brief // Find the sorted indices of a given vector
* author : Yotam Gingold
* source : https://gist.github.com/yig/32fe51874f3911d1c612
*/
template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {
                return v[i1] < v[i2];
            });

    return idx;
}

class ClassGeneticAlg {
public:
    vector<string> _self_UID; //identifiant unique

    //NN parameters
    /// \param _self_NHiddenNeuron number of neurons in hidden layer
    /// \param _self_NGeneration number of AG generation
    /// \param _self_NGeneration_refreshLog every _self_NGeneration_refreshLog the logs are updated
    /// \param _self_NIndividuals number of individuals
    /// \param _self_NKmer number of k-mers per individual
    /// \param _self_n_rotative_test number of tests after sampling train test populations
    /// \param _self_Method vector of method names
    /// \param _self_loading_pop_file file to load directly a pre-fixed population
    /// \param _self_loading_kmersamples_file file to load directly a pre-fixed population of kmers
    /// \param _self_loading_outer_file
    vector<uint64_t> _self_NHiddenNeuron;
    uint64_t _self_NGeneration;
    uint64_t _self_NGeneration_refreshLog;
    uint64_t _self_NIndividuals;
    uint64_t _self_NKmer; //per individual
    uint64_t _self_n_rotative_test;
    vector<string> _self_Method; //method to describe ML algorithm
    string _self_loading_pop_file;
    string _self_loading_kmersamples_file;
    string _self_loading_outer_file;

    /// \param _self_resRandTree vector of fitness for each individuals, computed by ML layer
    /// \param _self_resRandTreeOuter vector of fitness for each individuals, computed by ML layer on outers if _self_computeOutersForPopulation istrue
    /// \param _self_resRandTreePvalue pvalues of fitness
    /// \param _self_resRandTreeRank ranks of fitness
    /// \param _self_lresRandTreeRank inverse ranks of fitness
    vector<long double> _self_resRandTree;
    vector<long double> _self_resRandTreeOuter;
    vector<long double> _self_resRandTreePvalue;
    vector<uint64_t> _self_resRandTreeRank;
    vector<uint64_t> _self_lresRandTreeRank;

    /// \param _self_algorithmType type of algorithm to see GA evolve
    string _self_algorithmType;

    /// \param _self_scoremaxtreshold if the fitness is upper this threshold the GA stops
    /// \param _self_tauxMutation1
    /// \param _self_tauxMutation2
    /// \param _self_tauxMutation3
    /// \param _self_tauxTranslocation1
    /// \param _self_tauxTranslocation2
    /// \param _self_tauxTranslocation3
    float _self_scoremaxtreshold;
    float _self_tauxMutation1;
    float _self_tauxMutation2;
    float _self_tauxMutation3;
    float _self_tauxTranslocation1;
    float _self_tauxTranslocation2;
    float _self_tauxTranslocation3;

    /// \param  _self_Pcmin minimal probability for crossover, see IAGA_probability_Ravindran method
    /// \param  _self_Pcmax maximal probability for crossover, see IAGA_probability_Ravindran method
    /// \param  _self_Pmmin minimal probability for mutation, see IAGA_probability_Ravindran method
    /// \param  _self_Pmmax maximal probability for mutation, see IAGA_probability_Ravindran method
    /// \param  _self_lambda coeff as constant, see IAGA_probability_Ravindran method
    float _self_Pcmin;
    float _self_Pcmax;
    float _self_Pmmin;
    float _self_Pmmax;
    float _self_lambda;

    /// \param _self_percentage_outter percentage of samples that are seen as outers
    /// \param _self_percentage_test  percentage of samples that are seen as test
    /// \param _self_kill_ratio  percentage of individuals that are about to be killed at each iteration according their rank
    float _self_percentage_outter;
    float _self_percentage_test;
    float _self_kill_ratio;

    /// \param _self_NElite number of elites (the N better fitness) that can't be modified
    /// \param _self_kSelection
    /// \param _self_restart_AG
    /// \param _self_loghistorycount
    uint64_t _self_NElite;
    uint64_t _self_kSelection;
    uint64_t _self_restart_AG;
    uint64_t _self_loghistorycount;

    /// \param _self_computeOutersForPopulation true/false if the good candidates must be stored. Required for convergence test
    /// \param _self_thresholdTestGoodIndividual the threshold on tests to consider an individual as good
    /// \param _self_thresholdOuterGoodIndividual the threshold on outers to consider an individual as good
    bool _self_computeOutersForPopulation; //envoi et recuperation des outers pour memoriser les bon elements
    float _self_thresholdTestGoodIndividual;
    float _self_thresholdOuterGoodIndividual;
    bool _self_computeConvergence;

    bool _self_compute_all_outters;
    bool _self_mutationmode_1kmeronly;
    bool _self_doLog;
    bool _self_dooutter;

    /// \param _self_killElite boolean to indicate if elite must be killed
    /// \param _self_NGeneration_killElite if _self_killElite is true, every _self_NGeneration_killElite the elite is killed
    bool _self_killElite;
    uint64_t _self_NGeneration_killElite;

    vector<string> _self_pathLog;
    string _self_pathData;

    /// \param _self_population the objet that deals with kmers quantification
    /// \param _self_fifo the object that deals with communications with ML backend
    /// \param _self_log boolean if considering output as logs
    /// \param _self_Nlogs number of methods that share the same outers and split in test/outers
    /// \param _self_ActiveLogs need of evolving the respective method (if convergence is not obtained)
    ClassPopulation _self_population;
    ClassFIFO _self_fifo;
    vector<ClassLog> _self_log;
    uint64_t _self_Nlogs; //number of ML models at the same time
    vector<bool> _self_ActiveLogs;

    ClassGeneticAlg(string pUID, string ppathdata, string ppathlog, string phidneuron, string pmethod, uint64_t pgen = 10, uint64_t pind = 10,
        uint64_t pkmer = 50, float pscoremax = 1.1, float ptmut1 = 0.25, float ptmut2 = 0.05, float ptmut3 = 0, float pttrans1 = 0.8, float pttrans2 = 0.5, float pttrans3 = 0.1, uint64_t pelite = 1, float kSelectionratio = 0, bool pdolog = true) {
        _self_compute_all_outters=false;
        _self_mutationmode_1kmeronly=false;
        _self_killElite=false;


        _self_NGeneration = pgen; // number og generations
        _self_NIndividuals = pind; // number of individuals
        _self_NKmer = pkmer; // each individual has _self_NKmer kmers
        _self_scoremaxtreshold = pscoremax; // stop if scoremax is reached [0,inf]
        _self_tauxMutation1 = ptmut1; // mutation rate 1 [0,1]
        _self_tauxMutation2 = ptmut2; // mutation rate 2 [0,1]
        _self_tauxMutation3 = ptmut3; // mutation rate 3 [0,1]
        _self_tauxTranslocation1 = pttrans1; // translocation rate 1 [0,1]
        _self_tauxTranslocation2 = pttrans2; // translocation rate 2 [0,1]
        _self_tauxTranslocation3 = pttrans3; // translocation rate 3 [0,1]
        _self_NElite = pelite; //number of individuals that are considered as elite

        _self_kSelection = round(kSelectionratio * pind);
        _self_n_rotative_test=1;

        _self_pathData = ppathdata;
        //                                              => method as vector
        vector<string> contentMethod;
        boost::split(contentMethod, pmethod, boost::is_any_of(","));
        for (uint64_t ic = 0; ic < contentMethod.size(); ic++) {
            _self_Method.push_back(contentMethod[ic]);
        }
        //                                              => number of hidden neurons as vector
        vector<string> contentHidNeur;
        boost::split(contentHidNeur, phidneuron, boost::is_any_of(","));
        for (uint64_t ic = 0; ic < contentHidNeur.size(); ic++) {
            _self_NHiddenNeuron.push_back(stoi(contentHidNeur[ic]));
        }
        //                                              => path to Logs as vector
        vector<string> contentLog;
        boost::split(contentLog, ppathlog, boost::is_any_of(","));
        for (uint64_t ic = 0; ic < contentLog.size(); ic++) {
            _self_pathLog.push_back(contentLog[ic]);
        }

        ClassGeneticAlgConstruc2(pUID, pdolog);
    }

    ClassGeneticAlg(string pUID, string configfile, bool pdolog = true) {
        _self_n_rotative_test=1;
        _self_kill_ratio=0.5;
        _self_compute_all_outters=false;
        _self_mutationmode_1kmeronly=false;

        _self_algorithmType="IAGA";

        _self_NGeneration = 10; // number of generations
        _self_NGeneration_refreshLog=10000;// number of generations for log save frequency
        _self_NIndividuals = 200; // number of individuals
        _self_NKmer = 5; // each individual has _self_NKmer kmers
        _self_scoremaxtreshold = 1.1; // stop if scoremax is reached [0,inf]
        _self_tauxMutation1 = 0.25; // mutation rate 1 [0,1]
        _self_tauxMutation2 = 0.05; // mutation rate 2 [0,1]
        _self_tauxMutation3 = 0; // mutation rate 3 [0,1]
        _self_tauxTranslocation1 = 0.8; // translocation rate 1 [0,1]
        _self_tauxTranslocation2 = 0.5; // translocation rate 2 [0,1]
        _self_tauxTranslocation3 = 0.1; // translocation rate 3 [0,1]
        _self_NElite = 1; //number of individuals that are considered as elite

        _self_kSelection = 0;
        _self_percentage_outter=0.17;
        _self_percentage_test= 0.3;


        _self_killElite=false;
        _self_computeOutersForPopulation = false;
        _self_computeConvergence = false;

        readConfig(configfile);
        ClassGeneticAlgConstruc2(pUID, pdolog);


    }

    ClassGeneticAlg() {
    }

    void ClassGeneticAlgConstruc2(string pUID, bool pdolog = true) {
        ///Pour l'instant fixe;
        ///@todo : mise en place dans configuration
        cout << "---------ClassGeneticAlgConstruc2-------------------" << endl;
        cout << "Individuals=" << _self_NIndividuals << endl;

        cout << "Generation =" << _self_NGeneration << endl;
        cout << "kmer=" << _self_NKmer << endl;
        cout << "hiddenNeuron=";
        for (unsigned int i = 0; i < _self_NHiddenNeuron.size(); i++)
            cout << " " << _self_NHiddenNeuron[i];
        cout << endl << "kselection=" << _self_kSelection << endl;
        cout << "pathData=" << _self_pathData << endl;
        cout << "pathLog=";
        for (unsigned int i = 0; i < _self_pathLog.size(); i++)
            cout << " " << _self_pathLog[i];
        cout << endl << "method=";
        for (unsigned int i = 0; i < _self_Method.size(); i++)
            cout << " " << _self_Method[i];
        cout << endl << "elite=" << _self_NElite << endl;

        cout << "AlgorithmType=" << _self_algorithmType <<endl;
        /*cout << "mutationRate1=" << _self_tauxMutation1 << endl;
        cout << "mutationRate2=" << _self_tauxMutation2 << endl;
        cout << "mutationRate3=" << _self_tauxMutation3 << endl;
        cout << "TranslocationRate1=" << _self_tauxTranslocation1 << endl;
        cout << "TranslocationRate2=" << _self_tauxTranslocation2 << endl;
        cout << "TranslocationRate3=" << _self_tauxTranslocation3 << endl;*/
        cout << "scoreMaxTreshold=" << _self_scoremaxtreshold << endl;
        cout << "outterRate=" << _self_percentage_outter << endl;
        cout << "testRate=" << _self_percentage_test << endl;
        cout << "nRotativeTest=" << _self_n_rotative_test << endl;
        cout << "killRatio=" << _self_kill_ratio << endl;

        cout << "killelite=" << _self_killElite << endl;
        cout << "generationkillingelite=" << _self_NGeneration_killElite << endl;

        srand(time(NULL));
        _self_loghistorycount = 0;
        _self_restart_AG = 0;
        _self_doLog = pdolog;
        // several models at the same time can be asked => UID as vector
        vector<string> contentID;
        boost::split(contentID, pUID, boost::is_any_of(","));
        for (uint64_t ic = 0; ic < contentID.size(); ic++) {
            _self_UID.push_back(contentID[ic]);
        }
        _self_Nlogs = contentID.size();
        vector<bool> tmpAcLo(_self_Nlogs,true);
        _self_ActiveLogs = tmpAcLo;


        /************************************************
         * Create population and initiate table of results
        ************************************************/
        ClassPopulation tmp     = ClassPopulation();
        _self_population        = tmp;
        _self_resRandTree       = vector<long double> (_self_NIndividuals, 0.0);
        _self_resRandTreeOuter  = vector<long double> (_self_NIndividuals, 0.0);
        _self_resRandTreeRank   = vector<uint64_t> (_self_NIndividuals, 0);
        _self_lresRandTreeRank  = vector<uint64_t> (_self_NIndividuals, 0);
        _self_resRandTreePvalue = vector<long double> (_self_NIndividuals, 0.0);


        cerr << "ClassGeneticAlg::loading file" << endl;
        auto loadCSV_start = chrono::high_resolution_clock::now();

        /************************************************
         * Loading data
        ************************************************/
        if (_self_loading_kmersamples_file.empty())
            _self_population.loadCSV(_self_pathData, true);
        else
            _self_population.loadCSVwithConstraint(_self_pathData, _self_loading_kmersamples_file, 0, true);

        auto loadCSV_end = chrono::high_resolution_clock::now();
        string loadCSV_elapsed = to_string(chrono::duration_cast<std::chrono::milliseconds>(loadCSV_end - loadCSV_start).count());
        cerr << "loadCSV_elapsed time : " << loadCSV_elapsed << endl;
        //cerr <<"sizeof(_self_population) : "<<sizeof(_self_population.initial_data_x)*sizeof(_self_population.initial_data_x[0])*sizeof(_self_population.initial_data_x[0][0])<<endl;

        if (_self_UID.size() == 1) {
            if (_self_loading_pop_file.empty()) {
                cerr << "ClassGeneticAlg::creating table" << endl;
                _self_population.createTable(_self_NKmer, _self_NIndividuals);
            }
            else {
                cout << "ClassGeneticAlg::read individu table from file" << _self_loading_pop_file << endl;
                _self_population.readTable(_self_NKmer, _self_NIndividuals, _self_loading_pop_file);

                //Get restart number from last existing folder
                string cmd = "ls -d " + _self_pathLog[0] + "Dir/*/| grep -o  [0-9]*_[0-9]* > tmp" + _self_UID[0];
                system(cmd.c_str());

                ifstream file_afi((string("tmp") + _self_UID[0]).c_str(), ios::in);
                string line;
                vector <string> content;
                vector <u_int64_t> startn;
                if (!file_afi) {
                    cerr << "Error while opening " << ("tmp" + _self_UID[0]).c_str() << " in read mode" << endl;
                    exit(EXIT_FAILURE);
                }
                while (getline(file_afi, line)) {
                    cout << line << endl;
                    boost::split(content, line, boost::is_any_of("_"));
                    startn.push_back(stold(content[0]));
                    cout << "filepath=" << content[0] << endl;
                }
                file_afi.close();
                remove((string("tmp") + _self_UID[0]).c_str());

                if (!startn.empty()) {
                    sort(startn.begin(), startn.end());
                    _self_restart_AG = startn.back() + 1;
                    cout << "-- Restart searching, number " << _self_restart_AG << " --" << endl;

                }

            }
        }
        else {
            cerr << "ClassGeneticAlg::read indiv table from file bestindiv.csv" << endl;
            _self_population.readTable(_self_NKmer, _self_NIndividuals, "bestindiv.csv");
        }
        _self_dooutter = false;

        /************************************************
         * Initiate the logs
        ************************************************/
        vector<string> lf; //list of features to add in logs
        lf.push_back("scoremin");
        lf.push_back("nbElem");
        lf.push_back("bestSequences");
        lf.push_back("scoreHDin");
        lf.push_back("scoremax");
        //lf.push_back("time");
        lf.push_back("bestIndices");
        lf.push_back("scoremoy10");
        lf.push_back("method");
        lf.push_back("executionTimeTotal");
        lf.push_back("executionTimeConsumer");
        lf.push_back("ML_write_elapsed");
        for (uint64_t ilog = 0; ilog < _self_Nlogs; ilog++) {
            ClassLog Kellogs = ClassLog();
            Kellogs.init(lf, _self_NGeneration);
            _self_log.push_back(Kellogs);
        }
    }

    void applyOutter(float percentage, bool debug = false) {
        cerr << "ClassGeneticAlg::applying outter" << endl;
        if ((percentage > 0) and (percentage < 1)) {
            _self_dooutter = true;
            if (_self_loading_outer_file.empty())
                _self_population.createPartitionOutter(percentage, debug);
            else
                _self_population.readPartitionOutter(_self_loading_outer_file, debug);

            _self_log[0].writeOuter(_self_pathLog[0] + "Dir/", _self_population.initial_data_state);
        } else {
            cerr << "Error percentage of outers must be >0 and <1" << endl;
            exit(EXIT_FAILURE);
        }
    }

    //one step
    void checkrequestdone(string uid){
        if (_self_compute_all_outters){
            string filefinal = uid+"_request";
            bool state = true;


            while(state){
                ifstream fai(filefinal.c_str(), ios::in);
                state = !fai.fail();
                if(state){
                    fai.close();
                    this_thread::sleep_for(std::chrono::milliseconds(100));
                }

            }
        }

    }

    /**
     * \method evolve makes a round of ML and applies modifications on population according results
     * */
    long double evolve(uint64_t kellogs, int64_t ntry) {
        //cerr << "ClassGeneticAlg::evolve" << endl;
        checkrequestdone(_self_UID[kellogs]);
        // envoi des donnees a la couche FIFO pour calcul
        //cerr << "ClassGeneticAlg::evolve::Writing ML" << endl;
        // waiting previous request (outter) to be treated

        auto ML_beforewrite = chrono::high_resolution_clock::now();
        //_self_population.loadKmerValues();
        _self_fifo.writeMLbinary(_self_UID[kellogs], _self_population.getTrainFifo(false), _self_population.getTestFifo(false), _self_NIndividuals, _self_population.getDimSampleAnalysis(), _self_NKmer);

        // recuperation des resultats dans ResRandTree
        //cerr << "ClassGeneticAlg::evolve::retrieving results" << endl;
        auto ML_start = chrono::high_resolution_clock::now();
        _self_resRandTree = _self_fifo.readML(_self_UID[kellogs]);
        auto ML_end = chrono::high_resolution_clock::now();
        string ML_elapsed = to_string(chrono::duration_cast<std::chrono::milliseconds>(ML_end - ML_start).count());
        string ML_write_elapsed = to_string(chrono::duration_cast<std::chrono::milliseconds>(ML_start - ML_beforewrite).count());

        //if _self_computeOutersForPopulation then compute the outer scores for each individual
        if(_self_computeOutersForPopulation){
            //cerr << "compute outers for all the population" << endl;
            _self_fifo.writeMLbinary(_self_UID[kellogs], _self_population.getTrainOutterFifo_allPop(false), _self_population.getTestOutterFifo_allPop(false), _self_NIndividuals, _self_population.initial_data_y.size(), _self_NKmer);
            _self_resRandTreeOuter = _self_fifo.readML(_self_UID[kellogs]);
        }


        // calcul de resRandTreeRank
        // cerr << "ClassGeneticAlg::evolve::Ranking" << endl;
        vector<long double> ResRandTreeSort(_self_resRandTree);
        sort(ResRandTreeSort.begin(), ResRandTreeSort.end());
        for (uint64_t i = 0; i < _self_NIndividuals; i++) {
            vector<long double>::iterator k = find(ResRandTreeSort.begin(), ResRandTreeSort.end(), _self_resRandTree[i]);
            auto position = k - ResRandTreeSort.begin();
            _self_resRandTreeRank[i] = ResRandTreeSort.size() - position;
        }

        //update the store of good individuals by population instance
        if(_self_computeOutersForPopulation){
            //cerr << "ClassGeneticAlg::evolve::Treating good individuals" << endl;
            _self_population.updateOuterTest(kellogs, _self_resRandTree, _self_resRandTreeOuter, (long double) _self_thresholdTestGoodIndividual, (long double) _self_thresholdOuterGoodIndividual, _self_NGeneration-ntry+1);
        }

        if (_self_compute_all_outters){
            //compute all outterscrores
            //rank lineaire
            //cout<<"linear sorting"<<endl;
            vector<long double> lResRandTreeSort (_self_resRandTree);
            //vector<uint64_t> ResRandTreeSortguyid (_self_resRandTreeRank);
            int ii=0;
            for (auto i: sort_indexes(_self_resRandTree)) {
              //cout <<"ResRandTreeSort"<< _self_resRandTree[i] <<" , i="<<i<< ",ii"<<ii<<endl;
                _self_lresRandTreeRank[ii]=i;//inverse rank
                ii++;
            }

            _self_fifo.writeLastMLbinary(_self_UID[kellogs], _self_population.getTrainOutterFifoCurrentG(_self_lresRandTreeRank, false),
                            _self_population.getTestOutterFifoCurrentG(_self_lresRandTreeRank, false),
                            _self_population.initial_data_y.size(), _self_NKmer);

        }

        // calcul des pvalues
        //cerr << "ClassGeneticAlg::evolve::Pvalues" << endl;
        long double valuemin = ResRandTreeSort[0];
        //cout<<"valuemin "<<valuemin<<endl;
        long double valuemiddle = ResRandTreeSort[uint64_t(_self_NIndividuals / 2)];
        long double gamma = valuemin - 2 * valuemiddle;
        long double alpha = 0.0;
        for (uint64_t guy = 0; guy < _self_NIndividuals; guy++) {
            alpha += _self_resRandTree[guy] + gamma;
        }
        alpha = 1 / alpha;
        long double beta = gamma*alpha;
        //cout <<"alpha:"<<alpha<<" ; beta:"<<beta<<endl;
        for (uint64_t guy = 0; guy < _self_NIndividuals; guy++) {
            _self_resRandTreePvalue[guy] += alpha * _self_resRandTree[guy] + beta;
        }

        long double result = 0.0;
        for (uint64_t guy = 0; guy < _self_NIndividuals; guy++) {
            if (_self_resRandTree[guy] > result) {
                result = _self_resRandTree[guy];
            }
        }

        //Logs
        if (_self_doLog) {
            long double top10mean = 0.0;
            for (uint64_t i = ResRandTreeSort.size() - 1; i > ResRandTreeSort.size() - 11; i--) {
                top10mean += ResRandTreeSort[i];
            }
            top10mean /= 10;
            _self_log[kellogs].addFeatureToToken("scoremoy10", to_string(top10mean));
            _self_log[kellogs].addFeatureToToken("executionTimeConsumer", ML_elapsed);
            _self_log[kellogs].addFeatureToToken("ML_write_elapsed", ML_write_elapsed);

        }

        return result;
    }


    /*!
    * \brief compute probabilities of mutation rate & crossover
    * based on publication : Improved adaptive genetic algorithm with sparsity constraint applied to thermal neutron CT reconstruction of two phase flow
    * Mingfei Yan et al 2018 Meas. Sci. Technol.
    * */
    float IAGA_probability(long double f, long double fmin, long double favg, long double fmax, float rate1, float rate2, float rate3){
        //calculates the probability of mutation or translocation of an Improved Adaptative Genetic Algorithm
        float p=0.5;
        if (f < favg){
            p=(rate1*(favg-f)+rate2*(f-fmin))/(favg-fmin);
        } else {
            p=(rate2*(fmax-f)+rate3*(f-favg))/(fmax-favg);
        }
        return p;
    }


    /*!
    * \brief compute probabilities of mutation rate & crossover
    * based on publication : A Novel Clinical Decision Support System Using Improved Adaptive Genetic Algorithm for the Assessment of Fetal Well-Being
    * Sindhu Ravindran, Asral Bahari Jambek, Hariharan Muthusamy, and Siew-Chin Neoh,
    * Computational and Mathematical Methods in Medicine 2015
    * */
    float IAGA_probability_Ravindran(long double f, long double fmax, long double favg, long double Pmax, long double Pmin, long double lambda){
        //calculates the probability of mutation or translocation of an Improved Adaptative Genetic Algorithm
        float p = Pmax;
        if(f>=favg){
            p = (Pmax - (Pmax - Pmin)) / (1 + exp(lambda * ((f-favg)/(fmax-favg))));
        }
        return p;
    }

    float AGA_probability(long double f, long double fmin, long double favg, long double fmax, float rate1){
        //calculates the probability of mutation or translocation of an Adaptative Genetic Algorithm
        float p=0.5;
        if (f < favg){
            p=rate1;
        } else {
            p=rate1*((fmax-f)/(fmax-favg));
        }
        return p;
    }

    void geneticDisorder(int64_t ntry) {
        /////////////////////////////////////////////////////////
        // deaths & births
        cerr << endl << "ClassGeneticAlg::geneticDisorder" << endl;
        int killcount = 0;
        for (uint64_t guy = 0; guy < _self_NIndividuals; guy++) {
            uint64_t rankguy = _self_resRandTreeRank[guy];
            if (rankguy > _self_NElite) {
                long double pval_cum = 0.0;
                for (uint64_t i = 0; i < _self_NIndividuals; i++) {
                    if (_self_resRandTreeRank[i] <= rankguy) {
                        pval_cum += _self_resRandTreePvalue[i];
                    }
                }

                long double russianRoulette = (long double) (rand() % 1000) / 1000;
                //cout<<guy<<","<<russianRoulette<<","<<pval_cum<<endl;
                //cout<<pval_cum<<","<<_self_resRandTree[guy]<<","<<_self_resRandTreeRank[guy]<<endl;
                if (russianRoulette <= pval_cum*(_self_kill_ratio*2)) {
                    _self_population.IndividualEraser(guy);
                    //_self_resRandTreeRank[guy]=_self_NIndividuals;
                    killcount++;
                    _self_resRandTree.erase(_self_resRandTree.begin()+guy);
                    _self_resRandTree.push_back(0);

                }
            }

        }


        //cout << "killcount indiv: " << killcount << endl;
        //resort population new indiv last
        vector<long double> ResRandTreeSort(_self_resRandTree);
        //ResRandTreeSort =_self_resRandTree;
        vector<uint64_t> ResRandTreeSortguyid(_self_resRandTreeRank);
        //ResRandTreeSortguyid =_self_resRandTreeRank;
        long double fmax = 0;
        long double favg = 0;
        long double fmin = 1;
        int64_t total = 0;
        int ii = 0;
        for (auto i : sort_indexes(ResRandTreeSort)) {
            ResRandTreeSortguyid[ii] = i;
            _self_resRandTreeRank[i] = _self_NIndividuals - ii;
            //cout <<"ResRandTreeSort:"<< ResRandTreeSort[i] <<" , i="<<i<<"_self_resRandTreeRank[i]"<<_self_resRandTreeRank[i]<< endl;
            if (_self_resRandTree[i]>0){ //removing all the newly created individuals from the average
                favg+=_self_resRandTree[i];
                total++;
            }
            if (_self_resRandTree[i]>fmax)
                fmax=_self_resRandTree[i];
            if (_self_resRandTree[i]<fmin && _self_resRandTree[i]>0)
               fmin=_self_resRandTree[i];

            ii++;
        }
        favg/=total;


        /////////////////////////////////////////////////////////
        // mutation
        //cerr << "ClassGeneticAlg::evolve::Mutating" << endl;
        float pm=_self_tauxMutation1; //mutation probability (of a kmer in a given individual)
        for (uint64_t guy = 0; guy < _self_NIndividuals; guy++) {
            //cout<<"guy : "<<guy<<" , resRandTreeRank: "<<_self_resRandTreeRank[guy]<<" , resRandTree: "<<_self_resRandTree[guy]<<endl;
            if (_self_resRandTreeRank[guy] > _self_NElite && _self_resRandTree[guy]>0) { //don't need to mutate if you're the elite or if you've just been created
                if (_self_algorithmType=="IAGA_Ravindran"){
                    pm = IAGA_probability_Ravindran(_self_resRandTree[guy], fmax, favg, _self_Pmmax, _self_Pmmin, _self_lambda);
                }
                else{
                    if (_self_algorithmType=="IAGA"){
                        //cerr << "mutation::IAGA done" << endl;
                        pm=IAGA_probability(_self_resRandTree[guy],fmin,favg,fmax,_self_tauxMutation1,_self_tauxMutation2,_self_tauxMutation3);
                    } else if (_self_algorithmType=="AGA"){
                        //cerr << "mutation::AGA done" << endl;
                        pm=AGA_probability(_self_resRandTree[guy],fmin,favg,fmax,_self_tauxMutation1);
                    } else {
                        //cerr << "mutation::GA done"<< endl;
                        pm=_self_tauxMutation1;
                    }
                }
                if (_self_mutationmode_1kmeronly){
                 // old mutation only 1 kmer mutated by individus
                    float russianRoulette = (float) (rand() % 10000) / 10000;
                    if (russianRoulette < pm)
                        _self_population.applyMutation_only1kmer(guy);
                }else{
                    _self_population.applyMutation(guy, pm);
                }
            }
        }



        /////////////////////////////////////////////////////////
        // translocation
        //cerr << "ClassGeneticAlg::evolve::Translocation" << endl;
        uint64_t guy;
        float pc=_self_tauxTranslocation1; //cross-over probability (on an individual)
        for (uint64_t i = 0; i < _self_NIndividuals; i++) {
            guy = ResRandTreeSortguyid[i];
            //cerr<<i<<" : "<<guy<<" , "<<_self_resRandTreeRank[guy]<<" , "<<_self_resRandTree[guy]<<endl;
            if (_self_resRandTreeRank[guy] > _self_NElite) {
                if (_self_resRandTree[guy]>0){//if the individual's just been created, keep default pc
                    //cerr << "_self_algorithmType:"<<_self_algorithmType<<" pc:"<<pc<<endl;

                    if (_self_algorithmType=="IAGA_Ravindran"){
                        pc = IAGA_probability_Ravindran(_self_resRandTree[guy], fmax, favg, _self_Pcmax, _self_Pcmin, _self_lambda);
                    }
                    else{
                        if (_self_algorithmType == "IAGA"){
                            //cerr << "translocation::IAGA done"<< endl;
                            pc=IAGA_probability(_self_resRandTree[guy],fmin,favg,fmax,_self_tauxTranslocation1,_self_tauxTranslocation2,_self_tauxTranslocation3);
                        } else if (_self_algorithmType=="AGA"){
                            //cerr << "translocation::AGA done" << endl;
                            pc=AGA_probability(_self_resRandTree[guy],fmin,favg,fmax,_self_tauxTranslocation1);
                        } else {
                            //cerr << "translocation::GA done"<< endl;
                            pc=_self_tauxTranslocation1;
                        }
                    }
                }

                float russianRoulette = (float) (rand() % 10000) / 10000;
                if (russianRoulette < pc) {
                    //choosing another guy
                    bool state = true;
                    while (state) {
                        uint64_t teammate;
                        if (_self_kSelection == 0) {
                            teammate = rand() % _self_NIndividuals;
                        } else {
                            uint64_t teammatetmp;
                            uint64_t selBestRank = _self_NIndividuals;
                            for (uint64_t k = 0; k < _self_kSelection; k++) {
                                teammatetmp = rand() % _self_NIndividuals;
                                if (_self_resRandTreeRank[teammatetmp] < selBestRank) {
                                    selBestRank = _self_resRandTreeRank[teammatetmp];
                                    teammate = teammatetmp;
                                }
                            }
                        }

                        //if(_self_resRandTreeRank[teammate] > _self_NElite){
                        //@todo: why an other russian roullette?
                        //float russianRoulette2 = (float)(rand()%10000)/10000;
                        //if(russianRoulette2<_self_tauxTranslocation1){
                        state = ! _self_population.applyTranslocation(guy, teammate);
                        //cout<<"teammate select: "<<teammate<<" rank :"<<_self_resRandTreeRank[teammate]<<endl;



                        //_self_population.checkifunique(guy);
                        //}
                        //}
                    }
                }
            }
        }


        /////////////////////////////////////////////////////////
        // delete elite
        //cerr << "ClassGeneticAlg::evolve::Killing elite" << endl;
        if(_self_killElite){
            if (((_self_NGeneration - (ntry - 1)) % _self_NGeneration_killElite) == 0){
                cerr << "ClassGeneticAlg::evolve::Killing elite" << endl;
                for (uint64_t guy = 0; guy < _self_NIndividuals; guy++) {
                    if (_self_resRandTreeRank[guy] < _self_NElite){
                        _self_population.IndividualEraser(guy);
                        _self_resRandTree.erase(_self_resRandTree.begin()+guy);
                        _self_resRandTree.push_back(0.0);
                    }
                }
            }
        }


    }

    //run
    /**
     * \brief start the genetic algorithm
     * */
    void run() {
        applyOutter(_self_percentage_outter);
        cerr << "ClassGeneticAlg::Run" << endl;
        uint64_t kreal = 0;

        //find the most variating kmers and output the sub matrix for pre-analysis
        //_self_population.MakeMatrixFromMostVariantKmers(_self_pathLog[0] + "Dir/mostVariantKmers.matrix",_self_pathLog[0] + "Dir/lessVariantKmers.matrix",1000);

        //for _self_NGeneration generation
        for (uint64_t ntry = _self_NGeneration; ntry > 0; ntry--) {

            bool stateEvolution=false; //stays at false if all the methods have converged
            for (uint64_t kellogs = 0; kellogs < _self_Nlogs; kellogs++) {
                if(_self_ActiveLogs[kellogs]){
                    stateEvolution = true;
                }
            }

            if(stateEvolution){
                auto cycle_start = chrono::high_resolution_clock::now();
                kreal++;

                /*********************************************
                * split train/test
                *********************************************/
                cerr << "ClassGeneticAlg::createPartitionTrainTest" << endl;
                _self_population.createPartitionTrainTest(_self_percentage_test, _self_n_rotative_test);
                //re-init _self_resRandTreePvalue
                for (uint64_t ii = 0; ii < _self_NIndividuals; ii++){
                    _self_resRandTreePvalue[ii] = 0.0;
                }

                /*********************************************
                * Evolution of population
                *********************************************/
                vector<long double> results(_self_Nlogs, 0.0);
                for (uint64_t kellogs = 0; kellogs < _self_Nlogs; kellogs++) {
                    if(_self_ActiveLogs[kellogs]){
                        //for each method that shares the same outers and test/train splits
                        results[kellogs] = evolve(kellogs, ntry);
                        if (kellogs == 0)
                            cerr << endl<<"try " << ntry << ":\n";
                        for (uint64_t pfb = 0; pfb < kellogs; pfb++)
                            cerr << "\t";
                        cerr << "winner test= " << setprecision(3) << results[kellogs] ;
                        //add the winner
                        for (uint64_t guy = 0; guy < _self_NIndividuals; guy++) {
                            if (_self_resRandTree[guy] == results[kellogs]) {
                                _self_population.addWinner(guy, kellogs);
                                guy = _self_NIndividuals;
                            }
                        }

                        //do Log
                        if (_self_doLog) {
                            bool indiceSecure = true;
                            _self_log[kellogs].addFeatureToToken("scoremax", to_string(results[kellogs]));
                            _self_log[kellogs].addFeatureToToken("nbElem", to_string(_self_NIndividuals));
                            long double scoremin = 1.0;
                            string bestindices = "";
                            string bestsequences = "";
                            for (uint64_t guy = 0; guy < _self_NIndividuals; guy++) {
                                if (_self_resRandTree[guy] < scoremin) {
                                    scoremin = _self_resRandTree[guy];
                                }
                                if (_self_resRandTree[guy] == results[kellogs]) {
                                    if (indiceSecure) {
                                        bestindices = to_string(_self_population.table[guy][0]);
                                        bestsequences = _self_population.initial_sequenceID[_self_population.table[guy][0]];
                                        for (uint64_t jj = 1; jj < _self_population.table[guy].size(); jj++) {
                                            bestindices += "," + to_string(_self_population.table[guy][jj]);
                                            bestsequences += "," + _self_population.initial_sequenceID[_self_population.table[guy][jj]];
                                        }
                                        indiceSecure = false;
                                    }
                                }
                            }
                            _self_log[kellogs].addFeatureToToken("scoremin", to_string(scoremin));
                            _self_log[kellogs].addFeatureToToken("bestIndices", bestindices);
                            _self_log[kellogs].addFeatureToToken("bestSequences", bestsequences);
                            _self_log[kellogs].addFeatureToToken("method", _self_Method[kellogs]);
                        } //if (_self_doLog)
                    }
                } //for(uint64_t kellogs = 0; kellogs<_self_Nlogs; kellogs++)

                _self_population.updateCounts(_self_resRandTree);
                //tests de fin
                long double valmax = 0.0;
                for (uint64_t kellogs = 0; kellogs < _self_Nlogs; kellogs++){
                    if(_self_ActiveLogs[kellogs]){
                        valmax += results[kellogs];
                    }
                }
                uint64_t nlogactive = 0;
                for (uint64_t kellogs = 0; kellogs < _self_Nlogs; kellogs++){
                    if(_self_ActiveLogs[kellogs]){
                        nlogactive++;
                    }
                }
                valmax /= nlogactive;
                //apply mutations, translocations, deads and killing elite
                geneticDisorder(ntry);

                auto cycle_end = chrono::high_resolution_clock::now();
                string cycle_elapsed = to_string(chrono::duration_cast<std::chrono::milliseconds>(cycle_end - cycle_start).count());
                for (uint64_t kellogs = 0; kellogs < _self_Nlogs; kellogs++) {
                    _self_log[kellogs].addFeatureToToken("executionTimeTotal", cycle_elapsed);
                }
                if (valmax >= _self_scoremaxtreshold) {
                    cout << endl << "DEBUG : _self_scoremaxtreshold exceeded" << endl;
                    ntry = -1;
                }
                else {
                    if ((_self_doLog) and (ntry > 1)) {
                        for (uint64_t kellogs = 0; kellogs < _self_Nlogs; kellogs++)
                            _self_log[kellogs].updateToken(ntry - 1);
                    }
                }
                if (((((_self_NGeneration - (ntry - 1)) % _self_NGeneration_refreshLog) == 0)&(ntry != _self_NGeneration)) | (ntry <= 1)) {
                    _self_loghistorycount++;
                    //cout<<"_self_loghistorycount"<<_self_loghistorycount<<"   "<<ntry<<endl;
                    for (uint64_t kellogs = 0; kellogs < _self_Nlogs; kellogs++) {
                        if(_self_ActiveLogs[kellogs]){
                            //cout<<"computeHDinScore"<<kellogs<<endl;
                            cerr << endl << "checkrequestdone";
                            checkrequestdone(_self_UID[kellogs]);
                            cerr << endl << "computeHDinScore";
                            computeHDinScore(kreal);

                            //@todo make dir
                            string newpathhist = _self_pathLog[kellogs] + "Dir/" + to_string(_self_restart_AG) + "_" + to_string(_self_loghistorycount) + "/";
                            cout << "New log in " << newpathhist << endl;
                            const int dir_err = mkdir(newpathhist.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                            if (-1 == dir_err){
                                cout<<"Error creating directory! : " << newpathhist << endl << "try again after sleep" << endl;
                                const int dir_err2 = system(("mkdir -p " + newpathhist).c_str());

                                sleep(2);
                                if (-1 == dir_err2)
                                {
                                cout<<"Error creating directory! : "<< newpathhist<<endl;
                                exit(1);
                                }

                            }

                            //writing files
                            cerr << endl << "writeContent";
                            _self_log[kellogs].writeContent(newpathhist + "resultAG.json");
                            cerr << endl << "writeOccKmer";
                            _self_log[kellogs].writeOccKmer(newpathhist + "resultAG.json", _self_population.initial_sequenceID, _self_population.initial_sequenceCounts, _self_population.initial_sequenceCountsScore);
                            cerr << endl << "writeTableOfIndividuals";
                            _self_log[kellogs].writeTableOfIndividuals(newpathhist, _self_population.table);
                            if(_self_computeOutersForPopulation){
                                cerr << endl << "writeGoodIndividuals";
                                _self_log[kellogs].writeGoodIndividuals(_self_pathLog[kellogs] + "Dir/goodIndividuals.csv", _self_population.initial_sequenceID,
                                                    _self_population.tableWinnersWithOuters_individuals[kellogs],
                                                    _self_population.tableWinnersWithOuters_iteration[kellogs],
                                                    _self_population.tableWinnersWithOuters_test[kellogs],
                                                    _self_population.tableWinnersWithOuters_outer[kellogs]);
                                _self_population.clearGoodIndividuals(kellogs);
                            }

                            
                            bool convergence_log_c = false;
                            if(((_self_NGeneration - ntry + 1) / 1000.0) >5){

                                if(_self_computeConvergence and _self_computeOutersForPopulation){
                                //if(false){
                                    //cout << endl << "DEBUG : computing convergence " << endl;
                                    convergence_log_c = _self_population.convergence(_self_pathLog[kellogs] + "Dir/goodIndividuals.csv", newpathhist + "HammingSuite.tsv");
                                    if(convergence_log_c){
                                        _self_ActiveLogs[kellogs] = false;
                                        cout << "DEBUG : convergence = " << convergence_log_c << endl;
                                    }
                                }
                            }
                            else{
                                //cout << endl << "DEBUG : not enough data to compute convergence " << endl;
                            }
                        }

                    }
                }
            } // if(stateEvolution)
        } //for(int64_t ntry = _self_NGeneration; ntry>0; ntry--)

        for (uint64_t kellogs = 0; kellogs < _self_Nlogs; kellogs++) {
            checkrequestdone(_self_UID[kellogs]);
            _self_fifo.killML(_self_UID[kellogs]);
        }


    }

    void computeHDinScore(uint64_t kreal) {
        //cerr << endl << "entry in computeHDinScore" << endl;
        //calcul des HDin sur l'ensemble des premiers de la classe
        for (uint64_t kellogs = 0; kellogs < _self_Nlogs; kellogs++) {
            bool statefinal = true;
            uint64_t kbegin;
            uint64_t kend;
            uint64_t kcount = 0; //_self_NGeneration_refreshLog*(_self_loghistorycount-1);

            vector<long double> lastVector;
            while (statefinal) {
                kbegin = _self_NGeneration_refreshLog * (_self_loghistorycount - 1)+ (kcount * _self_NIndividuals);
                kend = _self_NGeneration_refreshLog * (_self_loghistorycount - 1)+((kcount + 1) * _self_NIndividuals);
                //cout<<"kbegin"<<kbegin<<"kend"<<kend<<"kcount"<<kcount<<"kreal"<<kreal<<endl;


                if (kend >= kreal) {
                    kend = kreal;
                    statefinal = false;
                }


                //_self_population.loadKmerWinnersValues(kellogs, kbegin, kend);
                //cerr << "entry in computeHDinScore" << endl;
                _self_fifo.writeLastMLbinary(_self_UID[kellogs], _self_population.getTrainOutterFifo(kellogs, kbegin, kend, false), _self_population.getTestOutterFifo(kellogs, kbegin, kend, false),
                        _self_population.initial_data_y.size(), _self_NKmer);
                vector<long double> tmpVector = _self_fifo.readML(_self_UID[kellogs]);



                /*
                cerr << "group\treplicate\tstate\tdictocc" << endl;
                for (uint64_t d = 0; d < _self_population.initial_data_y.size(); d++) {
                    cerr << _self_population.initial_data_group[d] << "\t" <<
                            _self_population.initial_data_replicate[d] << "\t" <<
                            _self_population.initial_data_state[d] << "\t" <<
                            _self_population.dictocc[_self_population.initial_data_y[d]] << endl;
                }

                 */

                for (uint64_t jvm = 0; jvm < tmpVector.size(); jvm++) {
                    lastVector.push_back(tmpVector[jvm]);
                }
                kcount++;

            }


            if (_self_doLog) {
                vector<string> resRandTreeInString;
                for (uint64_t pp = 0; pp < lastVector.size(); pp++) {
                    resRandTreeInString.push_back(to_string(lastVector[pp]));
                }


                _self_log[kellogs].addFeatureToTokenMultiGeneration("scoreHDin", resRandTreeInString, _self_NGeneration - (_self_NGeneration_refreshLog * (_self_loghistorycount - 1)));
            }


        }
        //cerr << "out of computeHDinScore" << endl;
    }

    /*!
    * \brief function to read configuration file in yml format
    *
    *
    * */
    void readConfig(string conffile) {
    	cerr << "opening configuration file " << conffile << endl;
        string line;
        ifstream file_afi(conffile.c_str(), ios::in);

        if (!file_afi) {
            cerr << "Error while opening " << conffile << " in read mode" << endl;
            exit(EXIT_FAILURE);
        }
        //map<string, long double> catalog;
        vector <string> HiddenNeuron;


        while (getline(file_afi, line)) {
            if(!line.empty()){
            	cerr << "reading line : " << line << endl;

                vector<string> content;
                boost::split(content, line, boost::is_any_of("="));

                boost::trim(content[0]);
                if (content.size()>1){
                    boost::trim(content[1]);
                }else{
                    content.push_back("");
                }
                //cout<<"contents:"<<content[0]<<"||"<<content[1]<<endl;

                transform(content[0].begin(), content[0].end(), content[0].begin(), ::tolower);

                /******************************************************
                 * Population
                 * ***************************************************/
                if (content[0] == "individuals")
                    _self_NIndividuals = strtold(content[1].c_str(), NULL);

                if (content[0] == "generation")
                    _self_NGeneration = strtold(content[1].c_str(), NULL);

                if (content[0] == "kmer")
                    _self_NKmer = strtold(content[1].c_str(), NULL);


                if (content[0] == "kselection")
                    _self_kSelection = strtold(content[1].c_str(), NULL);

                /******************************************************
                 * Files
                 * ***************************************************/
                if (content[0] == "pathdata")
                    _self_pathData = content[1];

                if (content[0] == "pathlog"){
                    boost::split(_self_pathLog, content[1], boost::is_any_of(","));
                    for (unsigned int i = 0; i < _self_pathLog.size(); i++)
                        system(("mkdir -p " + _self_pathLog[i] + "Dir/").c_str());
                }

                /******************************************************
                 * ML method
                 * ***************************************************/
                if (content[0] == "hiddenNeuron") {
                    boost::split(HiddenNeuron, content[1], boost::is_any_of(","));
                    for (int i = 0; i < HiddenNeuron.size(); i++)
                        _self_NHiddenNeuron.push_back(strtold(HiddenNeuron[i].c_str(), NULL));
                }
                if (content[0] == "method")
                    boost::split(_self_Method, content[1], boost::is_any_of(","));
                if (content[0] == "elite")
                    _self_NElite = strtold(content[1].c_str(), NULL);

                /******************************************************
                 * Mutation & translocation
                 * ***************************************************/
                if (content[0] == "algorithmtype")
                    _self_algorithmType=content[1];
                if (content[0] == "mutationrate1")
                    _self_tauxMutation1 = strtof(content[1].c_str(), NULL);
                if (content[0] == "mutationrate2")
                    _self_tauxMutation2 = strtof(content[1].c_str(), NULL);
                if (content[0] == "mutationrate3")
                    _self_tauxMutation3 = strtof(content[1].c_str(), NULL);
                if (content[0] == "translocationrate1")
                    _self_tauxTranslocation1 = strtof(content[1].c_str(), NULL);
                if (content[0] == "translocationrate2")
                    _self_tauxTranslocation2 = strtof(content[1].c_str(), NULL);
                if (content[0] == "translocationrate3")
                    _self_tauxTranslocation3 = strtof(content[1].c_str(), NULL);
                if (content[0] == "pcmin")
                    _self_Pcmin = strtof(content[1].c_str(), NULL);
                if (content[0] == "pcmax")
                    _self_Pcmax = strtof(content[1].c_str(), NULL);
                if (content[0] == "pmmin")
                    _self_Pmmin = strtof(content[1].c_str(), NULL);
                if (content[0] == "pmmax")
                    _self_Pmmax = strtof(content[1].c_str(), NULL);
                if (content[0] == "lambda")
                    _self_lambda = strtof(content[1].c_str(), NULL);

                /******************************************************
                 * Genetic alg parameters
                 * ***************************************************/
                if (content[0] == "scoremaxthreshold")
                    _self_scoremaxtreshold = strtof(content[1].c_str(), NULL);
                if (content[0] == "outterrate")
                    _self_percentage_outter = strtof(content[1].c_str(), NULL);
                if (content[0] == "testrate")
                    _self_percentage_test = strtof(content[1].c_str(), NULL);
                if (content[0] == "generationrefreshlog")
                    _self_NGeneration_refreshLog = strtold(content[1].c_str(), NULL);
                if (content[0] == "startfrompopulation")
                    _self_loading_pop_file = content[1];
                if (content[0] == "startfromouter")
                    _self_loading_outer_file = content[1];
                if (content[0] == "startfromkmersubset")
                    _self_loading_kmersamples_file = content[1];
                if (content[0] == "nrotativetest")
                    _self_n_rotative_test = strtold(content[1].c_str(), NULL);
                if (content[0] == "killratio")
                    _self_kill_ratio = strtof(content[1].c_str(), NULL);
                if (content[0] == "computealloutters"){
                    if ((content[1]=="1")||(content[1]=="true"))
                        _self_compute_all_outters = true;
                    cout<<"_self_compute_all_outters : "<<_self_compute_all_outters<<endl;
                }
                if (content[0] == "mutationmode_1kmeronly"){
                    if ((content[1]=="1")||(content[1]=="true"))
                        _self_mutationmode_1kmeronly = true;
                    cout<<"_self_mutationmode_1kmeronly : "<<_self_mutationmode_1kmeronly<<endl;
                }

                //Killing elite
                if (content[0] == "killelite"){
                    if ((content[1]=="1")||(content[1]=="true"))
                        _self_killElite = true;
                }
                if (content[0] == "generationkillingelite")
                    _self_NGeneration_killElite = strtold(content[1].c_str(), NULL);


                /******************************************************
                * methods around good candidates
                ******************************************************/
                if (content[0] == "computeoutersforpopulation"){
                    if ((content[1]=="1")||(content[1]=="true")){
                        _self_computeOutersForPopulation = true;
                    }
                    else{
                        _self_computeOutersForPopulation = false;
                    }
                    cout << "_self_computeOutersForPopulation : " << _self_computeOutersForPopulation << " (store all good individuals)" << endl;
                }
                if (content[0] == "computeconvergence"){
                    if ((content[1]=="1")||(content[1]=="true")){
                        _self_computeConvergence = true;
                    }
                    else{
                        _self_computeConvergence = false;
                    }
                    if(!_self_computeOutersForPopulation){
                        _self_computeConvergence = false;
                    }
                    cout << "_self_computeConvergence : " << _self_computeConvergence << " (convergence of genetic algorithm)" << endl ;
                }


                if (content[0] == "thresholdindividualtest")
                    _self_thresholdTestGoodIndividual = strtof(content[1].c_str(), NULL);
                if (content[0] == "thresholdindividualouter")
                    _self_thresholdOuterGoodIndividual = strtof(content[1].c_str(), NULL);


            }
            /* if (catalog.find(content[0]) == catalog.end())
                 catalog[content[0]] = 0;
             catalog[content[0]] += strtold(content[1].c_str(), NULL);*/
        }
        cerr << "DEBUG : _self_computeOutersForPopulation = " << _self_computeOutersForPopulation << endl;
        cerr << "DEBUG : _self_thresholdTestGoodIndividual = " << _self_thresholdTestGoodIndividual << endl;
        cerr << "DEBUG : _self_thresholdOuterGoodIndividual = " << _self_thresholdOuterGoodIndividual << endl;

        file_afi.close();
    }
};

int main(int argc, char *argv[]) {
    //phidneuron = 10, pgen = 70, pind = 20, pkmer = 50

    string pUID, ppathdata = "", ppathlog = "", phidneuron = "5", pmethod = "MLP", configfile = "";
    uint64_t pgen = 100, pind = 10, pkmer = 5;
    float pscoremax = 1.1, kselratio = 0.25;

    ClassGeneticAlg prod;



    if (argc < 8) {
        if (argc == 3) {
            pUID = argv[1]; //vector comma separated
            configfile = argv[2];

            prod = ClassGeneticAlg(pUID, configfile);

        }
        else {
            cout << "ClassGeneticAlg have missing arguments";
            return EXIT_FAILURE;
        }

    }
    else {
        pUID = argv[1]; //vector comma separated
        ppathdata = argv[2];
        ppathlog = argv[3]; //vector comma separated
        phidneuron = argv[4]; //vector comma separated
        pmethod = argv[5]; //vector comma separated
        pgen = stoi(argv[6]);
        pind = stoi(argv[7]);
        pkmer = stoi(argv[8]);

        prod = ClassGeneticAlg(pUID, ppathdata, ppathlog, phidneuron, pmethod, pgen, pind, pkmer, pscoremax, kselratio);
        prod._self_percentage_outter = 0.17;
        prod._self_percentage_test = 0.3;
    }
    if (argc > 9)
        pscoremax = stof(argv[9]);



    //ClassGeneticAlg prod = ClassGeneticAlg(pUID,  ppathdata,  ppathlog,  phidneuron , pmethod ,  pgen ,  pind ,pkmer,pscoremax,kselratio);
    /*
    cout << "_self_NIndividuals:" << prod._self_NIndividuals << endl;

    cout << "_self_NGeneration =" << prod._self_NGeneration << endl;
    cout << "kmer=" << prod._self_NKmer << endl;
    cout << "hiddenNeuron=";
    for (int i = 0; i < prod._self_NHiddenNeuron.size(); i++)
        cout << " " << prod._self_NHiddenNeuron[i] << endl;
    cout << endl << "kselection=" << prod._self_kSelection << endl;
    cout << "pathData=" << prod._self_pathData << endl;
    cout << "pathLog=";
    for (int i = 0; i < prod._self_pathLog.size(); i++)
        cout << " " << prod._self_pathLog[i] << endl;
    cout << endl << "method=";
    for (int i = 0; i < prod._self_Method.size(); i++)
        cout << " " << prod._self_Method[0] << endl;
    cout << endl << "elite=" << prod._self_NElite << endl;
    cout << "mutationRate=" << prod._self_tauxMutation1 << endl;
    cout << "TranslocationRate1=" << prod._self_tauxTranslocation1 << endl;
    cout << "TranslocationRate2=" << prod._self_tauxTranslocation2 << endl;
    cout << "scoreMaxTreshold=" << prod._self_scoremaxtreshold << endl;
    cout << "outterRate=" << prod._self_percentage_outter << endl;
    cout << "testRate=" << prod._self_percentage_test << endl;

    cout << "startFromPopulation=" << prod._self_loading_pop_file << endl;
    cout << "startFromOuter=" << prod._self_loading_outer_file << endl;
    cout << "startfromKmerSubset=" << prod._self_loading_kmersamples_file << endl;
    cout << "nRotativeTest=" << prod._self_n_rotative_test << endl;
    */
    prod.run();
    cerr << endl;

    return EXIT_SUCCESS;
}
