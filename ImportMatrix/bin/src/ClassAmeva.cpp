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
#include <omp.h>


using namespace std;

//!
//! \file   IndexingKMers.cpp
//! \brief  The purpose is to handle a list of kmer tabulated files (ASCII)
//!         to create a matrix of the most explicative kmers
//! \author Aubin.Thomas William.Ritchie
//! \version    1.0
//! \date   Sept 2017



//!
//! \brief Naive class to deal with discretization of continuous variables according supervised MDL algorithm
//! author : Aubin Thomas
//!
class AmevaDecomposition{
public:
    uint64_t _self_Nmax; //the maximum number of thresholds
    vector <long double> _self_thresholds; // in paper : L
    vector <long double> _self_candidates; // in paper : B
    vector <bool> _self_candidates_indicator; // if the corresponding candidate has been chosen as threshold
    long double _self_GlobalAmeva; // in paper : GlobalAmeva

    AmevaDecomposition(uint64_t pNmax = 28){
        _self_Nmax = pNmax; //will be revised according the data as size-1 if smaller
        _self_GlobalAmeva = 0;
    }

    long double ComputeAmevaScore(vector<long double> pvector_values, vector<string> pvector_classes, vector<string> pvector_classes_listing, vector<long double> thresholds_tmp, bool pdebug){
        uint64_t k = thresholds_tmp.size()-1;
        uint64_t l = pvector_classes_listing.size();
        uint64_t N = pvector_values.size();

        if(pdebug){
            cerr << "\t\t";
            cerr << "k : " << k << "\t";
            cerr << "l : " << l << "\t";
            cerr << "N : " << N << endl;
        }
        if(pdebug){
            cerr << "\t\tpvector_values : ";
            for(uint64_t idebug = 0; idebug < pvector_values.size(); idebug++){
                cerr << "  " << pvector_values[idebug];
            }
            cerr << endl;
            cerr << "\t\tpvector_classes : ";
            for(uint64_t idebug = 0; idebug < pvector_classes.size(); idebug++){
                cerr << "  " << pvector_classes[idebug];
            }
            cerr << endl;
            cerr << "\t\tpvector_classes_listing : ";
            for(uint64_t idebug = 0; idebug < pvector_classes_listing.size(); idebug++){
                cerr << "  " << pvector_classes_listing[idebug];
            }
            cerr << endl;
            cerr << "\t\tthresholds_tmp : ";
            for(uint64_t idebug = 0; idebug < thresholds_tmp.size(); idebug++){
                cerr << "  " << thresholds_tmp[idebug];
            }
            cerr << endl;

        }

        long double result = 0.0;
        for(uint64_t i = 0; i<l; i++){
            for(uint64_t j = 0; j<k; j++){
                long double ni  = 0;
                long double nj  = 0;
                long double nij = 0;

                for(uint64_t isample = 0; isample<N; isample++){
                    bool si = false;
                    bool sj = false;
                    if(pvector_classes[isample].compare(pvector_classes_listing[i])==0){
                        si = true;
                        ni++;
                    }
                    if((pvector_values[isample]>=thresholds_tmp[j]) and (pvector_values[isample]<=thresholds_tmp[j+1])){
                        sj = true;
                        nj++;
                    }
                    if((si==true) and (sj==true)){
                      nij++;
                    }
                }
                if(pdebug){
                  cerr << "\t\t " << i << "x" << j << "\t";
                  cerr << "ni : " << ni << "\t";
                  cerr << "nj : " << nj << "\t";
                  cerr << "nij : " << nij << endl;
                }
                result+=(nij*nij)/(ni*nj);
            }
        }
        result = (long double) N*(-1+result);
        result/= (long double)(k*(l-1));

        return result;
    }

    /** Discretization
    Discretization : apply Ameva discretization on vector of values according vector of classes
    **/
    vector<string> Discretization(vector<long double> pvector_values, vector<string> pvector_classes, vector<string> pvector_classes_listing, bool pdebug){
        if ((pvector_values.size()-1)<_self_Nmax) _self_Nmax = (pvector_values.size()-1);
        _self_GlobalAmeva = 0;

        if(pdebug){
            cerr << "_self_GlobalAmeva = " << _self_GlobalAmeva << endl;
            cerr << "pvector_values : ";
            for(uint64_t idebug = 0; idebug < pvector_values.size(); idebug++){
                cerr << "  " << pvector_values[idebug];
            }
            cerr << endl;
            cerr << "pvector_classes : ";
            for(uint64_t idebug = 0; idebug < pvector_classes.size(); idebug++){
                cerr << "  " << pvector_classes[idebug];
            }
            cerr << endl;
        }
        // get vector of values with order and unique procedures
        vector<long double> pvector_values_ord (pvector_values); // ordered vector
        sort(pvector_values_ord.begin(), pvector_values_ord.end());
        vector<long double> pvector_values_ord_unique; //distinct values of ordered values
        for(uint64_t i = 0; i < (pvector_values_ord.size()-1); i++){
            if(pvector_values_ord[i+1] != pvector_values_ord[i]){
                pvector_values_ord_unique.push_back(pvector_values_ord[i]);
            }
        }
        pvector_values_ord_unique.push_back(pvector_values_ord[(pvector_values_ord.size()-1)]);
        if(pdebug){
            cerr << "pvector_values_ord_unique : ";
            for(uint64_t idebug = 0; idebug < pvector_values_ord_unique.size(); idebug++){
                cerr << "  " << pvector_values_ord_unique[idebug];
            }
            cerr << endl;
        }

        // get init thresholds
        long double Do = pvector_values_ord_unique[0];
        long double Dk = pvector_values_ord_unique[pvector_values_ord_unique.size()-1];
        _self_thresholds.push_back(Do);
        _self_thresholds.push_back(Dk);
        if(pdebug){
            cerr << "_self_thresholds : ";
            for(uint64_t idebug = 0; idebug < _self_thresholds.size(); idebug++){
                cerr << "  " << _self_thresholds[idebug];
            }
            cerr << endl;
        }

        // prepare the candidate boundaries
        for(uint64_t i = 0; i < (pvector_values_ord_unique.size()-1); i++){
            long double candidate = (pvector_values_ord_unique[i]+pvector_values_ord_unique[i+1])/2;
            _self_candidates.push_back(candidate);
            _self_candidates_indicator.push_back(false); //by default all the candidate boudaries have not been chosen
        }
        if(pdebug){
            cerr << "_self_candidates : ";
            for(uint64_t idebug = 0; idebug < _self_candidates.size(); idebug++){
                cerr << "  " << _self_candidates[idebug];
            }
            cerr << endl;
        }

        // Rock N'Roll
        uint64_t k=1;
        bool stateRun = true;
        while(stateRun){
            if(pdebug){
              cerr << endl << "k = " << k << endl;
            }
            long double Ameva_max = 0;
            uint64_t rank_max = -1;

            for(uint64_t icandidate = 0; icandidate < _self_candidates.size(); icandidate++)
            {
                if(!_self_candidates_indicator[icandidate]){
                    //calcul de la valeur Ameva avec cet indicateur
                    vector<long double> thresholds_tmp (_self_thresholds);
                    thresholds_tmp.push_back(_self_candidates[icandidate]);
                    sort(thresholds_tmp.begin(), thresholds_tmp.end());
                    long double Ameva_c = ComputeAmevaScore(pvector_values, pvector_classes, pvector_classes_listing, thresholds_tmp, pdebug);
                    if(pdebug){
                      cerr << "\t\t" << "icandidate = " << icandidate << "\tAmeva_c = " << Ameva_c << endl;
                    }
                    if(Ameva_c>Ameva_max){
                        Ameva_max = Ameva_c;
                        rank_max  = icandidate;
                    }
                    if(pdebug){
                      cerr << "\t\t" << "Ameva_max = " << Ameva_max << "\trank_max = " << rank_max << endl;
                    }
                }
            }
            if(Ameva_max>_self_GlobalAmeva){
                _self_GlobalAmeva = Ameva_max;
                _self_candidates_indicator[rank_max] = true;
                _self_thresholds.push_back(_self_candidates[rank_max]);
                sort(_self_thresholds.begin(), _self_thresholds.end());
                if(pdebug){
                  cerr << "\t\t_self_thresholds : ";
                  for(uint64_t idebug = 0; idebug < _self_thresholds.size(); idebug++){
                      cerr << "  " << _self_thresholds[idebug];
                  }
                  cerr << endl;
                }
                k++;
            }
            else{
                stateRun = false;
            }
            if(k>_self_Nmax){
                stateRun = false;
            }
        }

        //transform into discrete states
        vector<string> result;
        if(pdebug){
          cerr << "\n\n_self_thresholds : ";
          for(uint64_t idebug = 0; idebug < _self_thresholds.size(); idebug++){
              cerr << "  " << _self_thresholds[idebug];
          }
          cerr << endl;
        }
        for(uint64_t i = 0; i < pvector_values.size(); i++){
            uint64_t rank=0;
            for(int j = 0; j<(_self_thresholds.size()-1); j++){
                if((pvector_values[i]>=_self_thresholds[j]) and (pvector_values[i]<=_self_thresholds[j+1])){
                    rank = j;
                    j = _self_thresholds.size();
                }
            }
            result.push_back(to_string(rank));
        }
        return result;
    }

}; //class


/*
int main(){
    vector<long double> test;
    test.push_back(1);
    test.push_back(55);
    test.push_back(2);
    test.push_back(345);
    test.push_back(512);
    test.push_back(404);
    test.push_back(2);
    test.push_back(3);
    test.push_back(5);
    test.push_back(4);
    test.push_back(23);
    test.push_back(2);

    vector<string> cat;
    cat.push_back("A");
    cat.push_back("A");
    cat.push_back("A");
    cat.push_back("A");
    cat.push_back("A");
    cat.push_back("A");
    cat.push_back("B");
    cat.push_back("B");
    cat.push_back("B");
    cat.push_back("B");
    cat.push_back("B");
    cat.push_back("B");

    vector<string> cat_listing;
    cat_listing.push_back("A");
    cat_listing.push_back("B");

    AmevaDecomposition obj = AmevaDecomposition(20);
    vector<string> result = obj.Discretization(test,cat,cat_listing,true);
    for(int i = 0;i < result.size(); i++){
        cout << "  " << test[i];
    }
    cout << endl;

    for(int i = 0;i < result.size(); i++){
        cout << "  " << cat[i];
    }
    cout << endl;

    for(int i = 0;i < result.size(); i++){
        cout << "  " << result[i];
    }
    cout << endl;
    return 0;
}
*/
