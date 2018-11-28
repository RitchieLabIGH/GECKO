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
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/policies/policy.hpp>
#include "ClassEntropy.cpp"

using namespace std;

void RemoveNonInformativeKmers(string pInpFile, string pResFile, uint64_t distMax=10, double pvalueMax=0.5, long double pMIthreshold=0.5, uint64_t psizeBuffer=100000, bool pdebug=false){
    ifstream file_ai(pInpFile.c_str(), ios::in);
    ofstream file_ao(pResFile.c_str(), ios::out);

    if(!file_ai){
      cerr << "Error" << endl;
      exit(0);
    }
    if(!file_ao){
      cerr << "Error" << endl;
      exit(0);
    }
    string line;
    string header1, header2;

    //get the 2 lines of headers
    getline(file_ai, header1); //the names of the samples, note it begins with \t
    getline(file_ai, header2); //the groups, note the first item is 'group'
    file_ao << header1 << endl << header2 << endl;

    /**********************************
    //get all the groups
    **********************************/
    vector<string> content00;
    boost::split(content00, header2, boost::is_any_of("\t"));
    long double nbSample = content00.size()-1; //number of samples for each kmer
    uint64_t sizecontent = content00.size(); //total number of column (samples + KmerID in fact)

    uint64_t nbSampleDec = uint64_t(nbSample/10);
    //cout << nbSampleDec  << " pour nbSample = " << nbSample << endl;

    map<string, uint64_t> mapGroupsToSize; //will get the counts of each group
    map<string, vector<bool> > mapGroupsToMask; //mask of presence/abscence for each group
    map<string, vector<uint64_t> > mapGroupsToMapToThresholds;

    for(uint64_t i =1; i<sizecontent; i++){
      if(mapGroupsToSize.find(content00[i])==mapGroupsToSize.end()){
        mapGroupsToSize[content00[i]] = 0;
        vector<bool> mask_c(nbSample,false);
        mapGroupsToMask[content00[i]] = mask_c;
      }
      mapGroupsToSize[content00[i]]+=1;
      mapGroupsToMask[content00[i]][i-1] = true;
    }

    
    /**********************************
    //definition variables
    **********************************/
    vector<int> basis(nbSample,0);
    vector< vector<int> > data; // bunch of kmers
    vector<bool> dataindic; // true/false if the kmer is flagged as valid or not (true as default)
    vector<string> dataKmerID;
    uint64_t numKMers = 0; // current number of kmers inside data
    uint64_t kw = 0; //number of kmers trully written in output

    data.push_back(basis);
    dataindic.push_back(true);
    dataKmerID.push_back("");
    numKMers++;
    
    while (getline(file_ai, line)){
        vector<string> content;
        boost::split(content, line, boost::is_any_of("\t"));
        vector<int> vect_states_c(basis);
        uint64_t apm = 0;
        dataKmerID.push_back(content[0]);
        for(uint64_t i = 1; i<sizecontent; i++){
          vect_states_c[i-1] = stoi(content[i]);
        }

        data.push_back(vect_states_c);
        dataindic.push_back(true);
        numKMers++;
        
        /**********************************
        //OK a bunch to treat
        **********************************/
        if(numKMers == psizeBuffer){
            EntropyDecomposition ED = EntropyDecomposition();
            uint64_t i = 0;
            for(uint64_t i =0; i<numKMers-1; i++){
              uint64_t kwd = 0;
              if(dataindic[i]){
                #pragma omp parallel for
                for(uint64_t j = i+1; j<numKMers; j++){
                  if(dataindic[j]){
                    long double distHamming = 0; ///distHamming is the cumulative distance between 'bits'
                    //computation of hamming distance
                    uint64_t k = 0;
                    while(k<nbSample){                        
                      if(data[i][k] != data[j][k]){
                          distHamming++;
                          if(distHamming>distMax){
                              goto speed1;
                          }
                      }
                      k++;
                    }
                    speed1:
                    if(distHamming<=distMax){
                      dataindic[j] = false;
                      kwd++;                  
                    } //if(distHamming<=distMax)

                    if(dataindic[j]){
                        //on essaie l'information mutuelle
                        if(ED.SymmetricUncertaintyWithInt(data[i], data[j], false)>pMIthreshold){
                            dataindic[j] = false;
                        }
                    }
                    
                    
                  } //if(dataindic[j])

                } //for(uint64_t j =i; j<numKMers; j++)
                //if (pdebug) cerr << "\t" << i << ", kwd = " << kwd << endl;
              } //if(dataindic[i])

            } //for(uint64_t i =0; i<numKMers; i++)
            for(uint64_t i =1; i<numKMers; i++){
              if(dataindic[i]){
                file_ao << dataKmerID[i];
                for(uint64_t j = 0; j<nbSample;j++){
                  file_ao << "\t" << data[i][j];
                }
                file_ao << endl;
                kw++;
              }
            }
            numKMers = 0;
            data.clear();
            dataindic.clear();
            dataKmerID.clear();
            data.push_back(basis);
            dataindic.push_back(true);
            dataKmerID.push_back("");
            numKMers++;
            cerr << "kw = " << kw << endl;
            //exit(0);
      } //if(numKMers == psizeBuffer) */
    } //while (getline(file_ai, line))
    
    
    
    EntropyDecomposition ED = EntropyDecomposition();
    uint64_t i = 0;
    for(uint64_t i =0; i<numKMers-1; i++){
      uint64_t kwd = 0;
      if(dataindic[i]){
        #pragma omp parallel for
        for(uint64_t j = i+1; j<numKMers; j++){
          if(dataindic[j]){
            vector<bool> status_c(nbSample,false); ///status_c is the flags of each possible mismatch
            long double distHamming = 0; ///distHamming is the cumulative distance between 'bits'
            //computation of hamming distance
            uint64_t k = 0;
            while(k<nbSample){                        
              if(data[i][k] != data[j][k]){
                  distHamming++;
                  status_c[k] = true;
                  if(distHamming>distMax){
                      goto speed2;
                  }
              }
              k++;
            }
            speed2:
            if(distHamming<=distMax){
              dataindic[j] = false;
              kwd++;
            } //if(distHamming<=distMax)

            if(dataindic[j]){
                //on essaie l'information mutuelle
                if(ED.SymmetricUncertaintyWithInt(data[i], data[j], false)>pMIthreshold){
                    dataindic[j] = false;
                }
            }
            
            
          } //if(dataindic[j])

        } //for(uint64_t j =i; j<numKMers; j++)
        //if (pdebug) cerr << "\t" << i << ", kwd = " << kwd << endl;
      } //if(dataindic[i])

    } //for(uint64_t i =0; i<numKMers; i++)
    for(uint64_t i =1; i<numKMers; i++){
      if(dataindic[i]){
        file_ao << dataKmerID[i];
        for(uint64_t j = 0; j<nbSample;j++){
          file_ao << "\t" << data[i][j];
        }
        file_ao << endl;
        kw++;
      }
    }
    numKMers = 0;
    data.clear();
    dataindic.clear();
    dataKmerID.clear();
    data.push_back(basis);
    dataindic.push_back(true);
    dataKmerID.push_back("");
    numKMers++;
    cerr << "kw = " << kw << endl;

    file_ai.close();
    file_ao.close();
}



int main(int argc, char** argv){
  /*EntropyDecomposition ED = EntropyDecomposition(10);
  vector<string> A;
  vector<string> B;
  A.push_back("21");A.push_back("1");A.push_back("33");A.push_back("2");A.push_back("64");A.push_back("4");A.push_back("43");A.push_back("2");A.push_back("15");A.push_back("1");
  B.push_back("10");B.push_back("01");B.push_back("10");B.push_back("50");B.push_back("60");B.push_back("40");B.push_back("30");B.push_back("20");B.push_back("10");B.push_back("10");

  cout << ED.SymmetricUncertainty(A,B,false) << endl;
*/
  //RemoveNonInformativeKmers("cancer/smockers/InOUtMatrix_sample100M_discret_AmevaOwn.txt", "cancer/smockers/InOUtMatrix_sample100M_discret_Ameva_OnlyInformative.txt", 18, 0.1, 100000, true);
  RemoveNonInformativeKmers(argv[1], argv[2], stoi(argv[3]), 0.3, 0.7, 500000, true);

  return 0;
}
