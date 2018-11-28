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
#include "ClassProfileOfKmers.cpp"




/*
################################################################################
################################################################################
*/

//! \fn long double fixthresholds(string plistpath="", uint64_t psizeOfSample=100000000)
//! \brief  to fix thresholds between signal and noise
//!         between all the data, takes the group with the biggest number of duplicates
//!         for each kmer it counts the number of 0
//!         if a kmer has a number of 0 upper than #duplicates-int64_t(#duplicates/3), the other values are considered as artefacts
//!         between all artefactual events, the biggest value associated with at least 2% of events is considered as the threshold
//! \param plistpath list of files of kmers in 3 col format
//! \param psizeOfSample the number of kmers that are taken randomly in first duplicate to compute the artefactual events (default 100M)
//! \return (long double) the value that is considered as a threshold between signal and noise
//!
long double fixthresholds(string plistpath="", uint64_t psizeOfSample=100000000){
    //__________________________________________________________________
    //récupérer depuis la liste le groupe possédant le plus de duplicats
    ifstream file_ai(plistpath.c_str(), ios::in);
    string groupmax="";
    int64_t occmax = 0;
    string rp, rn, rg;
    if(file_ai){
      map<string, int64_t> groupOccurences;
      while(file_ai >> rp >> rn >> rg){
        if(groupOccurences.find(rg)==groupOccurences.end()){
          groupOccurences[rg] = 0;
        }
        groupOccurences[rg] += 1;
      }
      //recherche maximum
      for(map<string, int64_t>::iterator goit = groupOccurences.begin(); goit != groupOccurences.end(); goit++){
        string class_c = goit->first;
        int64_t number_c = goit->second;
        if(number_c > occmax){
          occmax = number_c;
          groupmax = class_c;
        }
      }
    }
    else{
      cerr << "ERROR WHILE OPENING " << plistpath << " IN READ MODE" << endl;
      return -1;
    }
    file_ai.close();
    cerr << "fixthresholds: best group is " << groupmax << " with a number of occmax = " << occmax << " duplicates" << endl;

    //__________________________________________________________________
    //création de l'objet et sélectionner les fichiers correspondants
    int64_t limitOfFiles = 6;
    if(occmax> limitOfFiles){
        occmax = limitOfFiles;
    }
    ProfileOfKmers result = ProfileOfKmers(occmax, 30, psizeOfSample);
    int64_t kf=0;
    ifstream file_ai2(plistpath.c_str(), ios::in);
    while((file_ai2 >> rp >> rn >> rg) and (kf<occmax)){
      if(rg.compare(groupmax)==0){
        kf++;
        cerr << rp << ", " << rn << ", " << rg << " added" << endl;
        result.addFile(rp, rn, rg);
      }
    }
    file_ai2.close();

    //__________________________________________________________________
    //prendre un echantillon de ces données et sortir matrice
    result.getAleaData();

    //__________________________________________________________________
    //a travers les comptages regarder combien de comptages correspondent à des 0 dans les autres replicats
    int64_t numOf0required = occmax-int64_t((long double)occmax/3); //nombre de 0 pour dire que les autres occurences sont du bruit de fond
    map<long double, uint64_t> histogram; //frequence d'association bruit de fond / valeurs
    uint64_t allevent = 0; //toutes les fois que du bruit de fond est détecté

    //pour chaque kmer
    for (map < uint64_t, vector<long double> >::iterator it=result._data.begin(); it!=result._data.end(); ++it){
      //pour chaque duplicat
      uint64_t kmerc = it->first;
      vector<long double> datatmp;
      int64_t nbzerotmp = 0;
      for(int64_t i = 0; i<occmax; i++){
        if(result._data[kmerc][i]==0){
          nbzerotmp++;
        }
        else{
          datatmp.push_back(result._data[kmerc][i]);
        }
      }
      if(nbzerotmp>=numOf0required){
        for(int64_t j=0; j<datatmp.size(); j++){
          allevent++;
          if(histogram.find(datatmp[j])==histogram.end()){
            histogram[datatmp[j]]=0;
          }
          histogram[datatmp[j]]+=1;
        }
      }
    }

    //__________________________________________________________________
    //valeur max en dessus de 2% par rapport a allevent, la plus grande est la meilleure
    long double maxvalue = -1;
    for(map<long double, uint64_t>::iterator it = histogram.begin(); it!=histogram.end(); it++){
      long double valuec = it->first;
      uint64_t occc = it->second;
      cout << "HIST\t" << valuec << "\t" << occc << "\t" << allevent << "\t" << (long double) occc/allevent << endl;
      if((((long double) occc/allevent)>=0.02) and (valuec>maxvalue)){
        maxvalue = valuec;
      }
    }
    return maxvalue;
}


//! \fn void makeBigMatrix(string plist, string pout, long double pthreshold)
//! \brief import all kmers, and considering only counts >=pthreshold
//! \param plist list of files of kmers in 3 col format
//! \param pout the matrix file
//! \param pthreshold the threshold between signal and noise. Can be given by fixthresholds function
//!
void makeBigMatrix(string plist, string pout, long double pthreshold){
    //voir combien de fichiers composent la liste plist
    ifstream file_ai(plist.c_str(), ios::in);
    int64_t nbfiles = 0;
    string rp, rn, rg;
    if(file_ai){
      map<string, int64_t> groupOccurences;
      while(file_ai >> rp >> rn >> rg){
        nbfiles+=1;
      }
      file_ai.close();
    }
    else{
      exit(4651);
    }
    ProfileOfKmers result = ProfileOfKmers(nbfiles, 30, 0);
    result.addListOfFiles(plist);
    result.getMassiveData(pthreshold, pout, true);
}

//! \fn void getDiscretizationMDL(string pathListGroups, string pathInpFile, string pathResFile, uint64_t sizeBuff)
//! \brief apply a discretization of kmers along groups
//! \param plist list of files of kmers in 3 col format
//! \param pathInpFile matrix file that have been generated by makeBigMatrix
//! \param pathResFile output matrix
//! \param sizeBuff ????
//!
void getDiscretizationMDL(string plist, string pathInpFile, string pathResFile, uint64_t sizeBuff){
    //voir combien de fichiers composent la liste plist
    ifstream file_ai(plist.c_str(), ios::in);
    int64_t nbfiles = 0;
    string rp, rn, rg;
    if(file_ai){
      map<string, int64_t> groupOccurences;
      while(file_ai >> rp >> rn >> rg){
        nbfiles++;
      }
      file_ai.close();
      cerr << "Debug : nbfiles = " << nbfiles << endl;
    }
    else{
      exit(4651);
    }
    cerr << "Debug : create object" << endl;
    ProfileOfKmers result = ProfileOfKmers(nbfiles, 30, 0);
    result.addListOfFiles(plist);
    cerr << "Debug : launching discretization" << endl;
    result.discretizationMDL(pathInpFile, pathResFile, sizeBuff, false);
}


//! \fn void getDiscretizationAmeva(string pathListGroups, string pathInpFile, string pathResFile, uint64_t sizeBuff)
//! \brief apply a discretization of kmers along groups
//! \param plist list of files of kmers in 3 col format
//! \param pathInpFile matrix file that have been generated by makeBigMatrix
//! \param pathResFile output matrix
//! \param sizeBuff ????
//!
void getDiscretizationAmeva(string plist, string pathInpFile, string pathResFile, uint64_t sizeBuff){
  //voir combien de fichiers composent la liste plist
  ifstream file_ai(plist.c_str(), ios::in);
  int64_t nbfiles = 0;
  string rp, rn, rg;
  if(file_ai){
    map<string, int64_t> groupOccurences;
    while(file_ai >> rp >> rn >> rg){
      nbfiles++;
    }
    file_ai.close();
    cerr << "Debug : nbfiles = " << nbfiles << endl;
  }
  else{
    exit(4651);
  }
  cerr << "Debug : create object" << endl;
  ProfileOfKmers result = ProfileOfKmers(nbfiles, 30, 0);
  result.addListOfFiles(plist);
  cerr << "Debug : launching discretization" << endl;
  result.discretizationAmeva(pathInpFile, pathResFile, sizeBuff, false);
}


//! \fn void getDiscretizationR(string pathListGroups, string pathInpFile, string pathResFile, uint64_t sizeBuff)
//! \brief apply a discretization of kmers along groups
//! \param plist list of files of kmers in 3 col format
//! \param pathInpFile matrix file that have been generated by makeBigMatrix
//! \param pathResFile output matrix
//! \param sizeBuff ????
//!
void getDiscretizationR(string plist, string pathInpFile, string pathResFile, uint64_t sizeBuff){
    //voir combien de fichiers composent la liste plist
    ifstream file_ai(plist.c_str(), ios::in);
    int64_t nbfiles = 0;
    string rp, rn, rg;
    if(file_ai){
      map<string, int64_t> groupOccurences;
      while(file_ai >> rp >> rn >> rg){
        nbfiles++;
      }
      file_ai.close();
      cerr << "Debug : nbfiles = " << nbfiles << endl;
    }
    else{
      exit(4651);
    }
    cerr << "Debug : create object" << endl;
    ProfileOfKmers result = ProfileOfKmers(nbfiles, 30, 0);
    result.addListOfFiles(plist);
    cerr << "Debug : launching discretization" << endl;
    result.discretizationWithR(pathInpFile, pathResFile, sizeBuff, true);
}



void getReductionWithUncertainty(string plist, string pathInpFile, string pathResFile, long double threshold){
  ifstream file_ai(plist.c_str(), ios::in);
  int64_t nbfiles = 0;
  string rp, rn, rg;
  if(file_ai){
    map<string, int64_t> groupOccurences;
    while(file_ai >> rp >> rn >> rg){
      nbfiles++;
    }
    file_ai.close();
    cerr << "Debug : nbfiles = " << nbfiles << endl;
  }
  else{
    exit(4651);
  }
  ProfileOfKmers result = ProfileOfKmers(nbfiles, 30, 0);
  result.RemoveNonInformativeKmers(pathInpFile, pathResFile, threshold, false);
}

int main(int argc, char* argv[]){
    //fixthresholds();
    //string pathListFile = argv[1]; // list of path
    //string pathResFile = argv[2]; // output
    //makeBigMatrix(pathListFile,pathResFile);
    //testEntropyDecomposition();

    map<string, string> parameters;
    string mode(argv[1]);
    parameters["mode"] = mode;
    parameters["threshold"] = "0";
    parameters["sample"] = "1000000";
    parameters["buffer"] = "100000";

    for(int64_t kp = 2;kp<argc; kp+=2){
      string paramc = argv[kp];
      if(paramc.compare("-i")==0){
        parameters["input"] = argv[kp+1];
      }
      if(paramc.compare("-o")==0){
        parameters["output"] = argv[kp+1];
      }
      if(paramc.compare("-groups")==0){
        parameters["groups"] = argv[kp+1];
      }
      if(paramc.compare("-threshold")==0){
        parameters["threshold"] = argv[kp+1];
      }
      if(paramc.compare("-sample")==0){
        parameters["sample"] = argv[kp+1];
      }
      if(paramc.compare("-buffer")==0){
        parameters["buffer"] = argv[kp+1];
      }

    }

    cout << parameters["mode"] << "\t" << parameters["input"] << "\t" << parameters["output"] << endl;
    if(parameters["mode"].compare("importation")==0){
      //test
      makeBigMatrix(parameters["input"],parameters["output"],stoi(parameters["threshold"]));
    }
    if(parameters["mode"].compare("discretization")==0){
      getDiscretizationMDL(parameters["groups"], parameters["input"], parameters["output"], stoi(parameters["buffer"]));
    }
    if(parameters["mode"].compare("discretizationAmeva")==0){
      getDiscretizationAmeva(parameters["groups"], parameters["input"], parameters["output"], stoi(parameters["buffer"]));
    }
    if(parameters["mode"].compare("discretizationR")==0){
      getDiscretizationR(parameters["groups"], parameters["input"], parameters["output"], stoi(parameters["buffer"]));
    }
    if(parameters["mode"].compare("reduction")==0){
      getReductionWithUncertainty(parameters["groups"], parameters["input"], parameters["output"], 0.5);
    }
    if(parameters["mode"].compare("threshold")==0){
      fixthresholds(parameters["input"], stoi(parameters["sample"]));
    }


    return 0;
}
