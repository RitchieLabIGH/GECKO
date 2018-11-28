/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "ClassPopulation.cpp"
#include <sys/stat.h>
#include <unistd.h>

int main( int argc, char *argv[] ){
    //phidneuron = 10, pgen = 70, pind = 20, pkmer = 50

    
    string file_csv = "../BEAUTY_RNASeq_NORM_discret_AMEVA_OnlyInformative_HammingMI_TripleNegative.csv";
    //string path2occs = "../BEAUTYtrinegRes/figBEAUTY_trineg9_k25_mode0/BestIndiv20.csvforextractkm.count";
    string path2occs = "../BEAUTYtrinegRes/fig/countkmerall0.txt";
    //string file_csv = "../faketabcom.txt";
    //string path2occs = "../kmerlist.txt";
    int minocc =0;
    vector <string> vectPath2occs;
 
    if (argc > 1){
        file_csv  = argv[1] ; //vector comma separated
    }
    if (argc > 2){
        path2occs= argv[2];

    }
    if (argc > 3){
        minocc= stoi(argv[3]);
    }
    

    ClassPopulation popA = ClassPopulation();
    //popA.loadCSVwithConstraint(file_csv, path2occs, 0, true);
    popA.loadCSV(file_csv,  true);
    popA.loadIndivName(file_csv);
    boost::split(vectPath2occs, path2occs, boost::is_any_of(","));
    cout<< "Number of extraction to do "<<vectPath2occs.size()<<endl;
    for(unsigned int i=0;i<vectPath2occs.size();i++){
        //cout<< "Extracting  "<<vectPath2occs[i]<<"_minocc"<<minocc<<"_SampleMat.csv"<<endl;
        struct stat buffer;
        if (stat(vectPath2occs[i].c_str(), &buffer) == 0)
            if(minocc>1){
                cout<< "Extracting  "<<vectPath2occs[i]<<"_minocc"<<minocc<<"_SampleMat.csv"<<endl;

                popA.printDatas(vectPath2occs[i],vectPath2occs[i]+"_minocc"+std::to_string(minocc)+"_SampleMat.csv",minocc);
            }
            else{
                cout<< "Extracting  "<<vectPath2occs[i]<<"_SampleMat.csv"<<endl;

                popA.printDatas(vectPath2occs[i],vectPath2occs[i]+"_SampleMat.csv",minocc);
            }
        else
            cout<< "No file found  "<<vectPath2occs[i]<<endl;
    }



    return EXIT_SUCCESS;
}
