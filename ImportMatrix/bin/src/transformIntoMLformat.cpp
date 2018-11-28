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
/*
 *translate kmer matrix after extraction filter task to be read by the genetic algorithm
 * Arguments:
 * - path to the input file
 * - path for the output file
 * - maximum memory occupation for the matrix during the process in Go ( Default : 32 )
 */

int main(int argc, char* argv[]){
    string inputfile(argv[1]);
    string outputfile(argv[2]);

    float maxsizesubmat=32;
    if (argc>=4)
        maxsizesubmat=stof(argv[3]);
    vector<string> _Groups;
    vector<string> _Kmerid;
    vector<vector<float>> _Data;
    vector<ofstream*> tmp_out_file;
    
    
    bool breakmaxsize=false;
    bool firstwriting=true;
    bool gotline=true;
    uint linecount=0;

    string line;
    /*!
    * OPEN FILES the file
    * */
    ifstream in_file(inputfile, ios::in);
    if (!in_file) {
        cerr << "Could not open input file\n";
        return EXIT_FAILURE;
    }
    
    /*!
    * Reading the file
    * */


    getline(in_file, line);
    getline(in_file, line);
    boost::split(_Groups, line, boost::is_any_of("\t"));

    
    int iter=0;
    do{
        //read inputfile until linelimit of end of file
        breakmaxsize=false;
        
        do{
            if (getline(in_file, line)){

                vector<string> content;
                boost::split(content, line, boost::is_any_of("\t"));
                _Kmerid.push_back(content[0]);
                vector<float> tmpdata;
                for(uint igroup = 1; igroup<content.size(); igroup++){
                    tmpdata.push_back(stof(content[igroup],NULL));
                }
                _Data.push_back(tmpdata);
                linecount++;
                if (_Data.size()*_Groups.size()*sizeof(float)>=maxsizesubmat*800000000){
                    breakmaxsize=true;
                    
                }
                gotline=true;
            }else{
                gotline=false;
            }

            
        }while(gotline & (!breakmaxsize));
        
         /*!
        * Output
        * */
        cout<<"Write the "<<linecount<< " first k-mers of the matrix"<<endl;
        iter++;
        
        string str;
        long pos =0;

        for(uint64_t i = 0; i < _Groups.size() ; i++){
            if( firstwriting){
                tmp_out_file.push_back(new ofstream(outputfile+to_string(i)+".tmp", ios::app));
            //fstream tmp_out_file[i].open(outputfile+to_string(i)+".tmp", ios::app);//| fstream::app
            }
            if (!(*tmp_out_file[i])) {
                cerr << "Could not create output file\n";
                return EXIT_FAILURE;
            }
            if( firstwriting){
                (*tmp_out_file[i]) << _Groups[i];//<<" fingroup "
            }
           
            for(uint64_t j = 0; j < _Kmerid.size(); j++){
                if(i==0){
                    (*tmp_out_file[i]) << "," << _Kmerid[j];
                }
                else{
                    (*tmp_out_file[i]) << "," << _Data[j][i-1];
                }
            }
            //
        }
        
        _Data.clear();        
        firstwriting=false;
        _Kmerid.clear();




    }while(gotline );    
   
   for(uint64_t i = 0; i < _Groups.size() ; i++){
       tmp_out_file[i]->close();
   
   }
   in_file.close();
   
   //concatenate file
    cout << "Concatenate the final matrix file"<<endl;
    ofstream final_out_file(outputfile, ios::out );//| fstream::app
    if (!final_out_file) {
        cerr << "Could not create output file\n";
        return EXIT_FAILURE;
    }
   
        for(uint64_t i = 0; i < _Groups.size() ; i++){
            ifstream tmp_in_file(outputfile+to_string(i)+".tmp", ios::in);//| fstream::app
            getline(tmp_in_file, line);
            final_out_file<<line<<endl;
            tmp_in_file.close();
            if( remove( (outputfile+to_string(i)+".tmp").c_str() ) != 0 )
                perror( "Error deleting file" );
        }
    
    
    
    final_out_file.close();
/*
   
    for(uint64_t i = 0; i < _Groups.size() ; i++){
        cout << _Groups[i];
        for(uint64_t j = 0; j < _Kmerid.size(); j++){
            if(i==0){
                cout << "," << _Kmerid[j];
            }
            else{
                cout << "," << _Data[j][i-1];
            }
        }
        cout << endl;
    }
*/

//const char * add_string = argv[3];
//
//
//       for(uint64_t i = 0; i < _Groups.size() ; i++){
//            fstream out_file(outputfile+str(i)+."tmp", ios::app);//| fstream::app
//            
//            if( firstwriting){
//                out_file << _Groups[i];//<<" fingroup "
//            }else{
//                //copy outpuposition into inputposition
////                long pos = out_file.tellp();
//                out_file.seekg ( pos);
//                getline(out_file, str);
//                //copy input position into output position ,go back one for \n
//                pos = out_file.tellp();
//                out_file.seekp ( pos-1);
//
//                cout<<"second writing...reading ["<<str<<"]"<<endl;
//            }
//           
//            for(uint64_t j = 0; j < _Kmerid.size(); j++){
//                if(i==0){
//                    out_file << "," << _Kmerid[j];//<< "kmeridj="<<j
//                    //out_file.write( char("," + _Kmerid[j]));
//                }
//                else{
//                    out_file << "," << _Data[j][i-1];//<< "kmeridj="<<j
//                    //out_file.write( char("," + _Data[j][i-1]))
//                }
//                
//            }
//            if( firstwriting){
//                out_file <<endl;
//            }
//        }
//        firstwriting=false;
//        _Kmerid.clear();
//        _Data.clear();
//
//    }while(gotline );    
//   
//   out_file.close();
//   in_file.close();



    return EXIT_SUCCESS;
}