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
//#include <omp.h>


using namespace std;

//!
//! \file   IndexingKMers.cpp
//! \brief  The purpose is to handle a list of kmer tabulated files (ASCII)
//!         to create a matrix of the most explicative kmers
//! \author Aubin.Thomas William.Ritchie
//! \version    1.0
//! \date   may 2017



//!
//! \brief Naive class to deal with discretization of continuous variables according supervised MDL algorithm
//! author : Aubin Thomas
//!

class EntropyDecomposition {
public:
    uint64_t _Nmax;
    vector <long double> _Tthresholds;

    EntropyDecomposition(int64_t pNmax = 100) {
        _Nmax = pNmax;
    }

    /** Discretization
    Discretization : apply MDL discretization on vector of values according vector of classes
     **/
    vector<string> Discretization(vector<long double> pvector_values, vector<string> pvector_classes, bool pdebug) {
        //sort vectors
        vector<long double> pvector_values_cp(pvector_values); // copy that will be deleted along treatment
        vector<long double> pvector_values_ord(pvector_values); // sorted according values
        vector<string> pvector_classes_cp(pvector_classes);
        vector<string> pvector_classes_ord(pvector_classes); // sorted according vector of values (by association)
        map<string, long double> listClasses; // counts of classes
        uint32_t vsize = pvector_values.size();
        vector<string> result(vsize, "");
        //sort values
        sort(pvector_values_ord.begin(), pvector_values_ord.end());

        //get population of classes and sort classes according sort of values
        for (uint32_t i = 0; i < vsize; i++) {
            uint32_t vsize2 = pvector_values_cp.size(); //size that will decrease : vsize to 0
            //update of counts of each classes
            if (listClasses.find(pvector_classes[i]) != listClasses.end()) {
                listClasses[pvector_classes[i]] += 1;
            } else {
                listClasses[pvector_classes[i]] = 1;
            }

            //search of corresponding class
            for (uint32_t j = 0; j < vsize2; j++) {
                if (pvector_values_cp[j] == pvector_values_ord[i]) {
                    pvector_classes_ord[i] = pvector_classes_cp[j];
                    pvector_classes_cp.erase(pvector_classes_cp.begin() + j);
                    pvector_values_cp.erase(pvector_values_cp.begin() + j);
                    j = vsize2;
                }
            }
        }
        if (pdebug) {
            cerr << "initial values : ";
            for (uint32_t i = 0; i < vsize; i++) {
                cerr << " " << pvector_values_ord[i];
            }
            cerr << endl;
            cerr << "initial classes : ";
            for (uint32_t i = 0; i < vsize; i++) {
                cerr << " " << pvector_classes_ord[i];
            }
            cerr << endl;
        }

        //call to function to divide
        DivideIntoBinary(pvector_values_ord, pvector_classes_ord, listClasses, pdebug);

        //then give result in function of thresholds
        for (uint32_t i = 0; i < vsize; i++) {
            long double value_c = pvector_values[i];
            unsigned long long rank_c = 0;
            for (int64_t i = 0; i < _Tthresholds.size(); i++) {
                if (value_c >= _Tthresholds[i]) {
                    rank_c += 1;
                }
            }
            result[i] = to_string(rank_c);
        }
        if (pdebug) {
            for (uint32_t i = 0; i < vsize; i++) {
                cerr << result[i] << " " << endl;
            }
        }
        return result;
    }

    /** DivideIntoBinary
    DivideIntoBinary : iterative function to search for thresholds that maximize
     **/
    void DivideIntoBinary(vector<long double> pvector_values, vector<string> pvector_classes, map<string, long double> pmapclasses, bool pdebug) {
        if (pdebug) {
            cerr << endl << "DB values : ";
            for (uint32_t i = 0; i < pvector_values.size(); i++) {
                cerr << " " << pvector_values[i];
            }
            cerr << endl;
            cerr << "DB classes : ";
            for (uint32_t i = 0; i < pvector_values.size(); i++) {
                cerr << " " << pvector_classes[i];
            }
            cerr << endl;
        }

        //entropy0 : the initial entropy according classes
        long double entropy0 = 0;
        long double taille = pvector_values.size();
        for (map<string, long double>::iterator it = pmapclasses.begin(); it != pmapclasses.end(); ++it) {
            string class_c = it->first;
            long double number_c = it->second;
            long double ratio0 = (number_c / taille);
            entropy0 -= ratio0 * log2(ratio0);
        }
        if (pdebug) cerr << "\tentropy0 = " << entropy0 << endl;

        //search of threshold TA that maximize the gain of entropy
        long double GainMax = 0;
        uint32_t kmax = 0;
        long double entropyTA_up_Max;
        long double entropyTA_down_Max;

        long double TAprev = ((pvector_values[0] + pvector_values[1]) / 2) - 1;
        for (uint32_t k = 0; k < (pvector_values.size() - 1); k++) {
            long double TA = (pvector_values[k] + pvector_values[k + 1]) / 2;
            if (TA != TAprev) {

                long double entropyTA_up = 0; // entropy for events of A >= TA
                long double entropyTA_down = 0; // entropy for events of A < TA
                if (pdebug) cerr << "\tk = " << k << "\tTA = " << TA << endl;

                long double nb_A_up = 0;
                long double nb_A_down = 0;
                //computing the number of values >=TA and <TA
                for (uint32_t i = 0; i < pvector_values.size(); i++) {
                    if (pvector_values[i] >= TA) {
                        nb_A_up++;
                    } else {
                        nb_A_down++;
                    }
                }

                //computing entropies for >=TA and <TA by adding each part of each class
                for (map<string, long double>::iterator it = pmapclasses.begin(); it != pmapclasses.end(); ++it) {
                    string class_c = it->first;
                    long double nb_A_up_class_c = 0;
                    long double nb_A_down_class_c = 0;

                    for (uint32_t i = 0; i < pvector_values.size(); i++) {
                        if (pvector_classes[i].compare(class_c) == 0) {

                            if (pvector_values[i] >= TA) {
                                nb_A_up_class_c++;
                            } else {
                                nb_A_down_class_c++;
                            }
                        }
                    }
                    //if(pdebug) cerr << "\t\tclass_c = " << class_c << "\tnb_A_up_class_c = " << nb_A_up_class_c << "\tnb_A_down_class_c = " << nb_A_down_class_c << endl;


                    long double ratio_up = nb_A_up_class_c / nb_A_up;
                    long double ratio_down = nb_A_down_class_c / nb_A_down;
                    if ((ratio_up != 0) and (nb_A_up != 0)) entropyTA_up += -ratio_up * log2(ratio_up);
                    if ((ratio_down != 0) and (nb_A_down != 0)) entropyTA_down += -ratio_down * log2(ratio_down);

                }

                if (pdebug) cerr << "\t\tentropyTA_up = " << entropyTA_up << endl << "\t\tentropyTA_down = " << entropyTA_down << endl;
                //computing gain as entropy0-(ratio of up plus down)
                long double gain_c = entropy0 - ((nb_A_up / taille) * entropyTA_up)-((nb_A_down / taille) * entropyTA_down);
                if (pdebug) cerr << "\t\tgain_c = " << entropy0 << " - " << ((nb_A_up / taille) * entropyTA_up)+((nb_A_down / taille) * entropyTA_down) << " =  " << gain_c << endl;
                if (gain_c > GainMax) {
                    GainMax = gain_c;
                    kmax = k;
                    entropyTA_up_Max = entropyTA_up;
                    entropyTA_down_Max = entropyTA_down;
                }
                TAprev = TA;
            }
        }
        if (pdebug) cerr << "\tkmax = " << kmax << "\tentropyTA_up_Max = " << entropyTA_up_Max << "\tentropyTA_down_Max = " << entropyTA_down_Max << endl;

        //update pmapclasses for each subset
        map<string, long double> listClassesUp;
        map<string, long double> listClassesDown;
        long double TA = (pvector_values[kmax] + pvector_values[kmax + 1]) / 2;
        for (uint32_t i = 0; i < pvector_values.size(); i++) {
            if (pvector_values[i] >= TA) {
                if (listClassesUp.find(pvector_classes[i]) != listClassesUp.end()) {
                    listClassesUp[pvector_classes[i]] += 1;
                } else {
                    listClassesUp[pvector_classes[i]] = 1;
                }
            } else {
                if (listClassesDown.find(pvector_classes[i]) != listClassesDown.end()) {
                    listClassesDown[pvector_classes[i]] += 1;
                } else {
                    listClassesDown[pvector_classes[i]] = 1;
                }
            }
        }


        long double GainThreshold = (log2(taille - 1) / taille) + (log2(pow(3, pmapclasses.size()) - 2) - pmapclasses.size() * entropy0 + listClassesUp.size() * entropyTA_up_Max + listClassesDown.size() * entropyTA_down_Max) / taille;
        if ((taille <= 1) or (GainThreshold < 0)) GainThreshold = 1;

        if (pdebug) {
            cerr << "taille = " << taille << endl;
            cerr << "pmapclasses.size() = " << pmapclasses.size() << endl;
            cerr << "entropy0 = " << entropy0 << endl;
            cerr << "listClassesUp.size() = " << listClassesUp.size() << endl;
            cerr << "listClassesDown.size() = " << listClassesDown.size() << endl;
            cerr << "entropyTA_up_Max = " << entropyTA_up_Max << endl;
            cerr << "entropyTA_down_Max = " << entropyTA_down_Max << endl;
            cerr << "\tGainThreshold = " << GainThreshold << endl;
        }

        if ((GainMax > GainThreshold) and (_Tthresholds.size() <= _Nmax)) {
            _Tthresholds.push_back(TA);
            //decoupe vecteur selon TA
            vector<long double> pvector_values_Up;
            vector<string> pvector_classes_Up;
            vector<long double> pvector_values_Down;
            vector<string> pvector_classes_Down;

            for (uint32_t i = 0; i < pvector_values.size(); i++) {
                if (pvector_values[i] >= TA) {
                    pvector_values_Up.push_back(pvector_values[i]);
                    pvector_classes_Up.push_back(pvector_classes[i]);
                } else {
                    pvector_values_Down.push_back(pvector_values[i]);
                    pvector_classes_Down.push_back(pvector_classes[i]);
                }
            }

            DivideIntoBinary(pvector_values_Down, pvector_classes_Down, listClassesDown, pdebug);
            DivideIntoBinary(pvector_values_Up, pvector_classes_Up, listClassesUp, pdebug);
        }

    }

    /** what is entropy for a given value? **/
    long double NodeBinaryEntropy(vector<long double> pvec, long double pthreshold) {
        long double totalNumberOfValues = pvec.size();
        long double totalUp = 0;
        long double totalDown = 0;
        for (uint32_t i = 0; i < totalNumberOfValues; i++) {
            if (pvec[i] >= pthreshold) {
                totalUp += 1;
            }
        }
        totalDown = totalNumberOfValues - totalUp;
        long double fracUp = totalUp / totalNumberOfValues;
        long double fracDown = totalDown / totalNumberOfValues;
        return -fracUp * log2(fracUp) - fracDown * log2(fracDown);
    }

    /** SymmetricUncertainty between 2 vectors of equal size**/
    long double _old_flexible_SymmetricUncertainty(vector<string> pvectorA, vector<string> pvectorB, bool pdebug) {
        long double result = 0;
        long double HA = 0;
        long double HB = 0;
        long double HAB = 0;

        map<string, long double> mapA;
        map<string, long double> mapB;
        map<string, long double> mapAB;

        for (uint64_t i = 0; i < pvectorA.size(); i++) {
            if (mapA.find(pvectorA[i]) != mapA.end()) {
                mapA[pvectorA[i]] += 1;
            } else {
                mapA[pvectorA[i]] = 1;
            }
            if (mapB.find(pvectorB[i]) != mapB.end()) {
                mapB[pvectorB[i]] += 1;
            } else {
                mapB[pvectorB[i]] = 1;
            }
        }
        for (uint64_t i = 0; i < pvectorA.size(); i++) {

            string keyc = pvectorA[i] + "_" + pvectorB[i];
            if (mapAB.find(keyc) != mapAB.end()) {
                mapAB[keyc] += 1;
            } else {
                mapAB[keyc] = 1;
            }

        }
        if (pdebug)
            cerr << endl << "HA : " << endl;
        for (map<string, long double>::iterator ia = mapA.begin(); ia != mapA.end(); ia++) {
            long double ratioA = (ia->second) / (pvectorA.size());
            if (pdebug)
                cerr << ia->first << "\t" << ia->second << " / " << pvectorA.size() << " = " << ratioA << endl;
            if ((ratioA != 0) and (pvectorA.size() > 0))
                HA -= (ratioA * log2(ratioA));
            if (pdebug)
                cerr << "\tHA = " << HA << endl;
        }

        if (pdebug)
            cerr << endl << "HB : " << endl;
        for (map<string, long double>::iterator ib = mapB.begin(); ib != mapB.end(); ib++) {
            long double ratioB = (ib->second) / (pvectorB.size());
            if (pdebug)
                cerr << ib->first << "\t" << ib->second << " / " << pvectorB.size() << " = " << ratioB << endl;
            if ((ratioB != 0) and (pvectorB.size() > 0))
                HB -= (ratioB * log2(ratioB));
            if (pdebug)
                cerr << "\tHB = " << HB << endl;
        }

        if (pdebug)
            cerr << endl << "HAB : " << endl;
        for (map<string, long double>::iterator iab = mapAB.begin(); iab != mapAB.end(); iab++) {
            long double ratioAB = (iab->second) / (pvectorB.size());
            if (pdebug)
                cerr << iab->first << "\t" << iab->second << " / " << pvectorB.size() << " = " << ratioAB << endl;
            if ((ratioAB != 0) and (pvectorB.size() > 0))
                HAB -= (ratioAB * log2(ratioAB));
            if (pdebug)
                cerr << "\tHAB = " << HAB << endl;
        }
        result = 2 * (HA + HB - HAB) / (HA + HB);
        return result;

    }
    
    
    
    /** SymmetricUncertainty between 2 vectors of equal size**/
    long double _old_flexible_SymmetricUncertaintyWithInt(vector<int> pvectorA, vector<int> pvectorB, bool pdebug) {
        long double result = 0;
        long double HA = 0;
        long double HB = 0;
        long double HAB = 0;

        map<string, long double> mapA;
        map<string, long double> mapB;
        map<string, long double> mapAB;

        for (uint64_t i = 0; i < pvectorA.size(); i++) {
            string valA = to_string(pvectorA[i]);
            string valB = to_string(pvectorB[i]);
            string keyc = valA + "_" + valB;
            if (mapA.find(valA) != mapA.end()) {
                mapA[valA] += 1;
            } else {
                mapA[valA] = 1;
            }
            if (mapB.find(valB) != mapB.end()) {
                mapB[valB] += 1;
            } else {
                mapB[valB] = 1;
            }
        
            if (mapAB.find(keyc) != mapAB.end()) {
                mapAB[keyc] += 1;
            } else {
                mapAB[keyc] = 1;
            }

        }
        if (pdebug)
            cerr << endl << "HA : " << endl;
        for (map<string, long double>::iterator ia = mapA.begin(); ia != mapA.end(); ia++) {
            long double ratioA = (ia->second) / (pvectorA.size());
            if (pdebug)
                cerr << ia->first << "\t" << ia->second << " / " << pvectorA.size() << " = " << ratioA << endl;
            if ((ratioA != 0) and (pvectorA.size() > 0))
                HA -= (ratioA * log2(ratioA));
            if (pdebug)
                cerr << "\tHA = " << HA << endl;
        }

        if (pdebug)
            cerr << endl << "HB : " << endl;
        for (map<string, long double>::iterator ib = mapB.begin(); ib != mapB.end(); ib++) {
            long double ratioB = (ib->second) / (pvectorB.size());
            if (pdebug)
                cerr << ib->first << "\t" << ib->second << " / " << pvectorB.size() << " = " << ratioB << endl;
            if ((ratioB != 0) and (pvectorB.size() > 0))
                HB -= (ratioB * log2(ratioB));
            if (pdebug)
                cerr << "\tHB = " << HB << endl;
        }

        if (pdebug)
            cerr << endl << "HAB : " << endl;
        for (map<string, long double>::iterator iab = mapAB.begin(); iab != mapAB.end(); iab++) {
            long double ratioAB = (iab->second) / (pvectorB.size());
            if (pdebug)
                cerr << iab->first << "\t" << iab->second << " / " << pvectorB.size() << " = " << ratioAB << endl;
            if ((ratioAB != 0) and (pvectorB.size() > 0))
                HAB -= (ratioAB * log2(ratioAB));
            if (pdebug)
                cerr << "\tHAB = " << HAB << endl;
        }
        result = 2 * (HA + HB - HAB) / (HA + HB);
        return result;

    }

    long double SymmetricUncertainty(vector<string> pvectorA, vector<string> pvectorB, bool pdebug) {
        const short mlength = 29; // number of state max in the vector for 29 : the differents states have to be a string of number between 0 and 28
        long double result = 0;
        long double HA = 0;
        long double HB = 0;
        long double HAB = 0;
        long double ratioA, ratioB, ratioAB;

        int ia, ib;
        vector<long double> mapA(mlength, 0);
        vector<long double> mapB(mlength, 0);
        vector<long double> mapAB((mlength * mlength) + mlength, 0);
        if (pvectorA.size() != 0) {
            for (uint64_t i = 0; i < pvectorA.size(); i++) {
                ia = stoi(pvectorA[i]);
                ib = stoi(pvectorB[i]);
                mapA[ia]++;
                mapB[ib]++;
                mapAB[ia + (mlength * ib)]++;
            }

            for (uint64_t i = 0; i < mlength; i++) {
                if (mapA[i] != 0) {
                    ratioA = mapA[i] / pvectorA.size();
                    //cerr << "\nmapA[i] = " << mapA[i] << "\tpvectorA.size() = " << pvectorA.size()<<endl;
                    HA -= (ratioA * log2(ratioA));
                    //cerr << "\nHA = " << HA << "\tratio = " << ratioA<<endl;
                }
                if (mapB[i] != 0) {
                    ratioB = mapB[i] / pvectorA.size();
                    HB -= (ratioB * log2(ratioB));

                }
            }

            for (uint64_t i = 0; i < (mlength * mlength) + mlength; i++) {
                if (mapAB[i] != 0) {
                    ratioAB = mapAB[i] / pvectorA.size();
                    HAB -= (ratioAB * log2(ratioAB));
                }

            }

            result = 2 * (HA + HB - HAB) / (HA + HB);
            return result;
        } else {
            return 0;
        }

    }

    long double SymmetricUncertaintyWithInt(vector<int> pvectorA, vector<int> pvectorB, bool pdebug) {
        const short mlength = 100; // number of state max in the vector for 100 : the differents states have to be a string of number between 0 and 99
        long double result = 0;
        long double HA = 0;
        long double HB = 0;
        long double HAB = 0;
        long double ratioA, ratioB, ratioAB;

        int ia, ib;
        vector<long double> mapA(mlength, 0);
        vector<long double> mapB(mlength, 0);
        vector<long double> mapAB(mlength*mlength+2*mlength, 0);
        if (pvectorA.size() != 0) {
            for (uint64_t i = 0; i < pvectorA.size(); i++) {
                ia = pvectorA[i];
                ib = pvectorB[i];
                mapA[ia]++;
                mapB[ib]++;
                mapAB[ia+ib*mlength+mlength]++;
            }

            for (uint64_t i = 0; i < mlength; i++) {
                if (mapA[i] != 0) {
                    ratioA = mapA[i] / pvectorA.size();
                    //cerr << "\nmapA[i] = " << mapA[i] << "\tpvectorA.size() = " << pvectorA.size()<<endl;
                    HA -= (ratioA * log2(ratioA));
                    //cerr << "\nHA = " << HA << "\tratio = " << ratioA<<endl;
                }
                if (mapB[i] != 0) {
                    ratioB = mapB[i] / pvectorA.size();
                    //cerr << "\nmapB[i] = " << mapB[i] << "\tpvectorB.size() = " << pvectorB.size()<<endl;
                    HB -= (ratioB * log2(ratioB));
                    //cerr << "\nHB = " << HB << "\tratio = " << ratioB<<endl;

                }
            }

            for (uint64_t i = 0; i < (mlength * mlength) + 2*mlength; i++) {
                if (mapAB[i] != 0) {
                    ratioAB = mapAB[i] / pvectorA.size();
                    //cerr << "\nmapAB[i] = " << mapAB[i] << "\tpvectorB.size() = " << pvectorB.size()<<endl;
                    HAB -= (ratioAB * log2(ratioAB));
                    //cerr << "\nHB = " << HB << "\tratio = " << ratioB<<endl;
                }

            }

            result = 2 * (HA + HB - HAB) / (HA + HB);
            return result;
        } else {
            return 0;
        }

    }



	long double SymmetricUncertaintyWithIntLimited(vector<int> pvectorA, vector<int> pvectorB, int plimitation, bool pdebug) {
		const short mlength = 25; // number of state max in the vector for 100 : the differents states have to be a string of number between 0 and 99
		long double result = 0;
		long double HA = 0;
		long double HB = 0;
		long double HAB = 0;
		long double ratioA, ratioB, ratioAB;

		int ia, ib;
		vector<long double> mapA(mlength, 0);
		vector<long double> mapB(mlength, 0);
		vector<long double> mapAB(mlength*mlength + 2 * mlength, 0);
		if (pvectorA.size() != 0) {
			for (uint64_t i = 0; i < plimitation; i++) {
				ia = pvectorA[i];
				ib = pvectorB[i];
				mapA[ia]++;
				mapB[ib]++;
				mapAB[ia + ib * mlength + mlength]++;
			}

			for (uint64_t i = 0; i < mlength; i++) {
				if (mapA[i] != 0) {
					ratioA = mapA[i] / plimitation;
					//cerr << "\nmapA[i] = " << mapA[i] << "\tpvectorA.size() = " << pvectorA.size()<<endl;
					HA -= (ratioA * log2(ratioA));
					//cerr << "\nHA = " << HA << "\tratio = " << ratioA<<endl;
				}
				if (mapB[i] != 0) {
					ratioB = mapB[i] / plimitation;
					//cerr << "\nmapB[i] = " << mapB[i] << "\tpvectorB.size() = " << pvectorB.size()<<endl;
					HB -= (ratioB * log2(ratioB));
					//cerr << "\nHB = " << HB << "\tratio = " << ratioB<<endl;

				}
			}

			for (uint64_t i = 0; i < (mlength * mlength) + 2 * mlength; i++) {
				if (mapAB[i] != 0) {
					ratioAB = mapAB[i] / plimitation;
					//cerr << "\nmapAB[i] = " << mapAB[i] << "\tpvectorB.size() = " << pvectorB.size()<<endl;
					HAB -= (ratioAB * log2(ratioAB));
					//cerr << "\nHB = " << HB << "\tratio = " << ratioB<<endl;
				}

			}

			result = 2 * (HA + HB - HAB) / (HA + HB);
			return result;
		}
		else {
			return 0;
		}

	}

}; //class1



/*

int main(){
    vector<string> A;
    vector<string> B;
    EntropyDecomposition test = EntropyDecomposition();
    A.push_back("1");A.push_back("1");A.push_back("2");A.push_back("3");A.push_back("2");A.push_back("2");A.push_back("3");A.push_back("2");A.push_back("1");A.push_back("1");A.push_back("2");A.push_back("1");
    A.push_back("1");A.push_back("1");A.push_back("2");A.push_back("3");A.push_back("2");A.push_back("2");A.push_back("3");A.push_back("2");A.push_back("1");A.push_back("1");A.push_back("2");A.push_back("1");
    A.push_back("1");A.push_back("1");A.push_back("2");A.push_back("3");A.push_back("2");A.push_back("2");A.push_back("3");A.push_back("2");A.push_back("1");A.push_back("1");A.push_back("2");A.push_back("1");
    B.push_back("3");B.push_back("3");B.push_back("2");B.push_back("1");B.push_back("4");B.push_back("4");B.push_back("1");B.push_back("2");B.push_back("3");B.push_back("3");B.push_back("2");B.push_back("3");
    B.push_back("3");B.push_back("3");B.push_back("2");B.push_back("1");B.push_back("2");B.push_back("2");B.push_back("1");B.push_back("2");B.push_back("3");B.push_back("3");B.push_back("2");B.push_back("3");
    B.push_back("3");B.push_back("3");B.push_back("2");B.push_back("1");B.push_back("2");B.push_back("2");B.push_back("1");B.push_back("2");B.push_back("3");B.push_back("3");B.push_back("2");B.push_back("3");
    cout << test.SymmetricUncertainty(A,B,false) << endl;

    return EXIT_SUCCESS;
}
*/



