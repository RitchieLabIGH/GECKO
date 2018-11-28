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

using namespace std;


struct SeqAndCount {
	string sequence;
	uint64_t count;
};


bool compareByCount(const SeqAndCount &a, const SeqAndCount &b) {
	return a.count > b.count;
}




class ClassGoodIndividuals {
public:

	string _self_file; //the csv file aka goodCandidates.csv
	uint64_t _self_sizeKmer;

	ClassGoodIndividuals(string p_file, uint64_t p_sizeKmer) {
		_self_file = p_file;
		_self_sizeKmer = p_sizeKmer;
	}


	/********************************************************************************************
	* Helpers for computations
	* Computing score, sequence reformats, get couples of index,
	********************************************************************************************/
	/*!
	* \brief get all the scores from test and outer scores (see description of importCSVfile function)
	* p_score can be :
	*           - test
	*           - outer
	*           - prod : test*outer
	*           - adjust : test*outer-overfit
	* used in importCSVfile
	* */
	float getScores(float p_test, float p_outer, string p_score) {
		float overfit = p_test - p_outer;
		if (overfit < 0) {
			overfit = 0;
		}
		if (p_score.compare("test") == 0)
			return p_test;
		if (p_score.compare("outer") == 0)
			return p_outer;
		if (p_score.compare("pro") == 0)
			return p_test * p_outer;
		if (p_score.compare("adjust") == 0)
			return p_test * p_outer - overfit;
		return 0.0;
	}

	/*!
	* \brief reformat a list of IDs in ascending order
	* used in importCSVfile
	* */
	string reformatSequences(string p_seq) {
		vector<string> contentSeq;
		boost::split(contentSeq, p_seq, boost::is_any_of(","));
		sort(contentSeq.begin(), contentSeq.end());
		string result = contentSeq[0];
		for (uint i = 1; i < contentSeq.size(); i++) {
			result += "," + contentSeq[i];
		}
		return result;
	}

	/*!
	* \brief get all couples of sequences
	* used in computeNetwork
	* */
	vector<vector<string>> getCouplesOfIndex(string p_seq) {
		//ptrdiff_t pos = distance(Names.begin(), find(Names.begin(), Names.end(), old_name_));

		vector<vector<string>> result;
		vector<string> contentSeq;
		boost::split(contentSeq, p_seq, boost::is_any_of(","));

		for (uint i = 0; i < contentSeq.size(); i++) {
			string id1 = contentSeq[i];
			for (uint j = i + 1; j < contentSeq.size(); j++) {

				string id2 = contentSeq[j];
				vector<string> tmp;
				tmp.push_back(id1);
				tmp.push_back(id2);
				result.push_back(tmp);
			}
		}
		return result;
	}

	//! \fn uint64_t getIndice(string kmer, int64_t rank=0)
	//! \brief transform a sequence into an integer value that corresponds to {A,C,T,G} in base-4 algebra
	//! \return (uint64_t) indice that corresponds to sequence
	//!
	uint64_t getIndice(string kmer, int64_t rank = 0) {
		if (kmer.size() == rank) {
			return 0;
		}
		else {
			char character_c = kmer[kmer.size() - rank - 1];
			if (character_c == 'A') {
				return getIndice(kmer, rank + 1) + 0;
			}
			else {
				if (character_c == 'C') {
					return getIndice(kmer, rank + 1) + 1 * powl(4, rank);
				}
				else {
					if (character_c == 'T') {
						return getIndice(kmer, rank + 1) + 2 * powl(4, rank);
					}
					else {
						return getIndice(kmer, rank + 1) + 3 * powl(4, rank);
					}
				}
			}
		}
	}

	//! transform sequence ID to real nucleotidic sequence
	string idseq2fasta(uint64_t realid, uint64_t nbBasePerKmer) {
		string result = "";
		int64_t k = nbBasePerKmer - 1;
		vector<string> letters;
		letters.push_back("A");
		letters.push_back("C");
		letters.push_back("T");
		letters.push_back("G");

		while (k >= 0) {
			uint64_t kmin = pow(4, k);
			if (kmin <= realid) {
				uint64_t alpha = 3;
				bool statealpha = true;

				while (statealpha) {
					if (realid >= (alpha * kmin)) {
						realid -= (alpha * kmin);
						statealpha = false;
						result += letters[alpha];
					}
					else {
						alpha -= 1;
					}
				}
			}
			else {
				result += letters[0];
			}
			k -= 1;
		}
		return result;
	}


	/********************************************************************************************
	* Helpers for access to attributes
	* max iteration
	********************************************************************************************/
	/*!
	* \brief get the maximum iteration in the file
	* */
	uint64_t getMaxIteration() {
		ifstream file_ai(_self_file.c_str(), ios::in);
		if (!file_ai) {
			cerr << "Error while opening " << _self_file << " in read mode" << endl;
			exit(EXIT_FAILURE);
		}
		uint64_t read_iter;
		string read_sequences;
		float read_test, read_outer;
		uint64_t iter_max = 0;

		while (file_ai >> read_iter >> read_sequences >> read_test >> read_outer) {
			if (read_iter > iter_max) {
				iter_max = read_iter;
			}
		}
		return (iter_max);
	}


	/********************************************************************************************
	* Model estimation
	* Max values per iter, test vs outer, new kmers per iteration
	********************************************************************************************/
	/**
	 * \brief return the curve of best score
	 * */
	vector<float> getMaxValuePerIter(uint64_t numberOfIterationStep = 1, string p_score = "test", float p_valueFilter = 0.8) {
		uint64_t read_iter;
		string read_sequences;
		float read_test, read_outer;
		uint64_t maxIter = getMaxIteration();
		uint64_t maxBin = ceil((long double)maxIter / numberOfIterationStep) + 1;
		vector<float> vector_result(maxBin, 0);

		ifstream file_ai(_self_file.c_str(), ios::in);
		if (!file_ai) {
			cerr << "Error while opening " << _self_file << " in read mode" << endl;
			exit(EXIT_FAILURE);
		}

		while (file_ai >> read_iter >> read_sequences >> read_test >> read_outer ) {
			
				float score_current = getScores(read_test, read_outer, p_score);
				if (score_current >= p_valueFilter) {
					uint64_t bin_current = floor((long double)read_iter / numberOfIterationStep);

					if (vector_result[bin_current] < score_current) {
						vector_result[bin_current] = score_current;
					}
				}
			
		}

		return vector_result;
	}

	/**
	 * \brief return the curve of best score
	 * */
	vector<vector<float>> getTestOuterList(string p_score = "test", float p_valueFilter = 0.8) {
		uint64_t read_iter;
		string read_sequences;
		float read_test, read_outer;

		vector<vector<float>> vector_result;

		ifstream file_ai(_self_file.c_str(), ios::in);
		if (!file_ai) {
			cerr << "Error while opening " << _self_file << " in read mode" << endl;
			exit(EXIT_FAILURE);
		}

		while (file_ai >> read_iter >> read_sequences >> read_test >> read_outer) {
			
			float score_current = getScores(read_test, read_outer, p_score);
			if (score_current >= p_valueFilter) {
				vector<float> tmp;
				tmp.push_back(read_test);
				tmp.push_back(read_outer);
				vector_result.push_back(tmp);
			}
			
		}

		return vector_result;
	}

	/**
	 * \brief return the number of new kmers every numberOfIterationStep iterations
	 * */
	vector<uint64_t> getNewKmerPerIteration(uint64_t numberOfIterationStep = 10, string p_score = "test", float p_valueFilter = 0.8) {
		ifstream file_ai(_self_file.c_str(), ios::in);
		if (!file_ai) {
			cerr << "Error while opening " << _self_file << " in read mode" << endl;
			exit(EXIT_FAILURE);
		}
		uint64_t read_iter;
		string read_sequences;
		float read_test, read_outer;

		vector<uint64_t> vector_occurences;
		map<string, bool> map_kmerAlreadySeen;
		while (file_ai >> read_iter >> read_sequences >> read_test >> read_outer ) {
			if (getScores(read_test, read_outer, p_score) >= p_valueFilter) {
				//compute the number of new kmers
				uint64_t nbNewKmer = 0;
				vector<string> content;
				boost::split(content, read_sequences, boost::is_any_of(","));
				for (uint64_t i = 0; i < content.size(); i++) {
					if (map_kmerAlreadySeen.find(content[i]) == map_kmerAlreadySeen.end()) {
						map_kmerAlreadySeen[content[i]] = true;
						nbNewKmer++;
					}
				}
				//add this amount of new kmers to the appropriate bin
				if (nbNewKmer > 0) {
					uint64_t bin_current = floor((long double)read_iter / numberOfIterationStep);
					while (bin_current >= vector_occurences.size()) {
						vector_occurences.push_back(0);
					}
					vector_occurences[bin_current] += nbNewKmer;
				}
			}
			
		}
		file_ai.close();
		return vector_occurences;
	}

	/********************************************************************************************
	* Scores
	* Kendall + Hamming
	********************************************************************************************/
	/**
	 * \brief computes the Kendall rank correlation score for 2 SeqAndCount vector with same size
	 * */
	long double computeKendallCorrelation(vector<SeqAndCount> data1, vector<SeqAndCount> data2) {
		sort(data1.begin(), data1.end(), compareByCount);
		sort(data2.begin(), data2.end(), compareByCount);

		long double discord = 0;
		long double concord = 0;
		for (uint64_t i = 0; i < data1.size(); i++) {
			string seqA = data1[i].sequence;
			for (uint64_t j = i + 1; j < data1.size(); j++) {
				string seqB = data1[j].sequence;
				uint64_t kA = 0;
				uint64_t kB = 0;
				uint64_t mushroom = 2;
				for (uint64_t k = 0; k < data2.size(); k++) {
					if (data2[k].sequence.compare(seqA) == 0) {
						kA = k;
						mushroom--;
					}
					if (data2[k].sequence.compare(seqB) == 0) {
						kB = k;
						mushroom--;
					}
					if (mushroom == 0) {
						k = data2.size();
					}
				}

				if (kA > kB) {
					discord++;
				}
				else {
					concord++;
				}
			}

		}
		long double N = data1.size();
		long double tau = (concord - discord) / (N*(N - 1)*0.5);
		return tau;
	}

	/**
	 * \brief computes the Kendall distance between 2 SeqAndCount vector with same size
	 * */
	long double computeKendallDistance(vector<SeqAndCount> data1, vector<SeqAndCount> data2) {
		sort(data1.begin(), data1.end(), compareByCount);
		sort(data2.begin(), data2.end(), compareByCount);

		long double discord = 0;
		for (uint64_t i = 0; i < data1.size(); i++) {
			string seqA = data1[i].sequence;
			for (uint64_t j = i + 1; j < data1.size(); j++) {
				string seqB = data1[j].sequence;
				uint64_t kA = 0;
				uint64_t kB = 0;
				uint64_t mushroom = 2;
				for (uint64_t k = 0; k < data2.size(); k++) {
					if (data2[k].sequence.compare(seqA) == 0) {
						kA = k;
						mushroom--;
					}
					if (data2[k].sequence.compare(seqB) == 0) {
						kB = k;
						mushroom--;
					}
					if (mushroom == 0) {
						k = data2.size();
					}
				}

				if (kA > kB) {
					discord++;
				}
			}

		}
		long double N = data1.size();
		long double tau = discord / (N*(N - 1)*0.5);
		return tau;
	}

	/**
	 * \brief computes the Hamming distance for 2 SeqAndCount vector with same size
	 * */
	uint64_t computeHamming(vector<SeqAndCount> data1, vector<SeqAndCount> data2) {
		uint64_t distance = 0;
		//for(uint64_t i = 0; i < data1.size(); i++){
		//    cerr << data1[i].sequence << "\t" << data2[i].sequence << endl;
		//}

		map<string, bool> maptmp;
		for (uint64_t i = 0; i < data1.size(); i++) {
			maptmp[data1[i].sequence] = true;
		}
		for (uint64_t i = 0; i < data2.size(); i++) {
			if (maptmp.find(data2[i].sequence) == maptmp.end()) {
				distance++;
			}
		}
		//cerr << distance << endl;
		//exit(EXIT_FAILURE);
		return distance;
	}


	/********************************************************************************************
	* Histograms
	*
	********************************************************************************************/
	/**
	 * \brief get the histogram of counts of kmers. Cut the histogram at threshold_deriv
	 * return a vector of kmers with their count, sort in descending order
	* */
	vector<SeqAndCount> getPartialHistogram(string p_score = "test", float p_valueFilter = 0.8, long double threshold_deriv = 0.01) {
		ifstream file_ai(_self_file.c_str(), ios::in);
		if (!file_ai) {
			cerr << "Error while opening " << _self_file << " in read mode" << endl;
			exit(EXIT_FAILURE);
		}
		uint64_t read_iter;
		string read_sequences;
		float read_test, read_outer;
		map<string, uint64_t> map_kmerAlreadySeen;

		//read the kmers
		while (file_ai >> read_iter >> read_sequences >> read_test >> read_outer) {
				if (getScores(read_test, read_outer, p_score) >= p_valueFilter) {
					vector<string> content;
					boost::split(content, read_sequences, boost::is_any_of(","));
					for (uint64_t i = 0; i < content.size(); i++) {
						if (map_kmerAlreadySeen.find(content[i]) == map_kmerAlreadySeen.end()) {
							map_kmerAlreadySeen[content[i]] = 1;
						}
						else {
							map_kmerAlreadySeen[content[i]] += 1;
						}
					}
				}
		}
		file_ai.close();
		//cerr << "getPartialHistogram: found " << map_kmerAlreadySeen.size() << " elements" << endl;

		//get the counts and sort
		//cerr << "getPartialHistogram: start sorting " << endl;
		vector<SeqAndCount> values;
		for (auto it = map_kmerAlreadySeen.begin(); it != map_kmerAlreadySeen.end(); it++) {
			SeqAndCount tmp;
			tmp.sequence = it->first;
			tmp.count = it->second;
			values.push_back(tmp);
		}
		sort(values.begin(), values.end(), compareByCount);
		//cerr << "getPartialHistogram: end sorting " << endl;

		//now computes the derivative to get the stop value, and then delete the other values
		//cerr << "getPartialHistogram: start derivative " << endl;
		uint64_t i_stop = 0;
		bool state_deriv = false;
		if (values.size() > 10) {
			for (uint64_t i = 1; i < (values.size() - 1); i++) {
				long double deriv_c = (values[i + 1].count - values[i - 1].count) / 3;
				if (deriv_c <= threshold_deriv) {
					i_stop = i;
					i = values.size();
					state_deriv = true;
				}
			}
		}
		//cerr << "getPartialHistogram: end derivative " << endl;
		if (state_deriv) {
			values.erase(values.begin() + i_stop, values.end());
		}
		//cerr << "getPartialHistogram: end erasing " << endl;
		return values;
	}

	/**
	 * \brief get the histogram
	 * return a vector of kmers with their count, sort in descending order, limited to p_NFirstValues first values and iteration <= iter_max
	* */
	vector<SeqAndCount> getConstraintHeadOfCompleteHistogram(uint64_t iter_max, uint64_t p_NFirstValues = 10, string p_score = "test", float p_valueFilter = 0.8) {
		ifstream file_ai(_self_file.c_str(), ios::in);
		if (!file_ai) {
			cerr << "Error while opening " << _self_file << " in read mode" << endl;
			exit(EXIT_FAILURE);
		}
		uint64_t read_iter;
		string read_sequences;
		float read_test, read_outer;
		map<string, uint64_t> map_kmerAlreadySeen;

		//read the kmers
		while (file_ai >> read_iter >> read_sequences >> read_test >> read_outer) {
			if (getScores(read_test, read_outer, p_score) >= p_valueFilter) {
				if (read_iter <= iter_max) {
					vector<string> content;
					boost::split(content, read_sequences, boost::is_any_of(","));
					for (uint64_t i = 0; i < content.size(); i++) {
						if (map_kmerAlreadySeen.find(content[i]) == map_kmerAlreadySeen.end()) {
							map_kmerAlreadySeen[content[i]] = 1;
						}
						else {
							map_kmerAlreadySeen[content[i]] += 1;
						}
					}
				}
			}
		}
		file_ai.close();

		//get the counts and sort
		vector<SeqAndCount> values;
		for (auto it = map_kmerAlreadySeen.begin(); it != map_kmerAlreadySeen.end(); it++) {
			SeqAndCount tmp;
			tmp.sequence = it->first;
			tmp.count = it->second;
			values.push_back(tmp);
		}
		sort(values.begin(), values.end(), compareByCount);
		if (values.size() < p_NFirstValues) {
			p_NFirstValues = values.size();
		}
		values.erase(values.begin() + p_NFirstValues, values.end());
		return values;
	}


	/********************************************************************************************
	* Convergence
	*
	********************************************************************************************/
	bool getConvergence(string path_outputfile){
		vector<uint64_t> IX;
		vector<uint64_t> IY;
		uint64_t nbIterMax = getMaxIteration();
		ofstream file_ao(path_outputfile, ios::out);
		file_ao << "iteration\tdistance" << endl;

		cerr << "computing partial histogram + Hamming suite " << endl;
		vector<SeqAndCount> t1_ref = getPartialHistogram("test", 0); 
		cerr << "t1_ref size = " << t1_ref.size() << endl;

		for (uint64_t iter_c = 500; iter_c < (nbIterMax-500); iter_c += 500) {
			vector<SeqAndCount> t1_test = getConstraintHeadOfCompleteHistogram(iter_c, t1_ref.size(), "test", 0);//get the histogram at a given iteration for a given set of kmers
			uint64_t distance = computeHamming(t1_ref, t1_test);
			file_ao << iter_c << "\t" << distance << endl;
			IX.push_back(iter_c);
			IY.push_back(distance);
		}
		file_ao.close();


		bool result = false;

		for (uint64_t pos = 1; pos < (IX.size() - 1); pos++) {
			float slope = abs((float)(IY[pos + 1] - IY[pos - 1])) / ((float)(IX[pos + 1] - IX[pos - 1]));
			if (slope < 0.01) {
				result = true;
			}
			else {
				result = false;
			}
		}
		return result;
	}
};
