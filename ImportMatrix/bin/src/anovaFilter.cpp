#include "cephes/mconf.h"
#include <boost/algorithm/string.hpp>
#include <string>
#include <vector>
#include <map>
#include <numeric>
#include <math.h>
#include <algorithm>
#include <cmath>
#include <stdint.h>
#include <fstream>
#include <iostream>
#include <omp.h>

extern "C" {
double fdtrc(int ia, int ib, double x);
}


class Stats {
public:
	Stats();
	virtual ~Stats();

	static double long stdev(std::vector<double long> v) {
		return stdev(v, mean(v));
	}
	static double long stdev(std::vector<double long> v, double long mean) {
		std::vector<double> diff(v.size());
		std::transform(v.begin(), v.end(), diff.begin(),
				[mean](double x) {return x - mean;});
		double long sq_sum = std::inner_product(diff.begin(), diff.end(),
				diff.begin(), 0.0);
		return std::sqrt(sq_sum / v.size());
	}

	static double long mean(std::vector<double long> v) {
		double long sum = std::accumulate(v.begin(), v.end(), 0.0);
		return sum / (double long) v.size();
	}
	template <class T>
	static void rescale(std::vector<T> & values) {
		long double max = 0;
		for (auto v : values)
			max = v > max ? v : max;
		for (uint64_t i = 0; i < values.size(); i++) {
			values[i] = (values[i] * 100) / max;
		}
	}

	static std::vector<int64_t> distanceFromMean(std::string file,
			std::map<std::string, std::vector<uint64_t>> groupMap);
	static std::pair<double long, double long> f_test(std::vector<double long> & vect,
			std::vector<uint64_t> & groups,
			std::map<uint64_t, uint64_t> & grp_counts);
	static double long pearsonCorrelationCoefficient(std::vector<double long>,
			std::vector<double long>);
	static double long median(std::vector<double long>);

};


Stats::Stats() {
	// TODO Auto-generated constructor stub

}

Stats::~Stats() {
	// TODO Auto-generated destructor stub
}

///
/// Adapted from Python f_oneway function [sklearn]. It assumes the groups are ints in [0,n_of_groups)
///
std::pair<double long, double long> Stats::f_test(std::vector<double long> & vect,
		std::vector<uint64_t> & groups,
		std::map<uint64_t, uint64_t> & grp_counts) {

	uint64_t n_classes = grp_counts.size();
	double long n_samples = vect.size();
	//  n_samples_per_class = np.array([a.shape[0] for a in args]) -> grp_counts
	//  n_samples = np.sum(n_samples_per_class)  -> vect.size()
	double long ss_alldata = 0;
	double long summ_all = 0;
	std::vector<double long> sums_args(n_classes, 0.);

	for (uint64_t i = 0; i < vect.size(); i++) {
		ss_alldata += (vect[i] * vect[i]);
		sums_args[groups[i]] += vect[i];
		summ_all += (vect[i]);

	}

	double long square_of_sums_alldata = summ_all * summ_all;
	std::vector<double long> square_of_sums_args;
	double long ssbn = 0;
	for (uint64_t i = 0; i < n_classes; i++) {
		square_of_sums_args.push_back(sums_args[i] * sums_args[i]);
		ssbn += (square_of_sums_args[i] / grp_counts[i]);
	}

	ssbn -= (square_of_sums_alldata / n_samples);
	double long sstot = (ss_alldata - (square_of_sums_alldata / n_samples));
	double long sswn = sstot - ssbn;
	int dfbn = n_classes - 1;
	int dfwn = n_samples - n_classes;
	double long msb = ssbn / (double long) (n_classes - 1);
	double long msw = sswn / (double long) dfwn;
	double long f = msb / msw;
	double long prob = fdtrc(dfbn, dfwn, f);
	return std::pair<double long, double long>(f, prob);
}

double long sum(std::vector<double long> a) {
	double s = 0;
	for (int i = 0; i < a.size(); i++) {
		s += a[i];
	}
	return s;
}

double long mean(std::vector<double long> a) {
	return sum(a) / a.size();
}

double long sqsum(std::vector<double long> a) {
	double s = 0;
	for (int i = 0; i < a.size(); i++) {
		s += pow(a[i], 2);
	}
	return s;
}

std::vector<double long> operator-(std::vector<double long> a, double long b) {
	std::vector<double long> retvect;
	for (int i = 0; i < a.size(); i++) {
		retvect.push_back(a[i] - b);
	}
	return retvect;
}

std::vector<double long> operator*(std::vector<double long> a, std::vector<double long> b) {
	std::vector<double long> retvect;
	for (int i = 0; i < a.size(); i++) {
		retvect.push_back(a[i] * b[i]);
	}
	return retvect;
}

/// from https://codepad.co/snippet/MbadrcBL
double long Stats::pearsonCorrelationCoefficient(std::vector<double long> X,
		std::vector<double long> Y) {
	return sum((X - mean(X)) * (Y - mean(Y))) / (X.size() * stdev(X) * stdev(Y));
}

double long Stats::median(std::vector<double long> X) {
	std::sort(X.begin(), X.end());
	return X.size() % 2 == 0 ?
			X[X.size() / 2] : ((X[X.size() / 2] + X[(X.size() / 2) + 1]) / 2);
}





int main(int argc, char** argv) {
	if (argc < 4 ){
		std::cerr << "usage:\n\tanovaFilter <input matrix> <output matrix> <p-value threshold>\n\n";
		exit(1);
	}
	std::string inFname = std::string(argv[1]);
	std::string outFname = std::string(argv[2]);
	double long thr = std::stold(argv[3]);
	std::cout << "Input Matrix: " << inFname << "\n";
	std::cout << "Output Matrix: " << outFname << "\n";
	std::cout << "P-value Threshold: " << thr << "\n";
	std::string line;
	std::ofstream ostr(outFname);
	std::ifstream instr(inFname);
	std::vector<std::string> buffer, col_groups, col_names;
	uint64_t bufferSize = omp_get_max_threads() * 100000;
	bool running = true;
	getline(instr, line);
	boost::split(col_names, line, boost::is_any_of("\t"));
	ostr << line << "\n";
	col_names.erase(col_names.begin());
	getline(instr, line);
	ostr << line << "\n";
	boost::split(col_groups, line, boost::is_any_of("\t"));
	col_groups.erase(col_groups.begin());
	std::map < uint64_t, uint64_t > group_counts;
	std::map<std::string, uint64_t> int_groups;
	std::vector<uint64_t> groups;
	uint64_t ig=0;
	for (auto g : col_groups) {
		if (int_groups.count(g) == 0) {
			int_groups[g] = ig++;
			group_counts[int_groups[g]] = 0;
		}
		groups.push_back(int_groups[g]);
		group_counts[int_groups[g]]++;
	}

	while (running) {
		running = getline(instr, line) ? true : false;
		if (running) {
			buffer.push_back(line);
		}
		if (buffer.size() == bufferSize || !running) {
			uint64_t bs = buffer.size();
#pragma omp parallel for firstprivate(groups, group_counts, thr)
			for (uint64_t i = 0; i < bs; i++) {
				std::vector<std::string> content;
				boost::split(content, buffer[i], boost::is_any_of("\t"));
				std::vector<double long> values;
				for (uint64_t j=1; j<content.size(); j++) values.push_back(std::stold(content[j]));
				std::pair<double long, double long> res = Stats::f_test( values, groups, group_counts);
				if (res.second > thr ) {
					buffer[i]="";
				} else {
					buffer[i]+="\n";
				}
			}
			for (uint64_t i = 0; i < bs; i++)
				ostr << buffer[i];
			buffer.clear();
		}

	}
	ostr.close();
	std::cout << "Done.\n";
	return 0;

}
