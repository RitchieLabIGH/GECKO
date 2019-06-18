#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <map>
#include <algorithm>
#include <chrono>
#include <boost/algorithm/string.hpp>

class Split {
public:
	Split(std::string inName, std::string baseName) {
		_inFile = inName;
		_baseName = baseName;
		parseMatrixHeader(_inFile, col_names, col_groups);
	}
	virtual ~Split() {
	}
	;
	void flatSplitCross(uint64_t nrows, uint64_t rowsLimit);
	void meansGroupsSplit(uint64_t, uint64_t, uint64_t);

private:
	std::vector<std::vector<uint64_t>> makeSplitGroupsByMAD(uint64_t);
	std::vector<std::string> splitFile(std::string filename,
			std::string baseName, uint64_t nrows);
	uint64_t countFileLines(std::string & file);
	void parseMatrixHeader(std::string & file, std::vector<std::string> & names,
			std::vector<std::string> & groups);
	void printMatrixHeader(std::string & file,
			std::vector<std::string> & col_names,
			std::vector<std::string> & col_groups);
	std::string fillZero(uint64_t);
	std::string _inFile;
	std::string _baseName;
	std::vector<std::string> col_names;
	std::vector<std::string> col_groups;

};

uint64_t Split::countFileLines(std::string & file) {
	std::ifstream fin(file.c_str());
	if (!fin) {
		std::cerr << "countFileLines: Error while opening " << file
				<< " in read mode" << std::endl;
		exit (EXIT_FAILURE);
	}
	uint64_t count = 0;
	const int SZ = 1024 * 1024;
	std::vector<char> buff(SZ);
	fin.read(&buff[0], buff.size());
	const char * p;
	while (int cc = fin.gcount()) {
		p = &buff[0];
		for (int i = 0; i < cc; i++) {
			if (p[i] == '\n')
				count++;
		}
		fin.read(&buff[0], buff.size());
	}
	fin.close();
	return count;
}

void Split::flatSplitCross(uint64_t nrows, uint64_t maxrows) {
	std::vector < std::string > files = splitFile(_inFile, _baseName, nrows);
	uint64_t flen;
	for (auto it : files) {
		flen = countFileLines(it);
		if (flen > maxrows) {
			std::vector < std::string > files2 = splitFile(it, it + "_part",
					maxrows);
			remove(it.c_str());
		}
	}
}

std::vector<std::vector<uint64_t>> Split::makeSplitGroupsByMAD(
		uint64_t maxGroupSize) {
	std::vector < std::vector < uint64_t >> split_groups;
	std::map<std::string, std::vector<uint64_t>> groupMap;
	for (uint64_t g = 0; g < col_groups.size(); g++) {
		if (groupMap.count(col_groups[g]) == 0)
			groupMap[col_groups[g]] = std::vector<uint64_t>();
		groupMap[col_groups[g]].push_back(g);
	}

	std::vector<double long> distances(col_names.size(), 0);
	std::ifstream inFile;
	std::string line;
	inFile.open(_inFile.c_str(), std::ifstream::in);
	getline(inFile, line);
	if (line[0] == '\t') {
		getline(inFile, line); // consume header
	} else {
		inFile.close();
		inFile.open(_inFile.c_str(), std::ifstream::in);
	}
	std::vector < std::string > content;
	bool running = getline(inFile, line) ? true : false;
	uint64_t maxRows = 1000000;
	while (running) {
		maxRows--;
		boost::split(content, line, boost::is_any_of("\t"));
		std::vector<int> currRow;
		for (uint64_t i = 1; i < content.size(); i++) {
			currRow.push_back(std::stoi(content[i]));
		}
		for (auto gm : groupMap) {
			std::vector<int> data;
			for (uint64_t col : gm.second)
				data.push_back(currRow[col]);
			std::sort(data.begin(), data.end());
			float median =
					data.size() % 2 == 0 ?
							data[data.size() / 2] :
							((data[data.size() / 2]
									+ data[(data.size() / 2) + 1]) / 2);
			for (uint64_t i = 0; i < data.size(); i++)
				data[i] = std::abs(data[i] - median);
			std::sort(data.begin(), data.end());
			float MAD =
					data.size() % 2 == 0 ?
							data[data.size() / 2] :
							((data[data.size() / 2]
									+ data[(data.size() / 2) + 1]) / 2);
			for (uint64_t col : gm.second) {
				distances[col] += ((currRow[col] - median) / MAD);
			}
		}
		running = getline(inFile, line) ? true : false;
		if (maxRows == 0)
			running = false;
	}
	inFile.close();

	for (auto gm : groupMap) {
		std::sort(gm.second.begin(), gm.second.end(),
				[distances](uint64_t a, uint64_t b) {return distances[a] < distances[b];});
		split_groups.push_back(std::vector<uint64_t>());
		for (uint64_t col : gm.second) {
			if (split_groups.back().size() == maxGroupSize)
				split_groups.push_back(std::vector<uint64_t>());
			split_groups.back().push_back(col);

		}
	}
	for (uint64_t g = 0; g < split_groups.size(); g++) {
		std::cout << " G= " << g << "\n";
		for (uint64_t col : split_groups[g])
			std::cout << "\t" << col << "=" << col_groups[col] << " "
					<< col_names[col] << "\n";
	}
	return split_groups;
}


struct MatrixLine {
	std::string line;
	std::string file;
};


void Split::meansGroupsSplit(uint64_t maxGroupSize = 100, uint64_t rowlimit =
		1000000, uint64_t filesLimit = 5000) {
	auto start = std::chrono::system_clock::now();
	std::vector < std::vector < uint64_t >> split_groups = makeSplitGroupsByMAD(
			maxGroupSize);
	std::vector < std::vector<std::vector<uint64_t>> > buckets;
	std::ifstream inFile;
	std::string line;
	inFile.open(_inFile.c_str(), std::ifstream::in);
	// consume header
	getline(inFile, line);
	if (line[0] == '\t') {
		getline(inFile, line); // consume header
	} else {
		inFile.close();
		inFile.open(_inFile.c_str(), std::ifstream::in);
	}
	std::vector<std::map<uint64_t, uint64_t>> freqs;
	for (uint64_t i = 0; i < split_groups.size(); i++) {
		std::map < uint64_t, uint64_t > nmap;
		freqs.push_back(nmap);
	}
	std::vector < std::string > content;
	uint64_t mean;
	std::cout << "Creating the buckets...\n";
	uint64_t maxRows = 1000000;
	while (getline(inFile, line) && maxRows > 0) {
		maxRows--;
		boost::split(content, line, boost::is_any_of("\t"));
		for (uint64_t i = 0; i < split_groups.size(); i++) {
			mean = 0;
			for (uint64_t j = 0; j < split_groups[i].size(); j++) {
				mean += std::stoi(content[split_groups[i][j] + 1]);
			}
			freqs[i][mean]++;
		}
	}
	inFile.close();

	uint64_t max;

	int t;
	/// aggregate in order to be closer to the most frequent class

	for (uint64_t i = 0; i < freqs.size(); i++) {
		std::cout << "\n\n----- \ng " << i << "\n";
		max = 0;
		for (auto it : freqs[i]) {
			std::cout << "mean = " << it.first << " : freq = " << it.second
					<< "\n";
			max = max > it.second ? max : it.second;
		}

		std::vector < std::vector < uint64_t >> v0;
		buckets.push_back(v0);
		std::vector < uint64_t > v1;
		buckets[i].push_back(v1);
		t = 0;
		std::cout << "\nbucket: " << buckets[i].size() << " : ";
		for (auto it : freqs[i]) {
			if (t + it.second > max) {
				t = it.second;
				std::vector < uint64_t > v;
				buckets[i].push_back(v);
				std::cout << "\n bucket: " << buckets[i].size() << " : ";
			} else {
				t = t + it.second;
			}
			buckets[i][buckets[i].size() - 1].push_back(it.first);
			std::cout << " " << it.first << ", ";
		}
	}

	std::cout << "\nDone in " << std::chrono::duration_cast
			< std::chrono::seconds
			> (std::chrono::system_clock::now() - start).count()
					<< " s.\nSplitting the file...\n";
	std::cout.flush();
	start = std::chrono::system_clock::now();
	inFile.open(_inFile.c_str(), std::ifstream::in);
	// consume header
	getline(inFile, line);
	if (line[0] == '\t') {
		getline(inFile, line); // consume header
	} else {
		inFile.close();
		inFile.open(_inFile.c_str(), std::ifstream::in);
	}
	std::vector < uint64_t > means(split_groups.size());
	std::ofstream oFile;
	uint64_t lastF = 0;
	std::map < std::string, std::string > fnames;
	std::map < std::string, uint64_t > fcounts;
	std::vector < MatrixLine > buffer;
	uint64_t bufferSize = 100000;
	bool running = getline(inFile, line) ? true : false;
	while (running) {
		MatrixLine row;
		row.line=line;
		buffer.push_back(row);
		running = getline(inFile, line) ? true : false;
		if (!running || bufferSize == buffer.size()) {
#pragma omp parallel for private(content)
			for (uint64_t l = 0; l < buffer.size(); l++) {
				boost::split(content, buffer[l].line, boost::is_any_of("\t"));
				std::string of = "";
				for (uint64_t i = 0; i < split_groups.size(); i++) {
					means[i] = 0;
					for (uint64_t j = 0; j < split_groups[i].size(); j++) {
						means[i] += std::stoi(content[split_groups[i][j] + 1]);
					}
					bool found = false;
					for (uint64_t j = 0; j < buckets[i].size(); j++) {
						for (uint64_t k = 0; k < buckets[i][j].size(); k++) {
							if (means[i] == buckets[i][j][k]) {
								of += std::string(".") + std::to_string(j);
								found = true;
								break;
							}
						}
						if (found)
							break;
					}
					if (!found) { // the mean was not in the first maxRows rows, assign it to the last bucket
						of += std::string(".")
								+ std::to_string(buckets[i].size()-1);
					}
				}
				buffer[l].file = of;
			}

			for (uint64_t l = 0; l < buffer.size(); l++) {
				if (fnames.count(buffer[l].file) == 0) {
					if (lastF < filesLimit) {
						fnames[buffer[l].file] = _baseName + std::string("_")
								+ fillZero(lastF) + ".matrix";
						printMatrixHeader(fnames[buffer[l].file], col_names,
								col_groups);
						fcounts[fnames[buffer[l].file]] = 0;
						lastF++;
					} else {
						std::string minFile = "";
						for (auto el : fcounts) {
							if (fcounts.count(minFile) == 0
									|| fcounts[minFile] > el.second) {
								minFile = el.first;
							}
						}
						fnames[buffer[l].file] = minFile;
					}
				}
				fcounts[fnames[buffer[l].file]]++;
				oFile.open(fnames[buffer[l].file],
						std::ofstream::out | std::ofstream::app);
				oFile << buffer[l].line << "\n";
				oFile.close();
			}
			buffer.clear();
		}
	}
	inFile.close();
	for (auto it : fcounts) {
		if (it.second > rowlimit) {
			uint64_t balancedRowLimit = it.second
					/ (std::ceil((double) it.second / rowlimit));
			std::vector < std::string > files = splitFile(it.first,
					it.first + "_part", balancedRowLimit);
			remove(it.first.c_str());
		}
	}
	std::cout << "Done in " << std::chrono::duration_cast < std::chrono::seconds
			> (std::chrono::system_clock::now() - start).count() << " s.\n";
}

std::vector<std::string> Split::splitFile(std::string filename,
		std::string baseName, uint64_t nrows) {
	std::ifstream inFile(filename.c_str(), std::ifstream::in);
	std::string line;
	getline(inFile, line);
	if (line[0] == '\t') {
		getline(inFile, line); // consume header
	} else {
		inFile.close();
		inFile.open(filename.c_str(), std::ifstream::in);
	}
	uint64_t c = nrows;
	std::ofstream oFile;
	std::string fname;
	std::vector < std::string > files;
	while (getline(inFile, line)) {
		if (c == nrows) {
			c = 0;
			oFile.close();
			fname = baseName + "_" + fillZero(files.size()) + ".matrix";
			printMatrixHeader(fname, col_names, col_groups);
			oFile.open(fname, std::ofstream::app | std::ofstream::out);
			files.push_back(fname);
		}
		oFile << line << '\n';
		c++;
	}
	oFile.close();
	inFile.close();
	return files;
}

void Split::parseMatrixHeader(std::string & file,
		std::vector<std::string> & outNames,
		std::vector<std::string> & outGroups) {
	std::ifstream fin(file.c_str());
	if (!fin) {
		std::cerr << "IOTools::parseMatrixHeader : Error while opening " << file
				<< " in read mode" << std::endl;
		exit (EXIT_FAILURE);
	}
	std::string line;
	getline(fin, line);
	if (line[0] != '\t') {
		std::cerr << "IOTools::parseMatrixHeader: ERROR! wrong file " << file
				<< ", no header found\n";
		exit(-11);
	}
	boost::split(outNames, line, boost::is_any_of("\t"));
	outNames.erase(outNames.begin());
	getline(fin, line);
	boost::split(outGroups, line, boost::is_any_of("\t"));
	outGroups.erase(outGroups.begin());
	fin.close();
}

void Split::printMatrixHeader(std::string & file,
		std::vector<std::string> & col_names,
		std::vector<std::string> & col_groups) {
	std::ofstream of(file);
	if (!of) {
		std::cerr << "IOTools::printMatrixHeader: Error while opening " << file
				<< " in write mode" << std::endl;
		exit (EXIT_FAILURE);
	}
	for (auto it : col_names)
		of << '\t' << it;
	of << '\n' << "group";
	for (auto it : col_groups)
		of << '\t' << it;
	of << '\n';
	of.close();
}

std::string Split::fillZero(uint64_t n) {
	std::string out = "";
	for (uint64_t i = 10; i < 1000000; i *= 10) {
		if (n < i) {
			out += "0";
		}
	}
	return out + std::to_string(n);

}

int main(int argc, char** argv) {
	if (argc != 6 || argv[1] == "help") {
		std::cout << argv[0]
				<< " <input matrix> <base name> <max group size> <row limit> <files limit>\n"
						"\t- <input matrix>: the input matrix (samples x dimensions) with the first two rows of header.\n"
						"\t- <base name>: the base name of the submatrices. They will be in the form of: <base name>_NNNNN.matrix \n"
						"\t- <max group size>: the maximum number of sample that a group can contain.  \n"
						"\t- <row limit>: maximum number of rows that a submatrix can have.\n\t The submatrices that exceed this limit will be splitted in submatrices with homogeneous number of rows\n"
						"\t- <files limit>: maximum number of files (not considering the additional division caused by row limit exceed). Rare combinations will be placed in the same files. \n";
		exit(0);
	}
	std::string inName = std::string(argv[1]);
	std::string baseName = std::string(argv[2]);
	Split sp(inName, baseName);
	uint64_t maxGroupSize = std::stoi(argv[3]);
	uint64_t rowlimit = std::stoi(argv[4]);
	uint64_t fileLimit = std::stoi(argv[5]);
	sp.meansGroupsSplit(maxGroupSize, rowlimit, fileLimit);
	return 0;

}
