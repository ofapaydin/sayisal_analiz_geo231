#include "MatrixReader.h"
#include <string>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <iomanip>

using namespace std;

vector<vector<double>> MatrixReader::ReadMatris(string path)
{
	ifstream file;
	vector<vector<double>> data;
	string line;
	
	file.open(path);

	if (file.fail()) {
		cerr << "Matrix dosya hatasi." << endl;
		exit(1);
	}


	while (getline(file, line)) {
		istringstream reader(line);
		
		vector<double> lineData;
		
		string::const_iterator i = line.begin();

		while (!reader.eof()) {
			double val;
			reader >> val;
			
			if (reader.fail())
				break;

			lineData.push_back(val);
		}

		data.push_back(lineData);
	}

	return data;
}

vector<double> MatrixReader::ReadVector(string path)
{
	ifstream file;
	vector<double> data;
	string line;

	file.open(path);

	if (file.fail()) {
		cerr << "Vector dosya hatasi." << endl;
		exit(1);
	}


	while (getline(file, line)) {
		istringstream reader(line);

		string::const_iterator i = line.begin();

		while (!reader.eof()) {
			double val;
			reader >> val;

			if (reader.fail())
				break;

			data.push_back(val);
		}

	}

	return data;
}

