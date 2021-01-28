#pragma once
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

#ifndef MATRIX_READER
#define MATRIX_READER


class MatrixReader {


public:		
	vector<vector<double>> ReadMatris(string path);
	vector<double> ReadVector(string path);
};

#endif // !MATRIX_READER
