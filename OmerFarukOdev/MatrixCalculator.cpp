#include "MatrixCalculator.h"
#include <iostream>
#include <vector>

using namespace std;

double MatrixCalculator::DeterminantHesapla(vector<vector<double>> matris, int n)
{
	vector<vector<double>> temp(matris.size(), vector<double>(matris[0].size()));
	double determinant = 0;

	if (n == 0)
		n = static_cast<int>(matris.size());

	if (n == 2)
		return ((matris[0][0] * matris[1][1]) - (matris[1][0] * matris[0][1]));
	else {
		for (auto x = 0; x < n; x++) {
			auto subi = 0;
			for (auto i = 1; i < n; i++) {
				auto subj = 0;
				for (auto j = 0; j < n; j++) {
					if (j == x)
						continue;
					temp[subi][subj] = matris[i][j];

					subj++;
				}
				subi++;
			}
			determinant = determinant + (pow(-1, x) * matris[0][x] * DeterminantHesapla(temp, n - 1));
		}
	}

	return determinant;
}

double MatrixCalculator::IzHesapla(vector<vector<double>> matris)
{
	double sum = 0;
	const auto n = static_cast<int>(matris.size());

	for (auto i = 0; i < n; i++)
		sum += matris[i][i];

	return sum;
}

double MatrixCalculator::SutunNorm(vector<vector<double>> matris, int sutunNo) {
	double sum = 0;
	const auto n = static_cast<int>(matris.size());

	for (auto i = 0; i < n; i++)
		sum += pow(matris[i][sutunNo], 2);

	return sqrt(sum);
}

double MatrixCalculator::SatirNorm(vector<vector<double>> matris, int satirNo) {
	double sum = 0;
	const auto n = static_cast<int>(matris[0].size());

	for (auto i = 0; i < n; i++)
		sum += pow(matris[satirNo][i], 2);

	return sqrt(sum);
}

double MatrixCalculator::SpektralKondKatsayisi(vector<vector<double>> matris)
{
	const auto determinant = abs(this->DeterminantHesapla(matris));
	auto satirNormlariToplami = 0.0;
	const auto m = static_cast<int>(matris.size());

	for (auto i = 0; i < m; ++i)
		satirNormlariToplami += this -> SatirNorm(matris, i);

	return determinant / satirNormlariToplami;
}

double MatrixCalculator::OklidNorm(vector<vector<double>> matris) {
	double sum = 0;
	const auto m = static_cast<int>(matris.size());
	const auto n = static_cast<int>(matris[0].size());

	for (auto i = 0; i < n; i++)
		for (auto j = 0; j < m; j++)
			sum += pow(matris[i][j], 2);

	return sqrt(sum);
}


vector<vector<double>> MatrixCalculator::MatrisCarpimHesapla(vector<vector<double>> m1, vector<vector<double>> m2) {

	vector<vector<double>> carpim(m1.size(), vector<double>(m1[0].size()));
	const auto m = static_cast<int>(m1.size());
	const auto n = static_cast<int>(m1[0].size());

	for (auto i = 0; i < m; ++i)
		for (auto j = 0; j < n; ++j)
			for (auto k = 0; k < n; ++k)
				carpim[i][j] += m1[i][k] * m2[k][j];

	return carpim;
}

vector<vector<double>> MatrixCalculator::TranspozeHesapla(vector<vector<double>> matris) {

	vector<vector<double>> transpose(matris.size(), vector<double>(matris[0].size()));
	const auto m = static_cast<int>(matris.size());
	const auto n = static_cast<int>(matris[0].size());

	for (auto i = 0; i < m; ++i)
		for (auto j = 0; j < n; ++j) {
			transpose[j][i] = matris[i][j];
		}

	return transpose;
}

vector<vector<double>> MatrixCalculator::MartisNormlastir(vector<vector<double>> matris, double norm) {

	vector<vector<double>> normlastirilmis(matris.size(), vector<double>(matris[0].size()));
	const auto m = static_cast<int>(matris.size());
	const auto n = static_cast<int>(matris[0].size());

	for (auto i = 0; i < m; ++i) {
		for (auto j = 0; j < n; ++j) {
			normlastirilmis[i][j] = matris[i][j] / norm;
		}
	}

	return normlastirilmis;
}

void MatrixCalculator::MatrisYazdir(vector<vector<double>> matris) {
	const auto m = static_cast<int>(matris.size());
	const auto n = static_cast<int>(matris[0].size());

	for (auto i = 0; i < m; ++i) {
		for (auto j = 0; j < n; ++j) {
			cout << "\t" << matris[i][j] << "\t";
		}
		cout << "\n";
	}
}

vector<vector<double>> MatrixCalculator::GausMatrisTersi(vector<vector<double>> matris)
{
	auto height = static_cast<int>(matris.size());
	auto width = static_cast<int>(matris[0].size());

	vector<vector<double>> result(matris.size(), vector<double>(matris[0].size()));

	// identity matris oluþtur
	for (auto i = 0; i < width; i++)
		result[i][i] = 1;

	// Eþalon form
	for (auto j = 0; j < width; ++j) {
		// partial pivoting
		auto maxRow = j;
		for (auto i = j; i < height; ++i) {
			maxRow = matris[i][j] > matris[maxRow][j] ? i : maxRow;
		}

		matris[j].swap(matris[maxRow]);
		result[j].swap(result[maxRow]);

		// Satýr iþlemleri
		auto pivot = matris[j][j];
		auto& row1L = matris[j];
		auto& row1R = result[j];
		for (auto i = j + 1; i < height; ++i) {
			auto& row2L = matris[i];
			auto& row2R = result[i];
			auto temp = row2L[j];
			for (auto k = 0; k < width; ++k) {
				row2L[k] -= temp / pivot * row1L[k];
				row2R[k] -= temp / pivot * row1R[k];
			}
		}

		// Diyagonal elemanlar 1
		for (auto k = 0; k < width; ++k) {
			row1L[k] /= pivot;
			row1R[k] /= pivot;
		}
	}
	
	for (auto j = width - 1;; --j) {
		auto& row1L = matris[j];
		auto& row1R = result[j];
		for (auto i = 0; i < j; ++i) {
			auto& row2L = matris[i];
			auto& row2R = result[i];
			auto temp = row2L[j];
			for (auto k = 0; k < width; ++k) {
				row2L[k] -= temp * row1L[k];
				row2R[k] -= temp * row1R[k];
			}
		}
		
		if (j == 0) break;
	}
		
	return result;
}
