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

void MatrixCalculator::VektorYazdir(vector<double> vektor)
{
	const auto m = static_cast<int>(vektor.size());

	for (auto i = 0; i < m; ++i) {
		cout << "\t" << vektor[i] << "  ";
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

vector<vector<double>> MatrixCalculator::CholeskyFactorHesapla(vector<vector<double>> matris)
{
	const auto n = static_cast<int>(matris.size());
	
	vector<vector<double>> result(n, vector<double>(n));

	for (auto i = 0; i < n; i++) {
		for (auto j = 0; j <= i; j++) {		
						
			if (j == i)
			{
				double sum = 0;
				for (auto k = 0; k < j; k++)
					sum += result[j][k] * result[j][k];

				result[j][j] = sqrt(matris[j][j] - sum);
			}
			else
			{
				double sum = 0;
				for (auto k = 0; k < j; k++)
					sum += result[i][k] * result[j][k];

				result[i][j] = 1.0 / result[j][j] * (matris[i][j] - sum);
			}
		}
	}

	return result;
}

double MatrixCalculator::HadamardKatsayisiHesapla(vector<vector<double>> matris)
{
	const auto determinant = abs(this->DeterminantHesapla(matris));
	auto satirNormlariToplami = 0.0;
	const auto m = static_cast<int>(matris.size());

	for (auto i = 0; i < m; ++i)
		satirNormlariToplami += this->SatirNorm(matris, i);

	return determinant / satirNormlariToplami;	
}

vector<vector<double>> MatrixCalculator::Inverse(vector<vector<double>> matris)
{
	double d = 1.0 / this -> DeterminantHesapla(matris);
	vector<vector<double>> solution(matris.size(), std::vector<double>(matris.size()));

	for (size_t i = 0; i < matris.size(); i++) {
		for (size_t j = 0; j < matris.size(); j++) {
			solution[i][j] = matris[i][j];
		}
	}

	solution = this->TranspozeHesapla(this->CofactorHesapla(solution));

	for (size_t i = 0; i < matris.size(); i++) {
		for (size_t j = 0; j < matris.size(); j++) {
			solution[i][j] *= d;
		}
	}

	return solution;
}

vector<vector<double>> MatrixCalculator::CofactorHesapla(const vector<vector<double>> matris) {

	vector<vector<double>> solution(matris.size(), vector<double>(matris.size()));
	vector<vector<double>> subVect(matris.size() - 1, vector<double>(matris.size() - 1));

	for (auto i = 0; i < matris.size(); i++) {
		for (auto j = 0; j < matris[0].size(); j++) {
			auto p = 0;

			for (auto x = 0; x < matris.size(); x++) {
				if (x == i) {
					continue;
				}
				auto q = 0;

				for (size_t y = 0; y < matris.size(); y++) {
					if (y == j) {
						continue;
					}

					subVect[p][q] = matris[x][y];
					q++;
				}
				p++;
			}
			solution[i][j] = pow(-1, i + j) * this->DeterminantHesapla(subVect);
		}
	}
	return solution;
}

vector<double> MatrixCalculator::GaussJordanBilinmeyenVektorHesapla(vector<vector<double>> matris)
{
	vector<double> solution(vector<double>(matris.size()));
	vector<vector<double>> temp(matris.size(), vector<double>(matris.size()));
	const auto n = static_cast<int>(matris.size());
	
	for (auto j = 0; j < n; j++)
	{
		for (auto i = 0; i < n; i++)
		{
			if(i!=j)
			{
				const auto c = matris[i][j] / matris[j][j];
				
				for (auto k = 0; k < n; k++) {
					temp[i][k] = temp[i][k] - c * matris[j][k];
				}
			}
		}		
	}

	for (auto i = 0; i < n; i++)
		solution[i] = temp[i][0] / temp[i][i];

	return solution;
}
