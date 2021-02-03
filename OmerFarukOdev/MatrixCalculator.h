#pragma once
#include <vector>

using namespace std;

class MatrixCalculator
{
public:
	double DeterminantHesapla(vector<vector<double>> matris, int n = 0);
	static double IzHesapla(vector<vector<double>> matris);
	static double SatirNorm(vector<vector<double>> matris, int satirNo);
	static double SutunNorm(vector<vector<double>> matris, int sutunNo);
	static double OklidNorm(vector<vector<double>> matris);
	static vector<vector<double>> MatrisCarpimHesapla(vector<vector<double>> m1, vector<vector<double>> m2);
	static vector<vector<double>> TranspozeHesapla(vector<vector<double>> matris);
	static vector<vector<double>> MartisNormlastir(vector<vector<double>> matris, double norm);
	static void MatrisYazdir(vector<vector<double>> matris);
	static void VektorYazdir(vector<double> vektor);
	static vector<vector<double>> GausMatrisTersi(vector<vector<double>> matris);
	static vector<vector<double>> CholeskyFactorHesapla(vector<vector<double>> matris);
	double HadamardKatsayisiHesapla(vector<vector<double>> matris);
	vector<vector<double>> Inverse(vector<vector<double>> matris);
	vector<vector<double>> CofactorHesapla(const vector<vector<double>> matris);
	vector<double> GaussJordanBilinmeyenVektorHesapla(vector<vector<double>> matris);
};

