#include <iostream>
#include <array>
#include<math.h>
#include <vector>
#include <iomanip>
using namespace std;

const int M = 5;
const int N = 5;

array<double, M> b{ 0.7577, 0.7431, 0.3922, 0.6555, 0.1712 };

array<array<double, M>, N> a = {
	0.8147, 0.0975, 0.1576, 0.1419, 0.6557 ,
	0.9058, 0.2785, 0.9706, 0.4218, 0.0357 ,
	0.1270, 0.5469, 0.9572, 0.9157, 0.8491 ,
	0.9134, 0.9575, 0.4854, 0.7922, 0.9340 ,
	0.6324, 0.9649, 0.8003, 0.9595, 0.6787 };

double determinantHesapla(array<array<double, M>, N> matris, int n = N)
{
	double determinant = 0;
	array<array<double, M>, N> temp{};

	if (n == 2)
		return ((matris[0][0] * matris[1][1]) - (matris[1][0] * matris[0][1]));
	else {
		for (int x = 0; x < n; x++) {
			int subi = 0;
			for (int i = 1; i < n; i++) {
				int subj = 0;
				for (int j = 0; j < n; j++) {
					if (j == x)
						continue;
					temp[subi][subj] = matris[i][j];

					subj++;
				}
				subi++;
			}
			determinant = determinant + (pow(-1, x) * matris[0][x] * determinantHesapla(temp, n - 1));
		}
	}
	return determinant;
}

double izHesapla(array<array<double, M>, N> matris) {
	double sum = 0;

	for (int i = 0; i < N; i++)
		sum += matris[i][i];

	return sum;
}

double satirNorm(array<array<double, M>, N> matris, int satirNo) {
	double sum = 0;

	for (int i = 0; i < N; i++)
		sum += pow(matris[satirNo][i], 2);

	return sqrt(sum);
}

double sutunNorm(array<array<double, M>, N> matris, int sutunNo) {
	double sum = 0;

	for (int i = 0; i < N; i++)
		sum += pow(matris[i][sutunNo], 2);

	return sqrt(sum);
}

double oklidNorm(array<array<double, M>, N> matris) {
	double sum = 0;

	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
			sum += pow(matris[i][j], 2);

	return sqrt(sum);
}

array<array<double, M>, N> matrisCarpimHesapla(array<array<double, M>, N> m1, array<array<double, M>, N> m2) {
	array<array<double, M>, N> carpim{};

	for (int i = 0; i < M; ++i)
		for (int j = 0; j < N; ++j)
			for (int k = 0; k < N; ++k)
				carpim[i][j] += m1[i][k] * m2[k][j];

	return carpim;
}

array<array<double, M>, N> transpozeHesapla(array<array<double, M>, N> matris) {
	array<array<double, M>, N> transpose{};

	for (int i = 0; i < M; ++i)
		for (int j = 0; j < N; ++j) {
			transpose[j][i] = a[i][j];
		}

	return transpose;
}

array<array<double, M>, N> martisNormlastirVeYazdir(array<array<double, M>, N> matris, double norm) {
	array<array<double, M>, N> normlastirilmis;
	
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			normlastirilmis[i][j] = matris[i][j] / norm;
		}
	}

	return normlastirilmis;
}

void matrisYazdir(array<array<double, M>, N> matris) {
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			cout << "\t" << matris[i][j] << "\t";
		}
		cout << "\n";
	}
}


double spektralKondKatsayisi(array<array<double, M>, N> matris) {
	double determinant = abs(determinantHesapla(matris));
	double satirNormlariToplami = 0;

	for (int i = 0; i < M; ++i)
		satirNormlariToplami += satirNorm(matris, i);

	return determinant / satirNormlariToplami;

}

int main()
{
	cout << fixed << setprecision(7);

	cout << "DET(A): " << determinantHesapla(a) << "\n\n";
	cout << "IZ(ATA): " << izHesapla(matrisCarpimHesapla(transpozeHesapla(a), a)) << "\n\n";
	cout << "SATIR N(a1): " << satirNorm(a, 0) << "\n";
	cout << "SATIR N(a2): " << satirNorm(a, 1) << "\n";
	cout << "SATIR N(a3): " << satirNorm(a, 2) << "\n";
	cout << "SATIR N(a4): " << satirNorm(a, 3) << "\n";
	cout << "SATIR N(a5): " << satirNorm(a, 4) << "\n\n";
	cout << "SUTUN N(a1): " << sutunNorm(a, 0) << "\n";
	cout << "SUTUN N(a2): " << sutunNorm(a, 1) << "\n";
	cout << "SUTUN N(a3): " << sutunNorm(a, 2) << "\n";
	cout << "SUTUN N(a4): " << sutunNorm(a, 3) << "\n";
	cout << "SUTUN N(a5): " << sutunNorm(a, 4) << "\n\n";
	cout << "OKLID N(A): " << oklidNorm(a) << "\n\n";
	cout << "IZ(ATA)^1/2: " << sqrt(izHesapla(matrisCarpimHesapla(transpozeHesapla(a), a))) << "     =>    N(A) = IZ(ATA)^1/2" << "\n\n";
	
	cout << "Oklid Normuna Gore A matrisi  \n";
	matrisYazdir(martisNormlastirVeYazdir(a, oklidNorm(a)));
	cout << "\n\n";

	cout << "SPEKTRAL KS(A): " << spektralKondKatsayisi(a) << "\n\n";

	system("pause>0");
}
