#include <iostream>
#include <array>
#include<math.h>
#include <vector>
#include <iomanip>
#include "MatrixReader.h"
#include "OmerFarukOdev.h"
#include<algorithm>
#include<locale.h>

using namespace std;
//
//const int M = 5;
//const int N = 5;
//
//array<double, M> b{ 0.7577, 0.7431, 0.3922, 0.6555, 0.1712 };
//
//array<array<double, M>, N> a = {
//	0.8147, 0.0975, 0.1576, 0.1419, 0.6557 ,
//	0.9058, 0.2785, 0.9706, 0.4218, 0.0357 ,
//	0.1270, 0.5469, 0.9572, 0.9157, 0.8491 ,
//	0.9134, 0.9575, 0.4854, 0.7922, 0.9340 ,
//	0.6324, 0.9649, 0.8003, 0.9595, 0.6787 };


double determinantHesapla(vector<vector<double>> matris, int n = 0)
{
	vector<vector<double>> temp(matris.size(), vector<double>(matris[0].size()));
	double determinant = 0;

	if (n == 0)
		n = matris.size();

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

double izHesapla(vector<vector<double>> matris) {
	double sum = 0;
	double n = matris.size();

	for (int i = 0; i < n; i++)
		sum += matris[i][i];

	return sum;
}

double satirNorm(vector<vector<double>> matris, int satirNo) {
	double sum = 0;
	int n = matris[0].size();

	for (int i = 0; i < n; i++)
		sum += pow(matris[satirNo][i], 2);

	return sqrt(sum);
}

double sutunNorm(vector<vector<double>> matris, int sutunNo) {
	double sum = 0;
	int n = matris.size();

	for (int i = 0; i < n; i++)
		sum += pow(matris[i][sutunNo], 2);

	return sqrt(sum);
}

double oklidNorm(vector<vector<double>> matris) {
	double sum = 0;
	int m = matris.size();
	int n = matris[0].size();

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			sum += pow(matris[i][j], 2);

	return sqrt(sum);
}

vector<vector<double>> matrisCarpimHesapla(vector<vector<double>> m1, vector<vector<double>> m2) {

	vector<vector<double>> carpim(m1.size(), vector<double>(m1[0].size()));
	int m = m1.size();
	int n = m1[0].size();

	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			for (int k = 0; k < n; ++k)
				carpim[i][j] += m1[i][k] * m2[k][j];

	return carpim;
}

vector<vector<double>> transpozeHesapla(vector<vector<double>> matris) {

	vector<vector<double>> transpose(matris.size(), vector<double>(matris[0].size()));
	int m = matris.size();
	int n = matris[0].size();

	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j) {
			transpose[j][i] = matris[i][j];
		}

	return transpose;
}

vector<vector<double>> martisNormlastirVeYazdir(vector<vector<double>> matris, double norm) {

	vector<vector<double>> normlastirilmis(matris.size(), vector<double>(matris[0].size()));
	int m = matris.size();
	int n = matris[0].size();

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			normlastirilmis[i][j] = matris[i][j] / norm;
		}
	}

	return normlastirilmis;
}

void matrisYazdir(vector<vector<double>> matris) {
	int m = matris.size();
	int n = matris[0].size();

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			cout << "\t" << matris[i][j] << "\t";
		}
		cout << "\n";
	}
}


double spektralKondKatsayisi(vector<vector<double>> matris) {
	double determinant = abs(determinantHesapla(matris));
	double satirNormlariToplami = 0;
	int m = matris.size();

	for (int i = 0; i < m; ++i)
		satirNormlariToplami += satirNorm(matris, i);

	return determinant / satirNormlariToplami;

}

int main()
{
	setlocale(LC_ALL, "Turkish");
	cout << fixed << setprecision(7);

	string matrisFilePath = "matrix5.txt";
	string vektorFilePath = "vector5.txt";
	int menu;

	//cout << "Matris dosya yolunu giriniz: ";
	//cin >> matrisFilePath;

	//cout << "Vector dosya yolunu giriniz: ";
	//cin >> vektorFilePath;

	MatrixReader ornek;
	vector<vector<double>> matris = ornek.ReadMatris(matrisFilePath);
	vector<double> vektor = ornek.ReadVector(vektorFilePath);
	int m = matris.size();


	cout << "####### Aþaðýdaki menüden bir matris iþlemi seçiniz #######\n";
	cout << "\t 1. A Determinant hesapla. \n";
	cout << "\t 2. ATA Matris izini hesapla. \n";
	cout << "\t 3. A Matris satir normlarini hesapla. \n";
	cout << "\t 4. A Matris sütun normlarini hesapla. \n";
	cout << "\t 5. A Matris Oklid normlarini hesapla. \n";
	cout << "\t 6. A N(A) = (iz(ATA))^1/2 esit mi. \n";
	cout << "\t 7. A Matrisi Oklid normuna gore normlastir. \n";
	cout << "\t 8. A Matrisin ozdegerlerini hesapla. \n";
	cout << "\t 9. A Matrisin Spektral (Todd) sart sayisini hesapla ve kararsizligini yorumla. \n";
	cout << "\t 10. A Matrissinin Hadamard sart sayisini hesapla ve kararsizligini yorumla. \n";
	cout << "\t 11. Kramer kurali ile A Matrisinin tersini hesapla. \n";
	cout << "\t 12. Pivotlama ile matris tersi. \n";
	cout << "\t 13. Gauss ile matris tersi. \n";
	cout << "\t 14. Gauss algoritmasi ile x bilinmeyenler verktorunu hesapla. \n";
	cout << "\t 15. Gauss-Jordan yontemi ile x bilinmeyenler verktorunu hesapla. \n";
	cout << "\t 16. Modernlestirilmis Gauss Algoritmasi ile x bilinmeyenler vektorunu hesapla. \n";
	cout << "\t 17. Cholesky yontemi ile x bilinmeyenler vektorunu hesapla. \n";
	cout << "\t 18. Modernlestirilmis Gauss Algoritmasi ile ATA matrisinin tersini hesapla. \n";
	cout << "\t 19. Cholesky yontemi ile ATA matrisinin tersini hesapla. \n";

	while (true)
	{

	cout << "Menü seçiminiz ";
	cin >> menu;



		if (menu == 1)
			cout << "DET(A): " << determinantHesapla(matris) << "\n\n";

		else if (menu == 2)
			cout << "IZ(ATA): " << izHesapla(matrisCarpimHesapla(transpozeHesapla(matris), matris)) << "\n\n";

		else if (menu == 3) {
			for (int i = 0; i < m; i++)
			{
				cout << "SATIR N(a" << i << "): " << satirNorm(matris, i) << "\n";
			}
			cout << "\n";
		}

		else if (menu == 4) {
			for (int i = 0; i < m; i++)
			{
				cout << "SUTUN N(a" << i << "): " << sutunNorm(matris, i) << "\n";
			}

			cout << "\n";
		}

		else if (menu == 5)
			cout << "OKLID N(A): " << oklidNorm(matris) << "\n\n";

		else if (menu == 6)
			cout << "IZ(ATA)^1/2: " << sqrt(izHesapla(matrisCarpimHesapla(transpozeHesapla(matris), matris))) << "     =>    N(A) = IZ(ATA)^1/2" << "\n\n";

		else if (menu == 7) {
			cout << "Oklid Normuna Gore A matrisi  \n";
			matrisYazdir(martisNormlastirVeYazdir(matris, oklidNorm(matris)));
			cout << "\n\n";
		}

		else if (menu == 9)
			cout << "Hesaplama yapýlamadý " << "\n";

		else if (menu == 9)
			cout << "SPEKTRAL KS(A): " << spektralKondKatsayisi(matris) << "\n\n";

		else
		{
			cout << "Geçersiz bir iþlem. " << "\n";
		}

	}
	system("pause>0");
}
