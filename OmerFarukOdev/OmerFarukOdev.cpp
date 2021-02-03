#include <iostream>
#include <vector>
#include <iomanip>
#include "MatrixReader.h"
#include "OmerFarukOdev.h"

#include <algorithm>

#include "MatrixCalculator.h"

using namespace std;

static void PrintMenu()
{
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
	cout << "\t 17. CholeskyFactorHesapla yontemi ile x bilinmeyenler vektorunu hesapla. \n";
	cout << "\t 18. Modernlestirilmis Gauss Algoritmasi ile ATA matrisinin tersini hesapla. \n";
	cout << "\t 19. CholeskyFactorHesapla yontemi ile ATA matrisinin tersini hesapla. \n";
}

int main()
{
	setlocale(LC_ALL, "Turkish");
	cout << fixed << setprecision(7);

	string matrixFilePath = "matrix5.txt";
	string vectorFilePath = "vector5.txt";
	int menu;

	//cout << "Matris dosya yolunu giriniz: ";
	//cin >> matrisFilePath;

	//cout << "Vector dosya yolunu giriniz: ";
	//cin >> vektorFilePath;

	MatrixReader matrix_reader;
	MatrixCalculator matrix_calculator;

	auto matrix = matrix_reader.ReadMatris(matrixFilePath);
	auto vektor = matrix_reader.ReadVector(vectorFilePath);
	const auto size = static_cast<int>(matrix.size());
	
	PrintMenu();

	while (true)
	{
		cout << endl <<  "Ýþlem seçiminiz: ";
		cin >> menu;

		if (menu == 1)
		{
			cout << "DET(A): " << matrix_calculator.DeterminantHesapla(matrix);
		}
		else if (menu == 2)
		{
			auto t = matrix_calculator.TranspozeHesapla(matrix);
			auto c = matrix_calculator.MatrisCarpimHesapla(t, matrix);

			cout << "IZ(ATA): " << matrix_calculator.IzHesapla(c);
		}
		else if (menu == 3) {
			for (auto i = 0; i < size; i++)
			{
				auto satir_norm = matrix_calculator.SatirNorm(matrix, i);
				cout << "SATIR N(a" << i << "): " << satir_norm;
				cout << endl;
			}
		}
		else if (menu == 4) {
			for (auto i = 0; i < size; i++)
			{
				auto sutun_norm = matrix_calculator.SutunNorm(matrix, i);
				cout << "SUTUN N(a" << i << "): " << sutun_norm;
				cout << endl;
			}
		}
		else if (menu == 5)
		{
			cout << "OKLID N(A): " << matrix_calculator.OklidNorm(matrix);
		}
		else if (menu == 6)
		{
			auto na = matrix_calculator.OklidNorm(matrix);
			auto t = matrix_calculator.TranspozeHesapla(matrix);
			auto m = matrix_calculator.MatrisCarpimHesapla(t, matrix);
			auto iz = matrix_calculator.IzHesapla(m);
			auto izATA = sqrt(iz);

			cout << "IZ(ATA)^1/2 = " << izATA << "   ve   N(A) = " << na << " olduðu için IZ(ATA)^1/2 = N(A) olur.";
		}
		else if (menu == 7) {
			cout << "Oklid Normuna Gore A matrisi  \n";
			auto oklidNorm = matrix_calculator.OklidNorm(matrix);
			auto normlastirilmis = matrix_calculator.MartisNormlastir(matrix, oklidNorm);
			matrix_calculator.MatrisYazdir(normlastirilmis);
		}
		else if (menu == 8)
		{
			auto ozdegerler = matrix_calculator.OzdegerleriHesapla(matrix);
			auto s = static_cast<int>(ozdegerler.size());
			
			cout << "Özdeðerler:";

			cout << " [ ";

			for (auto  i = 0; i< s; i++)
				cout << ozdegerler[i] << " ";

			cout << " ] ";
		}
		else if (menu == 9)
		{
			auto condA = matrix_calculator.SpektralKondKatsayisi(matrix);
			string karar = condA > pow(10, 3) ? "Kararsýz" : "Kararlý";
			string karakter = condA > pow(10, 3) ? ">" : "<";
			
			cout << "cond(A) = " << condA  << " ve cond(A) " << karakter << " 10^3 olduðundan " << " matris " << karar << endl;
		}
		else if (menu == 10)
		{
			auto condA = matrix_calculator.HadamardKatsayisiHesapla(matrix);
			string karar = condA < pow(10,-2) ? "Kararsýz" : "Kararlý";
			string karakter = condA > pow(10, 3) ? ">" : "<";
			
			cout << "Hadamard katsayýsý: " << condA << " ve cond(A) " << karakter<< endl;
			cout << "Matris " << karar << ".";
		}
		else if (menu == 11)
		{
			cout << "Hesaplama yapýlamadý!!!";
		}
		else if (menu == 12)
		{
			cout << "Hesaplama yapýlamadý!!!";
		}
		else if (menu == 13)
		{
			auto sonuc = matrix_calculator.GausMatrisTersi(matrix);

			matrix_calculator.MatrisYazdir(sonuc);

			auto m = matrix_calculator.Inverse(matrix);
			cout << endl;
			matrix_calculator.MatrisYazdir(m);
		}
		else if (menu == 14)
		{
			cout << "Hesaplama yapýlamadý!!!";
		}
		else if (menu == 15)
		{
			cout << "Hesaplama yapýlamadý!!!";
		}
		else if (menu == 16)
		{
			cout << "Hesaplama yapýlamadý!!!";
		}
		else if (menu == 17)
		{
			cout << "Hesaplama yapýlamadý!!!";
		}
		else if (menu == 18)
		{
			cout << "Hesaplama yapýlamadý!!!";
		}
		else if (menu == 19)
		{
			auto a = matrix;
			auto at = matrix_calculator.TranspozeHesapla(a);
			auto ata = matrix_calculator.MatrisCarpimHesapla(a, at);
			
			auto c = matrix_calculator.CholeskyFactorHesapla(ata);			
			auto ic = matrix_calculator.GausMatrisTersi(c);
			auto tc = matrix_calculator.TranspozeHesapla(c);
			auto tic = matrix_calculator.GausMatrisTersi(tc);
			auto r = matrix_calculator.MatrisCarpimHesapla(ic, tic);

			cout << endl << "ATA Matrisi:" << endl;
			matrix_calculator.MatrisYazdir(ata);
			
			cout << endl << "C Matrisi:" << endl;
			matrix_calculator.MatrisYazdir(c);
			
			cout << endl << "A' Matrisi" << endl;			
			matrix_calculator.MatrisYazdir(r);			
		}
		else
		{
			cout << "Geçersiz bir iþlem seçtiniz!!!";
		}
		
		cout << endl;
	}
}
