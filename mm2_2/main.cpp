#define _USE_MATH_DEFINES
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
using namespace std;

const int N = 50;
const int M = 50;
const double k = 20;
const double epsilon = 0.1;
const double Y0 = 0.5;

bool Gauss(double **Aa, double *Bb, double *x, const int & number) // ф-и€ решение системы методом √аусса
{
	try
	{
		double **A = (double**)malloc(number * sizeof(double*));
		for (int i = 0; i < number; i++)
			A[i] = (double*)malloc(number * sizeof(double));

		for (int i = 0; i < number; i++)
		for (int j = 0; j < number; j++)
			A[i][j] = Aa[i][j];

		double *B = (double*)malloc(number * sizeof(double));
		for (int i = 0; i < number; i++)
			B[i] = Bb[i];

		int i, j, k;
		double flag;
		bool b;
		bool NotHomogeneous = true;
		for (i = 0; i < number - 1; i++)
		{
			if (A[i][i] == 0) // прохождение по диагональным элементам c проверкой равенства нулю
			{
				j = i + 1;
				b = false;
				while (b == false)
				{
					if (!(A[j][i] == 0))
					{
						for (k = 0; k < number; k++)
						{
							flag = A[i][k];
							A[i][k] = A[j][k];
							A[j][k] = flag;
						}
						flag = B[i];
						B[i] = B[j];
						B[j] = flag;
						b = true;
					}
					j++;
					if (j > number)
					{
						cout << "It isn't homogeneous system" << endl;
						NotHomogeneous = false;
						return NotHomogeneous;
					}
				}
			}
			for (j = i + 1; j < number; j++) // j - номер строки
			{
				flag = A[j][i] / A[i][i];
				for (k = i; k < number; k++)
					A[j][k] = A[j][k] - A[i][k] * flag;
				B[j] = B[j] - B[i] * flag;
			}
		}

		for (i = 0; i < number; i++) // приведение к виду, где на диагонали сто€т единицы
		{
			flag = A[i][i];
			for (j = i; j < number; j++)
			{
				A[i][j] = A[i][j] / flag;
			}
			B[i] = B[i] / flag;
		}

		for (i = 0; i < number; i++)
			x[i] = 0;
		x[number - 1] = B[number - 1];
		for (i = number - 2; i >= 0; i--)
		{
			x[i] = B[i];
			for (j = i + 1; j < number; j++)
			{
				x[i] = x[i] - A[i][j] * x[j];
			}
		}
		return NotHomogeneous;
	}
	catch (exception& e)
	{
		cout << "ERROR at Gauss\n";
		cout << e.what() << endl;
		return false;
	}
}


double func(double y)
{
	int NN = 10000;
	double res = 0;
	for (int i = 0; i <= NN; i++)
		res += sin(i*Y0)*sin(i*y);
	res = res * 2 / (double)NN;
	return res;
}

 double phi(double y)  //сюда ставим свою функцию
{
	double res;
	res = y * (1 - y)*(0.5 - y)*(0.5 - y);
	//res = y*(cosh(0.5) - cosh(y - 0.5));
	return res;
}

void printanswer(double *A, int N, int M, FILE *pFile)
{
	int I = N * 1;
	for (int i = 0; i < M ; i++) 
	{
		for (int j = 1; j < N; j++)
		{
			fprintf(pFile, "%1.6lf	", A[I]);
			I++;
		}
		fprintf(pFile, "\n");
	}
	fprintf(pFile, "\n\n");
}


int main()
{
	FILE *pFile;

	double X, Y, delta_x, delta_y, coef1, coef2;
	double **matrix, *U, *F;

	X = 5;
	Y = 1;

	delta_x = X / N;
	delta_y = Y / M;

	coef1 = 1. / delta_x / delta_x;
	coef2 = 1. / delta_y / delta_y;

	matrix = (double **)malloc(sizeof(double *)* M * N);
	for (int i = 0; i < M * N; i++)
		matrix[i] = (double *)malloc(sizeof(double)* M * N);
	U = (double *)malloc(sizeof(double)* M * N);
	F = (double *)malloc(sizeof(double)* M * N);

	for (int i = 0; i < M * N; i++)
	for (int j = 0; j < M * N; j++)
		matrix[i][j] = 0;
	
	for (int i = 0; i < N*M; i++)
		matrix[i][i] = 1;

	for (int i = M; i < M * (N - 1); i++)
	{
		if ((i%N == 0) || ((i + 1) % N == 0))
		{
			matrix[i][i] = 1;
			matrix[i][i - 1] = matrix[i][i + 1] = 0;
			matrix[i][i - M] = matrix[i][i + M] = 0;
		}
		else
		{
			matrix[i][i] = -2. / delta_x / delta_x + -2. / delta_y / delta_y + k * k * (1 + epsilon * phi(i * delta_y));
			matrix[i][i - 1] = matrix[i][i + 1] = coef2;
			matrix[i][i - M] = matrix[i][i + M] = coef1;
		}
		
	}
	for (int i = M * (N - 1); i < M * N; i++)
		matrix[i][i - 1] = -1;

	for (int i = 0; i < N * M; i++)
		F[i] = 0;

	for (int i = 0; i < M; i++)
		F[i] = func(i * delta_y);

	Gauss(matrix, F, U, M * N);
	for (int i = 0; i < M; i++) cout << U[i*M] << "   " << U[(i + 1)*N - 1] << endl;
	fopen_s(&pFile, "answer.txt", "w");
	printanswer(U, M, N, pFile);
	fclose(pFile);
	system("plot_script.plt");
	return 0;
}