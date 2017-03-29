#define _USE_MATH_DEFINES
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

const int N = 50;
const int M = 50;
const double k = 20;
const double epsilon = 0.1;
const double Y0 = 0.5;

bool Gauss(vector<vector<double>> & A, vector<double> & B, vector<double> & x, const int & number) // ф-и¤ решение системы методом √аусса
{
    try
    {
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
                            swap(A[i][k], A[j][k]);
                        }
                        swap(B[i], B[j]);
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
                    A[j][k] -= A[i][k] * flag;
                B[j] -= B[i] * flag;
            }
        }

        for (i = 0; i < number; i++) // приведение к виду, где на диагонали сто¤т единицы
        {
            flag = A[i][i];
            for (j = i; j < number; j++)
            {
                A[i][j] /= flag;
            }
            B[i] /= flag;
        }
        fill(x.begin(), x.begin() + number, 0);
        x[number - 1] = B[number - 1];
        for (i = number - 2; i >= 0; i--)
        {
            x[i] = B[i];
            for (j = i + 1; j < number; j++)
            {
                x[i] -= A[i][j] * x[j];
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
    const int NN = 10000;
    double res = 0;
#pragma omp parallel for reduction(+:res)
    for (int i = 0; i <= NN; i++)
        res += sin(i*Y0)*sin(i*y);
    res = res * 2 / (double)NN;
    return res;
}

inline double phi(double y)  //сюда ставим свою функцию
{
    return y*(cosh(0.5) - cosh(y - 0.5));
}

void printanswer(vector<double> & A, int N, int M, FILE *pFile)
{
    int I = N * 1;
    for (int i = 0; i < M; i++)
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
    vector<vector<double>> matrix(M*N);
    vector<double> U(M*N), F(M*N);

    X = 5;
    Y = 1;

    delta_x = X / N;
    delta_y = Y / M;

    coef1 = 1. / delta_x / delta_x;
    coef2 = 1. / delta_y / delta_y;

    for (int i = 0; i < M*N; i++)
    {
        matrix[i] = vector<double>(M*N);
        matrix[i][i] = 1;
    }

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

    fill(F.begin(), F.end(), 0);

    for (int i = 0; i < M; i++)
        F[i] = func(i * delta_y);

    Gauss(matrix, F, U, M * N);
    for (int i = 0; i < M; i++)
        cout << U[i*M] << "   " << U[(i + 1)*N - 1] << endl;
    fopen_s(&pFile, "answer.txt", "w");
    printanswer(U, M, N, pFile);
    fclose(pFile);
    system("plot_script.plt");
    return 0;
}