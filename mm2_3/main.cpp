#define _USE_MATH_DEFINES
#include <cstdio>
#include <cmath>
#include <exception>
#include <algorithm>
#include <vector>
using namespace std;

// #define DEBUG_MATRIX
// #define DEBUG_STEP

#define EPS 0.0001

// ������� ���������
bool converge(vector<double> & xk, vector<double> & xkp, int n)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
    {
        norm += (xk[i] - xkp[i])*(xk[i] - xkp[i]);
    }
    if (sqrt(norm) >= EPS)
        return false;
    return true;
}

/*
��� ������, ���:
a[n][n] - ������� �������������
x[n], p[n] - ������� � ���������� �������
b[n] - ������� ������ ������
��� ������������� ������� ������������ �
������ ���� ���������� � �������� ���������,
����� � ������ x[n] ������� ��������� ���������
����������� ������� ������� (��������, ��� ����)
*/
void gauss_zeidel(int n, vector<vector<double>> & a, vector<double> &  b, vector<double> & x, vector<double> & p)
{
    do
    {
        copy(x.begin(), x.end(), p.begin());  
        for (int i = 0; i < n; i++)
        {
            double var = 0;
            for (int j = 0; j < i; j++)
                var += (a[i][j] * x[j]);
            for (int j = i + 1; j < n; j++)
                var += (a[i][j] * p[j]);
            x[i] = (b[i] - var) / a[i][i];
        }
    } while (!converge(x, p, n));
}

inline double phi(double x)
{
    return sin(x) * exp(-x);
}

void fill_rhs(vector<double> & rhs, int nx, int y, vector<double> & u, double dy)
{
    rhs[0] = 0;
    rhs[nx - 1] = 0;
    for (int i = 1; i < nx - 1; i++)
    {
        rhs[i] = (2 * u[(y - 1)*nx + i] - u[(y - 2)*nx + i]) / dy / dy;
    }
}

int main()
{
    double width = 2;
    double height = 20;

    int nx = 100; // etta
    int ny = 200; // tau
    int nxy = nx*ny;

    double dx = width / (nx - 1.0);
    double dy = height / (ny - 1.0);

    // ���� ����������� ���� �����������
    double k = 5;

    // �������� ������ ������� - �����? ����� :3
    vector<double> u(nxy),new_u(nx),temp_u(nx),rhs(nx);
    vector<vector<double>> u_matrix(nx);
 
    for (int i = 0; i < nx; i++)
        u_matrix[i] = vector<double>(nx);

    FILE* result_f, *animation_f;
    fopen_s(&result_f, "velocity.txt", "w");
    fopen_s(&animation_f, "animation_data.txt", "w");

    // �� ������� ��������, ��� u(x, 0) = phi � u_y(x, 0) = 0
    // ������ ������ ���� �� y ����� phi(x)
    for (int i = 0; i < nx; i++)
    {
        u[i] = u[i + nx] = phi(i * dx);
    }

    // ���������� ��������� ��������� � ����
    for (int i = 0; i < nx; i++)
    {
        fprintf(animation_f, "%lf %lf\n", i*dx, u[i]);
    }
    fprintf(animation_f, "\n\n");

    for (int i = 0; i < nx; i++)
    {
        fprintf(result_f, "%lf ", u[i]);
    }
    fprintf(result_f, "\n");


    //���� ��� ��������� ������� �������� ���������� ����� 
    // * * *
    //   *
    //   *  
    u_matrix[0][0] = 1;
    u_matrix[nx - 1][nx - 1] = 1;
    for (int i = 1; i < nx - 1; i++)
    {
        u_matrix[i][i] = 1 / dy / dy + 2 / dx / dx + k;
        u_matrix[i][i - 1] = -1 / dx / dx;
        u_matrix[i][i + 1] = -1 / dx / dx;
    }

#ifdef DEBUG_MATRIX
    FILE* matrix_f;
    fopen_s(&matrix_f, "matrix.txt", "w");
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nx; j++)
            fprintf(matrix_f, "%.10f ", u_matrix[i][j]);
        fprintf(matrix_f, "\n");
    }
    fclose(matrix_f);
#endif // DEBUG_MATRIX

#ifdef DEBUG_STEP
    FILE* rhs_f;
    fopen_s(&rhs_f, "rhs.txt", "w");

    FILE* debug_step_f;
    fopen_s(&debug_step_f, "step.txt", "w");
#endif // DEBUG_STEP

    for (int y = 2; y < ny; y++)
    {
        printf("Step #%d of %d\n", y + 1, ny);

        //��������� ������ ����� (������� �� ��������)
        fill_rhs(rhs, nx, y, u, dy);

        // ������ ������ ������������ �������
        copy(rhs.begin(), rhs.end(), new_u.begin());  

        gauss_zeidel(nx, u_matrix, rhs, new_u, temp_u);

        // �������� ���������� ������� � �������
        for (int i = 0; i < nx; i++)
            u[y*nx + i] = new_u[i];

        for (int i = 0; i < nx; i++)
        {
            fprintf(animation_f, "%lf %lf\n", i*dx, u[y*nx + i]);
        }
        fprintf(animation_f, "\n\n");

        if (y % 2 == 0)
        {
            for (int i = 0; i < nx; i++)
            {
                fprintf(result_f, "%lf ", u[y*nx + i]);
            }
            fprintf(result_f, "\n");
        }

#ifdef DEBUG_STEP
        fprintf(rhs_f, "RHS step %i\n", step);
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
                fprintf(rhs_f, "%lf ", rhs[i*ny + j]);
            fprintf(rhs_f, "\n");
        }
        fprintf(rhs_f, "\n");

        fprintf(debug_step_f, "Step %i\nU\n", step);
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
                fprintf(debug_step_f, "%lf ", u[i*ny + j]);
            fprintf(debug_step_f, "\n");
    }
        fprintf(debug_step_f, "\nV\n");
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
                fprintf(debug_step_f, "%lf ", v[i*ny  j]);
            fprintf(debug_step_f, "\n");
        }
        fprintf(debug_step_f, "\n");
#endif // DEBUG_STEP
    }

    fclose(result_f);
    fclose(animation_f);   
#ifdef DEBUG_STEP
    fclose(debug_step_f);
    fclose(rhs_f);
#endif //DEBUG_STEP
	system("plot_script.plt");
	system("animation_script.plt");
	system("profile.gif");
}