#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <exception>
#include <algorithm>

// #define DEBUG_MATRIX
// #define DEBUG_STEP

#define EPS 0.0001

// Условие окончания
bool converge(double *xk, double *xkp, int n)
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
Ход метода, где:
a[n][n] - Матрица коэффициентов
x[n], p[n] - Текущее и предыдущее решения
b[n] - Столбец правых частей
Все перечисленные массивы вещественные и
должны быть определены в основной программе,
также в массив x[n] следует поместить начальное
приближение столбца решений (например, все нули)
*/
void gauss_zeidel(int n, double** a, double* b, double* x, double* p){
	do
	{
		for (int i = 0; i < n; i++)
			p[i] = x[i];

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

void print_gnuplot_script(int ny){

}

bool Gauss(double **A, double *B, double *x, const int number) // ф-ия решение системы методом Гаусса
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
						printf("It isn't homogeneous system\n");
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

		for (i = 0; i < number; i++) // приведение к виду, где на диагонали стоят единицы
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
	catch (std::exception e)
	{
		printf("ERROR at Gauss: %s\n", e.what());
		return false;
	}
}

inline double phi(double x){
	return x*x/exp(x);
}

void fill_rhs(double* rhs, int nx, int y, double* u, double dy){
	rhs[0] = 0;
	rhs[nx - 1] = 0;
	for (int i = 1; i < nx - 1; i++){
		rhs[i] = (2 * u[(y - 1)*nx + i] - u[(y - 2)*nx + i]) / dy / dy;
	}
}

int main(){

	double width = 2;
	double height = 20;

	int nx = 100; // etta
	int ny = 200; // tau
	int nxy = nx*ny;

	double dx = width / (nx - 1.0);
	double dy = height / (ny - 1.0);

	double k = 1;

	double* u = new double[nxy];

	double* new_u = new double[nx];
	double* temp_u = new double[nx];

	double* rhs = new double[nx];

	for (int i = 0; i < nxy; i++)
		u[i] = 0;

	double** u_matrix = new double*[nx];
	for (int i = 0; i < nx; i++)
		u_matrix[i] = new double[nx];

	FILE* result_f, *animation_f;
	fopen_s(&result_f, "velocity.txt", "w");
	fopen_s(&animation_f, "animation_data.txt", "w");

	// по условию известно, что u(x, 0) = phi и u_y(x, 0) = 0
	// значит первые слоя по y равны phi(x)
	for (int i = 0; i < nx; i++){
		u[i] = u[i + nx] = phi(i * dx);
		printf("[%lf %lf] ", i*dx, phi(i*dx));
	}

	for (int i = 0; i < nx; i++){
		fprintf(animation_f, "%lf %lf\n", i*dx, u[i]);
	}
	fprintf(animation_f, "\n\n");

	for (int i = 0; i < nx; i++){
		fprintf(result_f, "%lf ", u[i]);
	}
	fprintf(result_f, "\n");

	for (int i = 0; i < nx; i++)
	for (int j = 0; j < nx; j++)
		u_matrix[i][j] = 0;

	u_matrix[0][0] = 1;
	u_matrix[nx - 1][nx - 1] = 1;
	for (int i = 1; i < nx - 1; i++){
		u_matrix[i][i] = 1 / dy / dy + 2 / dx / dx - k;
		u_matrix[i][i - 1] = - 1 / dx / dx;
		u_matrix[i][i + 1] = - 1 / dx /dx;
	}

#ifdef DEBUG_MATRIX
	FILE* matrix_f;
	fopen_s(&matrix_f, "matrix.txt", "w");
	for (int i = 0; i < nx; i++){
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

	for (int y = 2; y < ny; y++){
		printf("Step #%d of %d\n", y + 1, ny);
		fill_rhs(rhs, nx, y, u, dy);

		// задаем первое приближенное решение
		for (int i = 0; i < nx; i++)
			new_u[i] = rhs[i];

		gauss_zeidel(nx, u_matrix, rhs, new_u, temp_u);

		for (int i = 0; i < nx; i++)
			u[y*nx + i] = new_u[i];

		for (int i = 0; i < nx; i++){
			fprintf(animation_f, "%lf %lf\n", i*dx, u[y*nx + i]);
		}
		fprintf(animation_f, "\n\n");

		if (y % 2 == 0){
			for (int i = 0; i < nx; i++){
				fprintf(result_f, "%lf ", u[y*nx + i]);
			}
			fprintf(result_f, "\n");
		}

#ifdef DEBUG_STEP
		fprintf(rhs_f, "RHS step %i\n", step);
		for (int i = 0; i < nx; i++){
			for (int j = 0; j < ny; j++)
				fprintf(rhs_f, "%lf ", rhs[i*ny + j]);
			fprintf(rhs_f, "\n");
		}
		fprintf(rhs_f, "\n");

		fprintf(debug_step_f, "Step %i\nU\n", step);
		for (int i = 0; i < nx; i++){
			for (int j = 0; j < ny; j++)
				fprintf(debug_step_f, "%lf ", u[i*ny + j]);
			fprintf(debug_step_f, "\n");
		}
		fprintf(debug_step_f, "\nV\n");
		for (int i = 0; i < nx; i++){
			for (int j = 0; j < ny; j++)
				fprintf(debug_step_f, "%lf ", v[i*ny + j]);
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
}