#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <exception>
#include <algorithm>

// #define DEBUG_MATRIX
// #define DEBUG_STEP

#define EPS 0.0001

// ”словие окончани€
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
’од метода, где:
a[n][n] - ћатрица коэффициентов
x[n], p[n] - “екущее и предыдущее решени€
b[n] - —толбец правых частей
¬се перечисленные массивы вещественные и
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

bool Gauss(double **A, double *B, double *x, const int number) // ф-и€ решение системы методом √аусса
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
	catch (std::exception e)
	{
		printf("ERROR at Gauss: %s\n", e.what());
		return false;
	}
}

inline double U(double x){
	return x * sin(x);
}

inline double U_x(double x){
	return sin(x) + x * cos(x);
}

void fill_matrix_rhs(double** mt, double* rhs, int nx, int ny, double* u, double* v, double dx, double dy, double dt, double viscosity){
	int nxy = nx * ny;
	for (int i = 0; i < nxy; i++)
	for (int j = 0; j < nxy; j++)
		mt[i][j] = 0;

	for (int i = 0; i < nx; i++){
		// на нижней границе u = 0
		mt[i*ny][i*ny] = 1;
		rhs[i * ny] = 0;

		// на верхней границе u = U(x)
		mt[i*ny + ny - 1][i*ny + ny - 1] = 1;
		rhs[i*ny + ny - 1] = U(i*dx);
	}

	for (int j = 1; j < ny - 1; j++){
		// на левой границе u = U(0)
		mt[j][j] = 1;
		rhs[j] = U(0);

		// на правой границе u_x = 0
		mt[(nx - 1)*ny + j][(nx - 1)*ny + j] = 1;
		mt[(nx - 1)*ny + j][(nx - 2)*ny + j] = -1;
		rhs[(nx - 1)*ny + j] = 0;
	}

	double max_d1 = 0;
	double max_d2 = 0;
	int d1_x = 0, d1_y = 0, d2_x = 0, d2_y = 0;

	double coeff = viscosity / dy / dy;

	for (int i = 1; i < nx - 1; i++)
	for (int j = 1; j < ny - 1; j++){
		int eq_ind = i * ny + j;
		mt[eq_ind][eq_ind] = 1 / dt + 2 * coeff +u[eq_ind] / dx + v[eq_ind] / dy;
		mt[eq_ind][eq_ind + 1] = -coeff;
		mt[eq_ind][eq_ind - 1] = -coeff - v[eq_ind] / dy;
		mt[eq_ind][eq_ind - ny] = -u[eq_ind] / dx;
		//double u_u_x = -u[eq_ind] * (u[eq_ind] - u[eq_ind - ny]) / dx;
		//double v_u_y = -v[eq_ind] * (u[eq_ind] - u[eq_ind - 1]) / dy;
		rhs[eq_ind] = u[eq_ind] / dt + U(i*dx) * U_x(i*dx);

		/*if (abs(u_u_x) > max_d1){
			d1_x = i;
			d1_y = j;
			max_d1 = abs(u_u_x);
		}
		if (abs(v_u_y) > max_d2){
			d2_x = i;
			d2_y = j;
			max_d2 = abs(v_u_y);
		}*/
	}

	// printf("MAX ABS: u_x %lf [%d %d] u_y %lf [%d %d]\n", max_d1, d1_x, d1_y, max_d2, d2_x, d2_y);
}

int main(){

	double width = M_PI/2.0;
	double height = 0.2;
	double total_time = 2;
	int step_num = 100;
	double dt = total_time / step_num;

	int nx = 40;
	int ny = 40;
	int nxy = nx*ny;

	double dx = width / (nx - 1.0);
	double dy = height / (ny - 1.0);

	double viscosity = 0.001;

	double* u = new double[nxy];
	double* v = new double[nxy];

	double* new_u = new double[nxy];
	double* temp_u = new double[nxy];

	double* rhs = new double[nxy];

	for (int i = 0; i < nxy; i++){
		u[i] = 0;
		v[i] = 0;
	}

	double** u_matrix = new double*[nxy];
	for (int i = 0; i < nxy; i++)
		u_matrix[i] = new double[nxy];
	
	double** u_matrix_temp = new double*[nxy];
	for (int i = 0; i < nxy; i++)
		u_matrix_temp[i] = new double[nxy];

	// fill_matrix_rhs(u_matrix, rhs, nx, ny, u, v, dx, dy, dt, viscosity);

	FILE* result_f;
	fopen_s(&result_f, "velocity.txt", "w");

	FILE* velocity_x_f; 
	fopen_s(&velocity_x_f, "velocity_x.txt", "w");

	FILE* velocity_y_f;
	fopen_s(&velocity_y_f, "velocity_y.txt", "w");

#ifdef DEBUG_MATRIX
	FILE* matrix_f;
	matrix_f = fopen("matrix.txt", "w");
	for (int i = 0; i < nxy; i++){
		for (int j = 0; j < nxy; j++)
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

	for (int step = 0; step < step_num; step++){
		printf("Step #%d of %d\n", step+1, step_num);
		fill_matrix_rhs(u_matrix, rhs, nx, ny, u, v, dx, dy, dt, viscosity);

		/* bool code = Gauss(u_matrix, rhs, new_u, nxy);
		if (code == false){
			printf("Gauss ERROR\n");
			return -1;
		}*/

		// задаем первое приближенное решение
		for (int i = 0; i < nxy; i++)
			new_u[i] = rhs[i];

		gauss_zeidel(nxy, u_matrix, rhs, new_u, temp_u);

		for (int i = 0; i < nxy; i++)
			u[i] = new_u[i];

		// считаем значени€ v на новом шаге через уравнение u_x + v_y = 0
		// условие на левой границе v = 0
		for (int j = 0; j < ny; j++)
			v[j] = 0;
		// на верхней и нижней границах v = 0
		for (int i = 0; i < nx; i++){
			v[i*ny] = 0;
			v[i*ny + ny - 1] = 0;
		}

		for (int i = 1; i < nx; i++)
		for (int j = 1; j < ny - 1; j++){
			// на правой границе v_y = 0
			if (i == nx - 1)
				v[i * ny + j] = v[i * ny + j - ny];
			// u_x + v_y = 0
			else
				v[i*ny + j] = v[i*ny + j - 1] + dy / dx * (new_u[i * ny + j] - new_u[(i - 1) * ny + j]);
		}

		if (step % 3 == 0){
			for (int i = 0; i < nx; i++){
				for (int j = 0; j < ny; j++){
					fprintf(result_f, "%lf %lf %lf %lf\n", i*dx, j*dy, 0.2*u[i*ny + j], 0.2*v[i*ny + j] * height / width);
					fprintf(velocity_x_f, "%lf %lf %lf %lf\n", i*dx, j*dy, 0.2*u[i*ny + j], 0.0);
					fprintf(velocity_y_f, "%lf %lf %lf %lf\n", i*dx, j*dy, 0.0, 0.2*v[i*ny + j] * height / width);
				}
			}

			if (step < step_num - 1){
				fprintf(result_f, "\n\n");
				fprintf(velocity_x_f, "\n\n");
				fprintf(velocity_y_f, "\n\n");
			}
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
	fclose(velocity_x_f);
	fclose(velocity_y_f);

#ifdef DEBUG_STEP
	fclose(debug_step_f);
	fclose(rhs_f);
#endif //DEBUG_STEP
}