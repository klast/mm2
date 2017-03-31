#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <exception>
#include <algorithm>

// #define DEBUG_MATRIX
// #define DEBUG_STEP

using namespace std;

#define EPS 0.001
#define SPEED_EPS 0.01
// ����� ������ ������� ���� � ���������

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
��� ������, ���:
a[n][n] - ������� �������������
x[n], p[n] - ������� � ���������� �������
b[n] - ������� ������ ������
��� ������������� ������� ������������ �
������ ���� ���������� � �������� ���������,
����� � ������ x[n] ������� ��������� ���������
����������� ������� ������� (��������, ��� ����)
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

inline double U(double x){
	return x * sin(x);
}

inline double U_x(double x){
	return sin(x) + x * cos(x);
}

/*// ������� �� ����� ������� v = 0
for (int j = 0; j < ny; j++)
	v[j] = 0;
// �� ������� � ������ �������� v = 0
for (int i = 0; i < nx; i++){
	v[i*ny] = 0;
	v[i*ny + ny - 1] = 0;
}

for (int i = 1; i < nx; i++)
	for (int j = 1; j < ny - 1; j++){
		// �� ������ ������� v_y = 0
		if (i == nx - 1)
			v[i * ny + j] = v[i * ny + j - ny];
		// u_x + v_y = 0
		else
			v[i*ny + j] = v[i*ny + j - 1] + dy / dx * (u[i * ny + j] - u[(i - 1) * ny + j]);
	}*/

void fill_matrix_rhs(double** mt, double* rhs, int nx, int ny, double* u, double* v, double dx, double dy, double dt, double viscosity){
	int nxy = nx * ny;
	int eq_num = 2 * nxy;
	int v_offset = 1;
	for (int i = 0; i < eq_num; i++)
	for (int j = 0; j < eq_num; j++)
		mt[i][j] = 0;

	for (int i = 0; i < nx; i++){
		// �� ������ ������� u = 0
		mt[i*ny][i*ny] = 1;
		rhs[i * ny] = 0;

		// �� ������� ������� u = U(x)
		mt[i*ny + ny - 1][i*ny + ny - 1] = 1;
		rhs[i*ny + ny - 1] = U(i*dx);
	}

	for (int j = 1; j < ny - 1; j++){
		// �� ����� ������� u = U(0)
		mt[j][j] = 1;
		rhs[j] = U(0);

		// �� ������ ������� u_x = 0
		mt[(nx - 1)*ny + j][(nx - 1)*ny + j] = 1;
		mt[(nx - 1)*ny + j][(nx - 2)*ny + j] = -1;
		rhs[(nx - 1)*ny + j] = 0;
	}

	double coeff = viscosity / dy / dy;

	// ������������� ������� ��������� ������� 
	for (int i = 1; i < nx - 1; i++)
	for (int j = 1; j < ny - 1; j++){
		int eq_ind = i * ny + j;
		mt[eq_ind][eq_ind] = 1 / dt + 2 * coeff + (u[eq_ind] - u[eq_ind - ny]) / dx;
		mt[eq_ind][eq_ind + 1] = -coeff;
		mt[eq_ind][eq_ind - 1] = -coeff;
		mt[eq_ind][eq_ind + nxy] = (u[eq_ind] - u[eq_ind - 1]) / dy;
		rhs[eq_ind] = u[eq_ind] / dt + U(i*dx) * U_x(i*dx);
	}

	// ������� �� ����� ������� v = 0
	for (int j = 0; j < ny; j++){
		mt[nxy + j][nxy + j] = 1;
		rhs[nxy + j] = 0;
	}
	// �� ������� � ������ �������� v = 0
	for (int i = 0; i < nx; i++){
		mt[nxy + i*ny][nxy + i*ny] = 1;
		rhs[nxy + i*ny] = 0;

		mt[nxy + i*ny + ny - 1][nxy + i*ny + ny - 1] = 1;
		rhs[nxy + i*ny + ny - 1] = 0;
	}
	// �� ������ ������� v_x = 0
	for (int j = 1; j < ny - 1; j++){
		int eq_ind = nxy + (nx - 1) * ny + j;
		mt[eq_ind][eq_ind] = 1 / dx;
		mt[eq_ind][eq_ind - ny] = -1 / dx;
		rhs[eq_ind] = 0;
	}

	for (int i = 1; i < nx - 1; i++)
	for (int j = 1; j < ny - 1; j++){
		int eq_ind = i * ny + j + nxy;
		mt[eq_ind][eq_ind] = 1 / dy;
		mt[eq_ind][eq_ind - 1] = - 1 / dy;
		mt[eq_ind][eq_ind - nxy] = 1 / dx;
		mt[eq_ind][eq_ind - nxy - ny] = - 1 / dy;
		rhs[eq_ind] = 0;
	}
}

int main(){
	double width = M_PI/2.0;
	double height = 1.5;
	double total_time = 0.5;
	int step_num = 100;
	double dt = total_time / step_num;

	int nx = 41;
	int ny = 61;
	int nxy = nx*ny;
	int eq_num = nxy * 2;

	double dx = width / (nx - 1.0);
	double dy = height / (ny - 1.0);

	double viscosity = 0.001;

	double* uv = new double[eq_num];
	double* temp_uv = new double[eq_num];

	double* u = uv;
	double* v = &uv[nxy];

	double* rhs = new double[2*nxy];

	for (int i = 0; i < nxy; i++){
		u[i] = 0;
		v[i] = 0;
	}

	double** u_matrix = new double*[eq_num];
	for (int i = 0; i < eq_num; i++)
		u_matrix[i] = new double[eq_num];

	FILE* velocity_f;
	fopen_s(&velocity_f, "velocity.txt", "w");

	FILE* velocity_x_f; 
	fopen_s(&velocity_x_f, "velocity_x.txt", "w");

	FILE* velocity_y_f;
	fopen_s(&velocity_y_f, "velocity_y.txt", "w");

	FILE* delta_f;
	fopen_s(&delta_f, "delta.txt", "w");

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

		// ��������� ������� � ������ �����
		fill_matrix_rhs(u_matrix, rhs, nx, ny, u, v, dx, dy, dt, viscosity);

		// ������ ������ ������������ �������
		for (int i = 0; i < 2*nxy; i++)
			uv[i] = rhs[i];

		// ������� �������� u �� ����� ���� ����� ����
		gauss_zeidel(2*nxy, u_matrix, rhs, uv, temp_uv);

		// ������� �������� v �� ����� ���� ����� ��������� u_x + v_y = 0
		// ������� �� ����� ������� v = 0
		for (int j = 0; j < ny; j++)
			v[j] = 0;
		// �� ������� � ������ �������� v = 0
		for (int i = 0; i < nx; i++){
			v[i*ny] = 0;
			v[i*ny + ny - 1] = 0;
		}

		for (int i = 1; i < nx; i++)
		for (int j = 1; j < ny - 1; j++){
			// �� ������ ������� v_y = 0
			if (i == nx - 1)
				v[i * ny + j] = v[i * ny + j - ny];
			// u_x + v_y = 0
			else
				v[i*ny + j] = v[i*ny + j - 1] + dy / dx * (u[i * ny + j] - u[(i - 1) * ny + j]);
		}

		// �������� �� �������
		for (int i = 0; i < nx - 1; i++){
			double min_delta = 1;
			for (int j = 0; j < ny; j++){
				double delta = abs(U(i*dx) - u[i*ny + j]);
				delta /= abs(U(i*dx));
				if (delta < min_delta)
					min_delta = delta;
				if (delta < SPEED_EPS){
					fprintf(delta_f, "%lf\n", j*dy);
					break;
				}
			}
			/*if (min_delta >= SPEED_EPS){
				fprintf(delta_f, "%lf\n", dy*(ny-1));
			}*/
			//printf("Delta %lf\n", min_delta);
		}
		if (step < step_num - 1){
			fprintf(delta_f, "\n\n");
		}


		// ������ ������ ��� ������� ���������
		//if (step % 3 == 0){
			for (int i = 0; i < nx; i+=2){
				for (int j = 0; j < ny; j+=2){
					fprintf(velocity_f, "%lf %lf %lf %lf\n", i*dx, j*dy, u[i*ny + j], v[i*ny + j] * height / width);
					fprintf(velocity_x_f, "%lf %lf %lf %lf\n", i*dx, j*dy, u[i*ny + j], 0.0);
					fprintf(velocity_y_f, "%lf %lf %lf %lf\n", i*dx, j*dy, 0.0, v[i*ny + j] * height / width);
				}
			}

			if (step < step_num - 1){
				fprintf(velocity_f, "\n\n");
				fprintf(velocity_x_f, "\n\n");
				fprintf(velocity_y_f, "\n\n");
			}
	//	}

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

	fclose(delta_f);
	fclose(velocity_f);
	fclose(velocity_x_f);
	fclose(velocity_y_f);

#ifdef DEBUG_STEP
	fclose(debug_step_f);
	fclose(rhs_f);
#endif //DEBUG_STEP
}