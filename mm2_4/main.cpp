#define _USE_MATH_DEFINES
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

typedef vector<double> vector_double;

// #define DEBUG_MATRIX
// #define DEBUG_STEP



#define EPS 0.0001

// ����� ������ ������� ���� � ���������
bool converge(vector_double &xk, vector_double &xkp, int n)
{
	double norm = 0;
	for (int i = 0; i < n; i++)
		norm += (xk[i] - xkp[i])*(xk[i] - xkp[i]);
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
void gauss_zeidel(int n, vector<vector<double>> &a, vector_double &b, vector_double &x, vector_double &p)
{
	do
	{
		copy(x.begin(), x.end(), p.begin());
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

inline double U(double x)
{
	return cosh(x);
}

inline double U_x(double x) 
{
	return sinh(x);
}

void fill_matrix_rhs(vector<vector<double>> &mt, vector_double &rhs, int nx, int ny, vector_double &u, vector_double &v, double dx, double dy, double dt, double viscosity)
{
	int nxy = nx * ny;
	for (int i = 0; i < nxy; i++)
		fill(mt[i].begin(), mt[i].end(), 0);		

	for (int i = 0; i < nx; i++)
	{
		// �� ������ ������� u = 0
		mt[i*ny][i*ny] = 1;
		rhs[i * ny] = 0;

		// �� ������� ������� u = U(x)
		mt[i*ny + ny - 1][i*ny + ny - 1] = 1;
		rhs[i*ny + ny - 1] = U(i*dx);
	}

	for (int j = 1; j < ny - 1; j++)
	{
		// �� ����� ������� u = U(0)
		mt[j][j] = 1;
		rhs[j] = U(0);

		// �� ������ ������� u_x = 0
		mt[(nx - 1)*ny + j][(nx - 1)*ny + j] = 1;
		mt[(nx - 1)*ny + j][(nx - 2)*ny + j] = -1;
		rhs[(nx - 1)*ny + j] = 0;
	}

	double coeff = viscosity / dy / dy;

	// ����� �����, � �������� ���� ����� �������� ����������� u_x � u_y
	for (int i = 1; i < nx - 1; i++)
		for (int j = 1; j < ny - 1; j++)
		{
			int eq_ind = i * ny + j;
			mt[eq_ind][eq_ind] = 1 / dt + 2 * coeff + u[eq_ind] / dx + v[eq_ind] / dy;
			mt[eq_ind][eq_ind + 1] = -coeff;
			mt[eq_ind][eq_ind - 1] = -coeff - v[eq_ind] / dy;
			mt[eq_ind][eq_ind - ny] = -u[eq_ind] / dx;
			rhs[eq_ind] = u[eq_ind] / dt + U(i*dx) * U_x(i*dx);
		}
}

int main()
{
	// �������� �������� - ������ ��������, ��� ������� �������� ����� ���������� �������� � ����� ����� ���� ����������
	// �� ��� �� �����
	double viscosity = 0.01;

	// ������ � ������ �������
	double width = M_PI / 2.0;
	double height = 1.2;

	// ����� �����, ��� �������� ���������� ������ � ����� ��������� �����
	double total_time = 0.3;
	int step_num = 50;
	double dt = total_time / step_num;

	// ����� ��� ������������� x � ������ ������� total_time
	// ����� ������ ������ ��������, ����������� ���� slice.plt, ����� �� ���� ������ ������
	const int slices_num = 2;
	double slices_x[slices_num] = { 0.1, 0.2 };

	// ������� �����
	int nx = 41;
	int ny = 41;
	int nxy = nx*ny;

	double dx = width / (nx - 1.0);
	double dy = height / (ny - 1.0);

	// ������������ (u) � �������������� (v) ������������ ��������
	// ��� ������� (� ������ �����) ����������� �� ��� y � ������ �������
	// �.�. �������� u[0], u[5] ������������� ��������� (0,0) � (0,5) 
	vector_double u(nxy), v(nxy), temp_u(nxy), rhs(nxy);	

	// ��������� ������������� �������
	for (int i = 0; i < nxy; i++)
	{
		u[i] = 0;
		v[i] = 0;
	}

	// �� ����� ���������, ��� ������� ������ ������� �� ���������� �������������
	// ����� ���������� � ������� ���������� u0(x,y)=U(x) (���������������� ���� ����)
	
	/*for(int i=0; i<nx; i++)
		for (int j = 0; j < ny; j++)
		{
			u[i*ny + j] = U(i*dx);
		}
*/
	// ������� ������������� ����
	vector<vector<double>> u_matrix(nxy,vector_double(nxy));	

	// ���� � ������� ��� �������� ��������� ��������
	FILE* velocity_f;
	fopen_s(&velocity_f, "velocity.txt", "w");
	// ���� � ������� ��� �������� �������� �� x
	FILE* velocity_x_f;
	fopen_s(&velocity_x_f, "velocity_x.txt", "w");
	// ���� � ������� ��� �������� �������� �� y
	FILE* velocity_y_f;
	fopen_s(&velocity_y_f, "velocity_y.txt", "w");
	// ���� ��� ���������� ����� �� y ��� ������������� x � �������� ������ �������
	FILE* slice_f;

#ifdef DEBUG_MATRIX
	FILE* matrix_f;
	matrix_f = fopen("matrix.txt", "w");
	for (int i = 0; i < nxy; i++) 
	{
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

	for (int step = 0; step < step_num; step++)
	{
		printf("Step #%d of %d\n", step + 1, step_num);
		// ��������� ������� � ������ �����
		fill_matrix_rhs(u_matrix, rhs, nx, ny, u, v, dx, dy, dt, viscosity);

		// ������ ������ ������������ �������
		for (int i = 0; i < nxy; i++)
			u[i] = rhs[i];

		// ������� �������� u �� ����� ���� ����� ����
		gauss_zeidel(nxy, u_matrix, rhs, u, temp_u);

		// ������� �������� v �� ����� ���� ����� ��������� u_x + v_y = 0
		// ������� �� ����� ������� v = 0
		for (int j = 0; j < ny; j++)
			v[j] = 0;
		// �� ������� � ������ �������� v = 0
		for (int i = 0; i < nx; i++)
		{
			v[i*ny] = 0;
			v[i*ny + ny - 1] = 0;
		}

		for (int i = 1; i < nx; i++)
			for (int j = 1; j < ny - 1; j++)
			{
				// �� ������ ������� v_x = 0
				if (i == nx - 1)
					v[i * ny + j] = v[i * ny + j - ny];
				// u_x + v_y = 0
				else
					v[i * ny + j] = v[i * ny + j - 1] - dy / dx * (u[i * ny + j] - u[(i - 1) * ny + j]);
			}

		// �������� ���� � ������������� x ��� �������� ������� �������
		// ����� �������, ������ ����� ���������� � ���������� ���������� �������
		// � ���� ����, ��� ���� ����������������, � ����� � ��� ����������� �������� ������� (�.�. � ������ ������� total_time)
		fopen_s(&slice_f, "slice.txt", "w");
		for (int i = 0; i < ny; i++)
		{
			fprintf(slice_f, "%lf ", i*dy);
			for (int slice = 0; slice < slices_num; slice++)
			{
				fprintf(slice_f, "%lf ", u[int(nx * (slices_x[slice] / width))* ny + i]);
			}
			fprintf(slice_f, "\n");
		}
		fclose(slice_f);

		// ���� ����� �������������, �� ��� ������������� � �������� ������� ���� (��� ���� �������)
		// � ������� ����� ������ ����� ��������, ����� �������� ������ �������� ��� ����� � ���������� �������
		printf("%lf %lf\n", u[nx / 2 * ny + ny - 1], u[nx / 2 * ny + ny - 2]);

		// �������� ����� � ������ ���� ������� ����� ��� �������� ���� ���������
		for (int i = 0; i < nx; i += 4)
		{
			for (int j = 0; j < ny; j += 4)
			{
				fprintf(velocity_f, "%lf %lf %lf %lf\n", i*dx, j*dy, 0.2*u[i*ny + j], 0.2*v[i*ny + j] * height / width);
				fprintf(velocity_x_f, "%lf %lf %lf %lf\n", i*dx, j*dy, 0.2*u[i*ny + j], 0.0);
				fprintf(velocity_y_f, "%lf %lf %lf %lf\n", i*dx, j*dy, 0.0, 0.2*v[i*ny + j] * height / width);
			}
		}

		if (step < step_num - 1)
		{
			fprintf(velocity_f, "\n\n");
			fprintf(velocity_x_f, "\n\n");
			fprintf(velocity_y_f, "\n\n");
		}

#ifdef DEBUG_STEP
		fprintf(rhs_f, "RHS step %i\n", step);
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++)
				fprintf(rhs_f, "%lf ", rhs[i*ny + j]);
			fprintf(rhs_f, "\n");
		}
		fprintf(rhs_f, "\n");

		fprintf(debug_step_f, "Step %i\nU\n", step);
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++)
				fprintf(debug_step_f, "%lf ", u[i*ny + j]);
			fprintf(debug_step_f, "\n");
		}
		fprintf(debug_step_f, "\nV\n");
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++)
				fprintf(debug_step_f, "%lf ", v[i*ny + j]);
			fprintf(debug_step_f, "\n");
		}
		fprintf(debug_step_f, "\n");
#endif // DEBUG_STEP
	}

	fclose(velocity_f);
	fclose(velocity_x_f);
	fclose(velocity_y_f);

#ifdef DEBUG_STEP
	fclose(debug_step_f);
	fclose(rhs_f);
#endif //DEBUG_STEP
}