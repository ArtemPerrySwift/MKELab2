#include "splain.h"
#include <string>

int splain::read(real* X, real* F, int n, ifstream& in)
{
	int i;
	for (i = 0; i < n; i++)
		in >> X[i] >> F[i];

	return 0;
}
int splain::read(real* X, int n, ifstream& in)
{
	int i;
	for (i = 0; i < n; i++)
		in >> X[i];

	return 0;
}

void splain::fillX_m()
{
	m = n / 4;
	X_m = new real[m];
	X_m[0] = X_n[0] - abs(X_n[1]) * 0.001;
	for (int i = 0; i < m - 1; i++)
		X_m[i] = X_n[4 * i] + X_n[4 * i + 1];

	X_m[m - 1] = X_n[n - 1] + abs(X_n[n - 1]) * 0.001;
}

void splain::init(string file)
{
	//in.open(file1);
	ifstream in;
	in.open(file);
	in >> n;
	X_n = new real[n];
	F = new real[n];
	read(X_n, F, n, in);
	in.close();

	fillX_m();
	di = new real[2 * m];
	ia = new int[2 * m + 1];
	gg = new real[5 * m - 4];
	read(X_m, m, in);
	in.close();
	fill_ia(ia, m);

	int gg_size = ia[2 * m];
	int i;
	for (i = 0; i < gg_size; i++)
		gg[i] = 0;

	real* f = new real[2 * m];
	build_A(di, ia, gg, X_m, X_n, m, n);

	cout << endl;
	bouild_f(f, X_m, X_n, F, n, m);

	calc_LLT(gg, di, ia, 2 * m);
	y = f;
	calc_y(y, f, 2 * m, ia, di, gg);
	q = y;
	calc_x(q, y, 2 * m, ia, di, gg);
}

real splain::local_fun1(real x, real x_i, real h)
{
	return 1 - 3 * (x - x_i) * (x - x_i) / (h * h) + 2 * (x - x_i) * (x - x_i) * (x - x_i) / (h * h * h);
}

real splain::local_fun2(real x, real x_i, real h)
{
	return (x - x_i) - 2 * (x - x_i) * (x - x_i) / h + (x - x_i) * (x - x_i) * (x - x_i) / (h * h);
}

real splain::local_fun3(real x, real x_i, real h)
{
	return 3 * (x - x_i) * (x - x_i) / (h * h) - 2 * (x - x_i) * (x - x_i) * (x - x_i) / (h * h * h);
}

real splain::local_fun4(real x, real x_i, real h)
{
	return -(x - x_i) * (x - x_i) / h + (x - x_i) * (x - x_i) * (x - x_i) / (h * h);
}

void splain::fill_ia(int* ia, int m)
{
	ia[0] = 0;
	ia[1] = 0;
	ia[2] = 1;
	int i;
	int m2 = 2 * m;
	for (i = 2; i < m2; i += 2)
	{
		ia[i + 1] = ia[i] + 2;
		ia[i + 2] = ia[i + 1] + 3;
	}

}

void splain::build_local_m(real M[4][4], real x, real x_i, real h)
{
	real funct[4];
	int i, j;
	funct[0] = local_fun1(x, x_i, h);
	funct[1] = local_fun2(x, x_i, h);
	funct[2] = local_fun3(x, x_i, h);
	funct[3] = local_fun4(x, x_i, h);
	for (i = 0; i < 4; i++)
		for (j = i; j < 4; j++)
		{
			M[i][j] = funct[i] * funct[j];
			M[j][i] = M[i][j];
		}

}

void splain::bouild_f(real* f, real* X_m, real* X_n, real* F, int n, int m)
{
	int i, j;
	real x;
	real h;
	int m2 = 2 * m;
	for (i = 0; i < m2; i++)
		f[i] = 0;
	for (i = 0; i < n; i++)
	{
		x = X_n[i];
		for (j = 0; j < m && X_m[j] < x; j++);
		j--;
		h = X_m[j + 1] - X_m[j];
		f[2 * j] += F[i] * local_fun1(x, X_m[j], h);
		f[2 * j + 1] += F[i] * local_fun2(x, X_m[j], h);
		f[2 * j + 2] += F[i] * local_fun3(x, X_m[j], h);
		f[2 * j + 3] += F[i] * local_fun4(x, X_m[j], h);
	}
}

void splain::build_A(real* di, int* ia, real* gg, real* X_m, real* X_n, int m, int n)
{
	real M_local[4][4];
	int i, j;
	real x;
	real h;
	int gg_size = 5 * m - 4;
	int m2 = 2 * m;
	for (i = 0; i < gg_size; i++)
		gg[i] = 0;
	for (i = 0; i < m2; i++)
		di[i] = 0;
	for (i = 0; i < n; i++)
	{
		x = X_n[i];
		for (j = 0; j < m && X_m[j] < x; j++);
		if (j)
			j--;
		h = X_m[j + 1] - X_m[j];
		build_local_m(M_local, x, X_m[j], h);
		gg[ia[2 * j + 2] - 1] += M_local[0][1];
		gg[ia[2 * j + 3] - 2] += M_local[0][2];
		gg[ia[2 * j + 3] - 1] += M_local[1][2];
		gg[ia[2 * j + 4] - 3] += M_local[0][3];
		gg[ia[2 * j + 4] - 2] += M_local[1][3];
		gg[ia[2 * j + 4] - 1] += M_local[2][3];
		di[2 * j] += M_local[0][0];
		di[2 * j + 1] += M_local[1][1];
		di[2 * j + 2] += M_local[2][2];
		di[2 * j + 3] += M_local[3][3];
	}
	int m_1 = m - 1;
	for (i = 0; i < m_1; i++)
	{
		h = X_m[i + 1] - X_m[i];
		build_local_alpha(M_local, h);
		gg[ia[2 * i + 2] - 1] += ALPHA / 30 * M_local[0][1];
		gg[ia[2 * i + 3] - 2] += ALPHA / 30 * M_local[0][2];
		gg[ia[2 * i + 3] - 1] += ALPHA / 30 * M_local[1][2];
		gg[ia[2 * i + 4] - 3] += ALPHA / 30 * M_local[0][3];
		gg[ia[2 * i + 4] - 2] += ALPHA / 30 * M_local[1][3];
		gg[ia[2 * i + 4] - 1] += ALPHA / 30 * M_local[2][3];
		di[2 * i] += ALPHA / 30 * M_local[0][0];
		di[2 * i + 1] += ALPHA / 30 * M_local[1][1];
		di[2 * i + 2] += ALPHA / 30 * M_local[2][2];
		di[2 * i + 3] += ALPHA / 30 * M_local[3][3];

		build_local_betta(M_local, h);
		gg[ia[2 * i + 2] - 1] += BETTA * M_local[0][1];
		gg[ia[2 * i + 3] - 2] += BETTA * M_local[0][2];
		gg[ia[2 * i + 3] - 1] += BETTA * M_local[1][2];
		gg[ia[2 * i + 4] - 3] += BETTA * M_local[0][3];
		gg[ia[2 * i + 4] - 2] += BETTA * M_local[1][3];
		gg[ia[2 * i + 4] - 1] += BETTA * M_local[2][3];
		di[2 * i] += BETTA * M_local[0][0];
		di[2 * i + 1] += BETTA * M_local[1][1];
		di[2 * i + 2] += BETTA * M_local[2][2];
		di[2 * i + 3] += BETTA * M_local[3][3];
	}
}


void splain::build_local_alpha(real M[4][4], real h)
{
	M[0][0] = 36 / h;
	M[1][0] = M[0][1] = 3;
	M[2][0] = M[0][2] = -36 / h;
	M[3][0] = M[0][3] = 3;
	M[1][1] = 4 * h;
	M[2][1] = M[1][2] = -3;
	M[3][1] = M[1][3] = -h;
	M[2][2] = 36 / h;
	M[3][2] = M[2][3] = -3;
	M[3][3] = 4 * h;
}

void splain::build_local_betta(real M[4][4], real h)
{
	M[0][0] = 60 / (h * h * h);
	M[1][0] = M[0][1] = 30 / (h * h);
	M[2][0] = M[0][2] = -60 / (h * h * h);
	M[3][0] = M[0][3] = 30 / (h * h);
	M[1][1] = 16 / h;
	M[2][1] = M[1][2] = -30 / (h * h);
	M[3][1] = M[1][3] = 14 / h;
	M[2][2] = 60 / (h * h * h);
	M[3][2] = M[2][3] = -30 / (h * h);
	M[3][3] = 16 / h;
}

bool splain::calc_LLT(real* al, real* di, int* ia, int n)
{
	int i_beg, i_end;
	int j_beg, j_end;
	int i, j, ij, ij_0;
	int ik, kj;
	real sum, sum_di;
	int ind_i, ind_j, r;

	for (i = 0; i < n; i++)
	{
		i_beg = ia[i];
		i_end = ia[i + 1];

		ij = i_beg;
		ij_0 = i - (i_end - i_beg);

		sum_di = 0;


		for (j = ij_0; j < i; j++, ij++)
		{
			i_beg = ia[i];
			j_beg = ia[j];
			j_end = ia[j + 1];

			ind_i = i - (i_end - i_beg);
			ind_j = j - (j_end - j_beg);
			r = ind_i - ind_j;

			r > 0 ? j_beg += r : i_beg -= r;

			sum = 0;
			for (ik = i_beg, kj = j_beg; kj < j_end; ik++, kj++)
				sum += al[ik] * al[kj];

			al[ij] = (al[ij] - sum) / di[j];
			sum_di += al[ij] * al[ij];
		}

		if (di[i] <= sum_di)
		{
			cout << "Error : Invalid data. Matrix is not positively de-fined";
			return true;
		}


		di[i] = sqrt(di[i] - sum_di);
	}

	return false;
}

void splain::calc_y(real* y, real* f, int n, int* ia, real* di, real* al)
{
	real sum;
	int i_beg, i_end, j, ij;
	for (int i = 0; i < n; i++) // Ğåøàåì ÑËÀÓ L*y = F
	{
		sum = 0;
		i_beg = ia[i];
		i_end = ia[i + 1];
		j = i - (i_end - i_beg);
		for (ij = i_beg; ij < i_end; ij++, j++)
			sum += al[ij] * y[j];

		y[i] = (f[i] - sum) / di[i];

	}
}
void splain::calc_x(real* x, real* y, int n, int* ia, real* di, real* al)
{
	int i_beg, i_end, j, ij;
	for (int i = n - 1; i >= 0; i--) // ğåøàåì ÑËÀÓ L^T * x = y
	{
		x[i] = y[i] / di[i];
		i_beg = ia[i];
		i_end = ia[i + 1];
		j = i - (i_end - i_beg);
		for (ij = i_beg; ij < i_end; ij++, j++)
			y[j] -= al[ij] * x[i];

	}
}

real splain::count_splain(real x)
{
	int i;
	int m2 = 2 * m;
	for (i = 0; i < m2 && X_m[i] < x; i++);
	i--;
	real h = X_m[i + 1] - X_m[i];
	return q[2 * i] * local_fun1(x, X_m[i], h) + q[2 * i + 1] * local_fun2(x, X_m[i], h) + q[2 * i + 2] * local_fun3(x, X_m[i], h) + q[2 * i + 3] * local_fun4(x, X_m[i], h);
}