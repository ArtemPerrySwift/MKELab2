#include <cmath>
#include <iostream>
#include <iomanip>
#include "slae.h"
#include "array.h"

using namespace arrayspace;
namespace slae
{
	void SLAE::init(SparseMatrixAsym A)
	{
		b = new double[A.n];
		init(A, b);

	}

	void SLAE::init(SparseMatrixAsym A, double* b)
	{
		r = new double[A.n];
		z = new double[A.n];
		p = new double[A.n];
		f = new double[A.n];
		buf_v = new double[A.n];
		buf_v1 = new double[A.n];
		this->A = A;
		this->b = b;
		LU.copy(A);
	}

	int SLAE::count_LOS(real* x, int maxiter, real eps)
	{
		LU.decomp_mat_LU(A);
		bool is_L_diag_1 = false;
		bool is_U_diag_1 = true;

		int n = A.n;
		int i, j, k;
		real norm_ff = sqrt(scal(f, f, n));
		real err;

		bool fl = true;
		i = 0;
		while (fl)
		{
			A.mult(x, buf_v);
			for (j = 0; j < n; j++)
				buf_v[j] = f[j] - buf_v[j];

			calc_Lx(LU, r, buf_v, is_L_diag_1);
			calc_Ux(LU, z, r, is_U_diag_1);

			A.mult(z, buf_v);
			calc_Lx(LU, p, buf_v, is_L_diag_1);

			real a_k, b_k;
			real skal_pp;

			err = sqrt(scal(r, r, n)) / norm_ff;
			fl = i < maxiter&& err > eps;
			for (k = 0; fl && k < n; i++, k++)
			{

				skal_pp = scal(p, p, n);
				a_k = scal(p, r, n) / skal_pp;

				for (j = 0; j < n; j++)
				{
					x[j] += a_k * z[j];
					r[j] -= a_k * p[j];
				}

				calc_Ux(LU, buf_v, r, is_U_diag_1);

				A.mult(buf_v, buf_v1);
				calc_Lx(LU, buf_v, buf_v1, is_L_diag_1);

				b_k = -scal(p, buf_v, n) / skal_pp;
				calc_Ux(LU, buf_v1, r, is_U_diag_1);
				for (j = 0; j < n; j++)
				{
					z[j] = buf_v1[j] + b_k * z[j];
					p[j] = buf_v[j] + b_k * p[j];
				}
				err = sqrt(scal(r, r, n)) / norm_ff;
				fl = i < maxiter&& err > eps;
				cout << setprecision(14) << "iteration: " << i << "; err: " << err << endl;
				//itarations = i;
				//nev = err;
			}
		}
		return 0;
	}

	int SLAE::calc_Ux(SparseMatrixAsym& LU, real* x, real* f, bool is_diag_1)
	{
		int n = LU.n;

		real* U_di = LU.di;
		real* U_ggu = LU.ggu;
		int* U_ig = LU.ig;
		int* U_jg = LU.jg;

		int i, j;
		int i_beg, i_end;
		int ij;
		for (i = 0; i < n; i++)
			x[i] = f[i];

		if (is_diag_1)
			for (i = n - 1; i > -1; i--)
			{
				i_beg = U_ig[i];
				i_end = U_ig[i + 1];
				for (ij = i_beg; ij < i_end; ij++)
				{
					j = U_jg[ij];
					x[j] -= U_ggu[ij] * x[i];
				}
			}
		else
			for (i = n - 1; i > -1; i--)
			{
				i_beg = U_ig[i];
				i_end = U_ig[i + 1];
				x[i] /= U_di[i];
				for (ij = i_beg; ij < i_end; ij++)
				{
					j = U_jg[ij];
					x[j] -= U_ggu[ij] * x[i];
				}
			}

		return 0;
	}

	int SLAE::calc_Lx(SparseMatrixAsym& LU, real* x, real* f, bool is_diag_1)
	{
		int n = LU.n;

		real* L_di = LU.di;
		real* L_ggl = LU.ggl;
		int* L_ig = LU.ig;
		int* L_jg = LU.jg;

		int i, j;
		int i_beg, i_end;
		int ij;
		if (is_diag_1)
			for (i = 0; i < n; i++)
			{
				x[i] = f[i];
				i_beg = L_ig[i];
				i_end = L_ig[i + 1];
				for (ij = i_beg; ij < i_end; ij++)
				{
					j = L_jg[ij];
					x[i] -= L_ggl[ij] * x[j];
				}
			}
		else
			for (i = 0; i < n; i++)
			{
				x[i] = f[i];
				i_beg = L_ig[i];
				i_end = L_ig[i + 1];
				for (ij = i_beg; ij < i_end; ij++)
				{
					j = L_jg[ij];
					x[i] -= L_ggl[ij] * x[j];
				}
				x[i] /= L_di[i];
			}
		return 0;
	}

	void SLAE::setOneVariableSolve(int iVar, double varMean)
	{
		int* ig = A.ig;
		int* jg = A.jg;
		double* ggl = A.ggl;
		double* ggu = A.ggu;
		double* di = A.di;
		int n = A.n;

		varMean = 0;
		int ind_beg = ig[iVar];
		int ind_end = ig[iVar + 1];
		int i_ind, b_ind;

		int ind_current;
		for (i_ind = ind_beg; i_ind < ind_end; i_ind++)
		{
			b_ind = jg[i_ind];

			b[b_ind] -= ggu[i_ind] * varMean;
			ggu[i_ind] = 0;
		}

		ind_current = iVar;
		for (iVar++; iVar < n; iVar++)
		{
			ind_beg = ig[iVar];
			ind_end = ig[iVar + 1];
			for (i_ind = ind_beg; i_ind < ind_end; i_ind++)
			{
				if (jg[i_ind] == ind_current)
				{

					b[iVar] -= ggl[i_ind] * varMean;
					ggl[i_ind] = 0;
				}
			}
		}
		di[ind_current] = 1;
		b[ind_current] = varMean;
	}
}
