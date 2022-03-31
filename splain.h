#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;
typedef double real;

const real ALPHA = 0;
const real BETTA = 0;

int read(real* X, real* F, int n, ifstream& in);
int read(real* X, int n, ifstream& in);
real local_fun1(real x, real x_i, real h);
real local_fun2(real x, real x_i, real h);
real local_fun3(real x, real x_i, real h);
real local_fun4(real x, real x_i, real h);
void fill_ia(int* ia, int m);
void build_local_m(real M[4][4], real x, real x_i, real h);
void bouild_f(real* f, real* X_m, real* X_n, real* F, int n, int m);
void build_A(real* di, int* ia, real* gg, real* X_m, real* X_n, int m, int n);
void build_local_alpha(real M[4][4], real h);
void build_local_betta(real M[4][4], real h);
bool calc_LLT(real* al, real* di, int* ia, int n);
void calc_y(real* y, real* f, int n, int* ia, real* di, real* al);
void calc_x(real* x, real* y, int n, int* ia, real* di, real* al);
real count_splain(real x, real* X_m, real* q, int m);

struct splain
{
	int n;
	int m;
	real* X_n;
	real* X_m;
	real* F;
	real* di;
	int* ia;
	real* gg;
	real* q;
	real* y;


	int read(real* X, real* F, int n, ifstream& in);
	int read(real* X, int n, ifstream& in);
	real local_fun1(real x, real x_i, real h);
	real local_fun2(real x, real x_i, real h);
	real local_fun3(real x, real x_i, real h);
	real local_fun4(real x, real x_i, real h);
	void fill_ia(int* ia, int m);
	void build_local_m(real M[4][4], real x, real x_i, real h);
	void bouild_f(real* f, real* X_m, real* X_n, real* F, int n, int m);
	void build_A(real* di, int* ia, real* gg, real* X_m, real* X_n, int m, int n);
	void build_local_alpha(real M[4][4], real h);
	void build_local_betta(real M[4][4], real h);
	bool calc_LLT(real* al, real* di, int* ia, int n);
	void calc_y(real* y, real* f, int n, int* ia, real* di, real* al);
	void calc_x(real* x, real* y, int n, int* ia, real* di, real* al);
	real count_splain(real x);
	void init(string file1);
	void fillX_m();
};
/*
int main()
{
	int n;
	int m;
	real* X_n;
	real* X_m;
	real* F;
	ifstream in;
	real* di;
	int* ia;
	real* gg;
	real* q;
	real* y;
	in.open("input1.txt");
	in >> m;
	X_m = new real[m];
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
	in.open("input2.txt");
	in >> n;
	X_n = new real[n];
	F = new real[n];
	read(X_n, F, n, in);
	in.close();

	real* f = new real[2 * m];
	build_A(di, ia, gg, X_m, X_n, m, n);

	cout << endl;
	bouild_f(f, X_m, X_n, F, n, m);

	calc_LLT(gg, di, ia, 2 * m);
	y = f;
	calc_y(y, f, 2 * m, ia, di, gg);
	q = y;
	calc_x(q, y, 2 * m, ia, di, gg);

	cout << "f" << endl;
	for (real x = 1; x < 5; x += 0.2)
		cout << setprecision(12) << x << "  " << count_splain(x, X_m, q, m) << endl;;
	return 0;
}
*/


