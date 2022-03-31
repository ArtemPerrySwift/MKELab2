#pragma once
#include "mesh.h"
using namespace meshspace;

class FinalElemBase
{
public:
	static const int N_KNOTS_1D = 2;
	static const int N_KNOTS_2D = 4;
	const double LOCAL_MATRIX_G_1[N_KNOTS_1D][N_KNOTS_1D] = { 1, -1, -1, 1 };
	const double LOCAL_MATRIX_M_1[N_KNOTS_1D][N_KNOTS_1D] = { 2, 1, 2, 1 };

	void setFinalElement(FinalElement finalElement, KnotsStorage knotsStorage, Material material, double* q);
	void buildLocalM(double M[N_KNOTS_2D][N_KNOTS_2D], double gamma);
	void buildLocalG(double G[N_KNOTS_2D][N_KNOTS_2D]);
	void addNewt(double G[N_KNOTS_2D][N_KNOTS_2D], double B[N_KNOTS_2D]);
	void buildNewt(double G[N_KNOTS_2D][N_KNOTS_2D], double B[N_KNOTS_2D]);
	void buildLocalB(double B[N_KNOTS_2D]);
	void multiply_M_and_F(double M[N_KNOTS_2D][N_KNOTS_2D], double F[N_KNOTS_2D], double B[N_KNOTS_2D]);
	void fillF(double F[N_KNOTS_2D]);

	double basicFun2D_0(double x, double y);
	double basicFun2D_1(double x, double y);
	double basicFun2D_2(double x, double y);
	double basicFun2D_3(double x, double y);

	double basicGradFun2D_0(double x, double y);
	double basicGradFun2D_1(double x, double y);
	double basicGradFun2D_2(double x, double y);
	double basicGradFun2D_3(double x, double y);

	static double basicFun2D_0_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]);
	static double basicFun2D_1_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]);
	static double basicFun2D_2_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]);
	static double basicFun2D_3_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]);

	static double basicGradFun2D_0_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]);
	static double basicGradFun2D_1_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]);
	static double basicGradFun2D_2_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]);
	static double basicGradFun2D_3_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]);

	FinalElemBase();
	int globalInd[N_KNOTS_2D];

	double f(double x, double y);

	double lambda(double x, double y);

	double du_dx(double x, double y);

	double du_dy(double x, double y);

	double magn_ind(double x, double y);

	double count_mu(double x, double y);

	double dlambda_du(double u, double x, double y);

	double u(double x, double y);

private:
	FinalElement finalElement;
	Material material;
	double xKnot[N_KNOTS_1D];
	double yKnot[N_KNOTS_1D];
	Knot knots[N_KNOTS_2D];
	double q[N_KNOTS_2D];
	double GAUSSE_ROOT_X[2];
	double GAUSSE_ROOT_Y[2];
	bool isLinear;

	double(* grads[N_KNOTS_2D])(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]);
	double(* basicFun2D[N_KNOTS_2D])(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]);

	double hx, hy;
	double avrY, avrX;
	static double basicFun1D_0(double h, double sqrPoint, double arg);
	static double basicFun1D_1(double h, double sqrPoint, double arg);

	static const int GAUSSE_ROOT_NUM = 2;
	static const int LAST_LOCAL_IND = N_KNOTS_2D - 1;
	const double GAUSSE_ROOTS[GAUSSE_ROOT_NUM] = { 0.5773502691896257, -0.5773502691896257 };
	const double GAUSSE_W = 1;

	double count_dG_dq(int i, int r, int j);
	double countGausse2ForG(int i, int j, int l);

};

