#include "FinalElemetsBasic.h"
#include "DifferentEquParams.h"

FinalElemBase::FinalElemBase()
{
	grads[0] = &basicGradFun2D_0_;
	grads[1] = &basicGradFun2D_1_;
	grads[2] = &basicGradFun2D_2_;
	grads[3] = &basicGradFun2D_3_;

	basicFun2D[0] = &basicFun2D_0_;
	basicFun2D[1] = &basicFun2D_1_;
	basicFun2D[2] = &basicFun2D_2_;
	basicFun2D[3] = &basicFun2D_3_;
}

double FinalElemBase::countGausse2ForG(int i, int j, int l)
{
	int k, m;
	double sum = 0;
	for (k = 0; k < GAUSSE_ROOT_NUM; k++)
		for (m = 0; m < GAUSSE_ROOT_NUM; m++)
			sum += GAUSSE_W * GAUSSE_W * (*grads[i])(GAUSSE_ROOT_X[k], GAUSSE_ROOT_Y[m], hx, xKnot, hy, yKnot) * 
			(*grads[j])(GAUSSE_ROOT_X[k], GAUSSE_ROOT_Y[m], hx, xKnot, hy, yKnot) * (*basicFun2D[l])(GAUSSE_ROOT_X[k], GAUSSE_ROOT_Y[m], hx, xKnot, hy, yKnot);

	return hx * hy / 4;
}

double FinalElemBase::basicFun1D_0(double h, double sqrPoint, double arg) { return (sqrPoint - arg) / h; }
double FinalElemBase::basicFun1D_1(double h, double sqrPoint, double arg) { return (arg - sqrPoint) / h; }

double FinalElemBase::basicFun2D_0(double x, double y){ return basicFun1D_0(hx, xKnot[1], x) * basicFun1D_0(hy, yKnot[1], y);}
double FinalElemBase::basicFun2D_1(double x, double y){ return basicFun1D_1(hx, xKnot[0], x) * basicFun1D_0(hy, yKnot[1], y);}
double FinalElemBase::basicFun2D_2(double x, double y){ return basicFun1D_0(hx, xKnot[1], x) * basicFun1D_1(hy, yKnot[0], y);}
double FinalElemBase::basicFun2D_3(double x, double y){ return basicFun1D_1(hx, xKnot[0], x) * basicFun1D_1(hy, yKnot[0], y);}

double FinalElemBase::basicFun2D_0_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]) { return basicFun1D_0(hx, xKnot[1], x) * basicFun1D_0(hy, yKnot[1], y); }
double FinalElemBase::basicFun2D_1_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]) { return basicFun1D_1(hx, xKnot[0], x) * basicFun1D_0(hy, yKnot[1], y); }
double FinalElemBase::basicFun2D_2_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]) { return basicFun1D_0(hx, xKnot[1], x) * basicFun1D_1(hy, yKnot[0], y); }
double FinalElemBase::basicFun2D_3_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]) { return basicFun1D_1(hx, xKnot[0], x) * basicFun1D_1(hy, yKnot[0], y); }

double FinalElemBase::basicGradFun2D_0(double x, double y){return -basicFun1D_0(hy, yKnot[1], y) / hx - basicFun1D_0(hx, xKnot[1], x) / hy;}
double FinalElemBase::basicGradFun2D_1(double x, double y){return basicFun1D_0(hy, yKnot[1], y) / hx - basicFun1D_1(hx, xKnot[0], x) / hy;}
double FinalElemBase::basicGradFun2D_2(double x, double y){return -basicFun1D_1(hy, yKnot[0], y) / hx + basicFun1D_0(hx, xKnot[1], x) / hy;}
double FinalElemBase::basicGradFun2D_3(double x, double y){return basicFun1D_1(hy, yKnot[0], y) / hx + basicFun1D_1(hx, xKnot[0], x) / hy;}

double FinalElemBase::basicGradFun2D_0_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]) { return -basicFun1D_0(hy, yKnot[1], y) / hx - basicFun1D_0(hx, xKnot[1], x) / hy; }
double FinalElemBase::basicGradFun2D_1_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]) { return basicFun1D_0(hy, yKnot[1], y) / hx - basicFun1D_1(hx, xKnot[0], x) / hy; }
double FinalElemBase::basicGradFun2D_2_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]) { return -basicFun1D_1(hy, yKnot[0], y) / hx + basicFun1D_0(hx, xKnot[1], x) / hy; }
double FinalElemBase::basicGradFun2D_3_(double x, double y, double hx, double xKnot[N_KNOTS_1D], double hy, double yKnot[N_KNOTS_1D]) { return basicFun1D_1(hy, yKnot[0], y) / hx + basicFun1D_1(hx, xKnot[0], x) / hy; }


double FinalElemBase::count_dG_dq(int i, int r, int j)
{
	return dlambda_du(q[j], xKnot[j % 2], yKnot[j / 2]) * countGausse2ForG(i, r, j);
}

double FinalElemBase::f(double x, double y)
{
	return material.toku;
}

double FinalElemBase::lambda(double x, double y)
{
	return 1 / count_mu(x, y);
}

double FinalElemBase::du_dx(double x, double y)
{
	return (basicFun1D_0(hy, knots[0].y, y) * (q[1] - q[0]) + basicFun1D_1(hy, knots[3].y, y) * (q[1] - q[0])) / hx;
}

double FinalElemBase::du_dy(double x, double y)
{
	return (basicFun1D_0(hx, knots[0].x, x) * (q[1] - q[0]) + basicFun1D_1(hx, knots[1].x, x) * (q[1] - q[0])) / hy;
}

double FinalElemBase::magn_ind(double x, double y)
{
	return sqrt(pow(du_dx(x, y), 2) + pow(du_dy(x, y), 2));
}

double FinalElemBase::count_mu(double x, double y)
{
	return material.isMuConst ? material.mu : material.muSplain.count_splain(magn_ind(x, y));
}

double FinalElemBase::dlambda_du(double u, double x, double y)
{
	return 0;
}

void FinalElemBase::setFinalElement(FinalElement finalElement, KnotsStorage knotsStorage, Material material, double* q)
{
	this->material = material;
	this->finalElement = finalElement;
	knots[0] = knotsStorage.knots[finalElement.ver1];
	knots[1] = knotsStorage.knots[finalElement.ver2];
	knots[2] = knotsStorage.knots[finalElement.ver3];
	knots[3] = knotsStorage.knots[finalElement.ver4];

	/*Пока без понятия как там рополагаются x и y так что это место довольно опасное, если работать не будет, то  возможно ошибка здесь*/
	globalInd[0] = finalElement.ver1;
	globalInd[1] = finalElement.ver2;
	globalInd[2] = finalElement.ver3;
	globalInd[3] = finalElement.ver4;

	/*Упорядочиваем узлы для локальной нумерации*/
	int i, j;
	for (i = 1; i < N_KNOTS_2D; i++)
	{
		if (knots[i].x < knots[0].x && knots[i].y < knots[0].y)
		{
			swap(knots[i], knots[0]);
			swap(globalInd[i], globalInd[0]);
		}

		if (knots[i].x > knots[LAST_LOCAL_IND].x && knots[i].y > knots[LAST_LOCAL_IND].y)
		{
			swap(knots[i], knots[LAST_LOCAL_IND]);
			swap(globalInd[i], globalInd[LAST_LOCAL_IND]);
		}
	}

	if (knots[2].y < knots[1].y)
	{
		swap(knots[1], knots[2]);
		swap(globalInd[1], globalInd[2]);
	}

	xKnot[0] = knots[0].x;
	xKnot[1] = knots[1].x;

	yKnot[0] = knots[0].y;
	yKnot[1] = knots[2].y;

	for (i = 0; i < N_KNOTS_2D; i++){ this->q[i] = q[globalInd[i]]; }

	hx = xKnot[1] - xKnot[0];
	hy = yKnot[1] - yKnot[0];

	avrX = (knots[1].x + knots[0].x) / 2;
	avrY = (knots[2].y + knots[0].y) / 2;

	GAUSSE_ROOT_X[0] = avrX + GAUSSE_ROOTS[0] * hx / 2;
	GAUSSE_ROOT_X[1] = avrX + GAUSSE_ROOTS[1] * hx / 2;

	GAUSSE_ROOT_Y[0] = avrY + GAUSSE_ROOTS[0] * hy / 2;
	GAUSSE_ROOT_Y[1] = avrY + GAUSSE_ROOTS[1] * hy / 2;
}

void FinalElemBase::buildLocalM(double M[N_KNOTS_2D][N_KNOTS_2D], double gamma)
{
	int i, j;
	for (i = 0; i < N_KNOTS_2D; i++)
	{
		for (j = i; j < N_KNOTS_2D; j++)
		{
			M[i][j] = gamma*LOCAL_MATRIX_M_1[i % 2][j % 2] * LOCAL_MATRIX_M_1[i / 2][j / 2] * hx * hy / 36.0;
			M[j][i] = M[i][j];
		}
	}
}

void FinalElemBase::buildLocalG(double G[N_KNOTS_2D][N_KNOTS_2D])
{
	int i, j, l;
	double sum;
	for (i = 0; i < N_KNOTS_2D; i++)
	{
		for (j = i; j < N_KNOTS_2D; j++)
		{
			sum = 0;
			for (l = 0; l < N_KNOTS_2D; l++)
				sum += lambda(knots[l].x, knots[l].y) * countGausse2ForG(i, j, l);

			G[i][j] = sum;
			G[j][i] = G[i][j];
		}
	}
}

void FinalElemBase::multiply_M_and_F(double M[N_KNOTS_2D][N_KNOTS_2D], double F[N_KNOTS_2D], double B[N_KNOTS_2D])
{
	int i, j;
	double sum;
	for (i = 0; i < N_KNOTS_2D; i++)
	{
		B[i] = 0;
		for (j = 0; j < N_KNOTS_2D; j++)
			B[i] += M[i][j] * F[j];
	}
}

void FinalElemBase::buildLocalB(double B[N_KNOTS_2D])
{
	double M[N_KNOTS_2D][N_KNOTS_2D];
	double F[N_KNOTS_2D];
	buildLocalM(M, 1.0);
	fillF(F);
	multiply_M_and_F(M, F, B);
}

void FinalElemBase::fillF(double F[N_KNOTS_2D])
{
	for (int i = 0; i < N_KNOTS_2D; i++)
	    F[i] = DifferentEquParams::f(knots[i].x, knots[i].y);
}

void FinalElemBase::addNewt(double G[N_KNOTS_2D][N_KNOTS_2D], double B[N_KNOTS_2D])
{
	double sum;
	double sumB;
	int i, j, r;
	for (i = 0; i < N_KNOTS_2D; i++)
	{
		sumB = 0;
		for (j = 0; j < N_KNOTS_2D; j++)
		{
			sum = 0;
			for (r = 0; r < N_KNOTS_2D; r++)
			{
				sum += count_dG_dq(i, r, j)*q[r];
			}
			G[i][j] += sum;
			sumB += sum * q[j];
		}
		B[i] += sumB;
	}
}

void FinalElemBase::buildNewt(double G[N_KNOTS_2D][N_KNOTS_2D], double B[N_KNOTS_2D])
{
	buildLocalG(G);
	buildLocalB(B);
	addNewt(G, B);
}