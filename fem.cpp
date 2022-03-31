#include "fem.h"
#include "DifferentEquParams.h"
#include "array.h"
#include "GlobalMatAssembler.h"
#include <iostream>

using namespace arrayspace;

void FEM::init(Mesh mesh)
{
	this->mesh = mesh;
	SparseMatrixAsym A;
	A.buildPortrait(mesh);
	slae.init(A);
	q = new double[A.n];
	qPrev = new double[A.n];
	bForCheck = new double[A.n];
	
	descr = bForCheck;
	arrayspace::fill_vec(q, A.n, 0.0);

}

void FEM::addFirstConditions()
{
	int* Iknots = mesh.knotsWithFirstConditionStorage.IKnots;
	int kt1 = mesh.knotsWithFirstConditionStorage.kt1;

	Knot* knots = mesh.knotsStorage.knots;
	double uMean;
	int ind;
	for (int i = 0; i < kt1; i++)
	{
		ind = Iknots[i];
		uMean = DifferentEquParams::u1(knots[ind].x, knots[ind].y);
		slae.setOneVariableSolve(ind, uMean);
	}
}

void FEM::LinearTask()
{
	slae.A.fillMatrix(0.0);
	simpleAssembler.assembleGlobalMatrix(slae, mesh, q);
	addFirstConditions();
	slae.count_LOS(q, 1000, 1e-12);
}

void FEM::NonLinearTask()
{
	double nonLinearErrMax = 1e-12;
	double nonLinearErr = nonLinearErrMax * 10;
	int maxiter = 1000;
	double minSolutionChange = 1e-10;
	double solutionChange = minSolutionChange * 10;
	for (int i = 0; i < maxiter && nonLinearErr > nonLinearErrMax && solutionChange > minSolutionChange; i++)
	{
		arrayspace::copy(qPrev, q, slae.A.n);
		slae.A.fillMatrix(0.0);
		newtAssembler.assembleGlobalMatrix(slae, mesh, q);
		addFirstConditions();
		slae.count_LOS(q, 1000, 1e-12);
		nonLinearErr = countDescr();
		solutionChange = countChange();
		cout << "Iteration " << i << " Descr = " << nonLinearErr << " Solution Change = " << solutionChange << endl;
	}
}

double FEM::countDescr()
{
	simpleAssembler.assembleGlobalMatrix(slae, mesh, q);
	slae.A.mult(q, bForCheck);
	arrayspace::minus(slae.b, bForCheck, descr, slae.A.n);
	return abs(scal(descr, descr, slae.A.n) / scal(slae.b, slae.b, slae.A.n));
}

double FEM::countChange()
{
	double* dif = qPrev;
	arrayspace::minus(q, qPrev, dif, slae.A.n);

	return abs(scal(dif, dif, slae.A.n) / scal(q, q, slae.A.n));
}

