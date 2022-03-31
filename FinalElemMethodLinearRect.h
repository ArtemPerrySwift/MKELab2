#pragma once
#include "sparsematrix.h"
#include "mesh.h"
#include "DifferentEquParams.h"
#include "FinalElemetsBasic.h"

using namespace sparsematrix;
using namespace meshspace;

class FinalElemMethodLinearRect
{
	SparseMatrixAsym A;
	Mesh P1;

	void buildPortrait();
	void buildLocalM(int iFinalElement, double localM[FinalElemBase::N_KNOTS_2D][FinalElemBase::N_KNOTS_2D]);
	void buildLocalG(int iFinalElement, double localG[FinalElemBase::N_KNOTS_2D][FinalElemBase::N_KNOTS_2D]);
	void buildLocalGNewt(int iFinalElement, double localG[FinalElemBase::N_KNOTS_2D][FinalElemBase::N_KNOTS_2D]);
public:
	void countLinearTask(double*& q);
	void countNonLinearTask(double*& q);
	void addLocalMatrToGlobal();
};



