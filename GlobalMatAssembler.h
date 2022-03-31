#pragma once
#include "slae.h"
#include "mesh.h"
#include "FinalElemetsBasic.h"
using namespace slae;
using namespace meshspace;

class GlobalMatAssembler
{
public:
	void assembleGlobalMatrix(SLAE& A, Mesh& mesh, double* q);

protected:
	static const int N_KNOTS_1D = FinalElemBase::N_KNOTS_1D;
	static const int N_KNOTS_2D = FinalElemBase::N_KNOTS_2D;

	double localG[N_KNOTS_2D][N_KNOTS_2D], localM[N_KNOTS_2D][N_KNOTS_2D], localB[N_KNOTS_2D];

	void addLocalToGlobalAll(SLAE& slae, double G_local[N_KNOTS_2D][N_KNOTS_2D], double M_local[N_KNOTS_2D][N_KNOTS_2D], double localB[N_KNOTS_2D], int L[N_KNOTS_2D]);

	virtual void assembleLocal(FinalElemBase finalElemBase) = 0;
};

class NewtAssembler : public GlobalMatAssembler
{
	void assembleLocal(FinalElemBase finalElemBase) override;
};

class SimpleAssembler :public GlobalMatAssembler
{
	void assembleLocal(FinalElemBase finalElemBase) override;
};

