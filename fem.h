#pragma once
#include "slae.h"
#include "mesh.h"
#include "FinalElemetsBasic.h"
#include "GlobalMatAssembler.h"

using namespace meshspace;
using namespace slae;

class FEM
{
	SLAE slae;
	Mesh mesh;
	NewtAssembler newtAssembler;
	SimpleAssembler simpleAssembler;

	double* q, *qPrev;
	double* descr;
	double* bForCheck;
	void addFirstConditions();
	void LinearTask();
	void NonLinearTask();
	void init(Mesh mesh);

	/// <summary>
	/// Расчёт невязки
	/// </summary>
	double countDescr();
	double countChange();

public:
	double solution(double x, double y);
};

