#pragma once
#include "mesh.h"

namespace sparsematrix
{
	struct SparseMatrixSym;

	struct SparseMatrixPortrait
	{
		int* ig;
		int* jg;
		double* di;
		int n;
		SparseMatrixPortrait();
		SparseMatrixPortrait(SparseMatrixPortrait& portarait);
		double* operator * (double* v);
		void mult(double* v, double* res);
#ifdef MESH_IS_INCLUDE
		void buildPortrait(meshspace::Mesh mesh);
#endif
		virtual void fillMatrix(double mean) = 0;
	protected:
		void copy(SparseMatrixPortrait& M);
		virtual void multMatElemAndVecElem(int i, int j, int ij, double* v, double*res) = 0;

		

		//virtual void printInFullFormat() = 0;
	};

	struct SparseMatrixAsym : public SparseMatrixPortrait
	{
		double* ggl;
		double* ggu;

		SparseMatrixAsym();

		void copy(SparseMatrixAsym& M);

		void copy(SparseMatrixSym& M);

		int decomp_mat_LU(SparseMatrixAsym& LU);

		void fillMatrix(double mean) override;
	private:
		void multMatElemAndVecElem(int i, int j, int ij, double* v, double* res) override;
		void allocateMemoryForElems();
		
	};

	struct SparseMatrixSym : public SparseMatrixPortrait
	{
		double* gg;

		SparseMatrixSym();

		void copy(SparseMatrixSym& M);

		//operator SparseMatrixAsym() const;

		void fillMatrix(double mean) override;

	private:
		void multMatElemAndVecElem(int i, int j, int ij, double* v, double* res) override;
		void allocateMemoryForElems();
		
	};
}

