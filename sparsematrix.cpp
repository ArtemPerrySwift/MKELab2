#include "sparsematrix.h"
#include <utility>
#include "array.h"
#ifdef MESH_IS_INCLUDE
#include <set>
using namespace meshspace;
#endif
using namespace arrayspace;

namespace sparsematrix
{
	SparseMatrixPortrait::SparseMatrixPortrait()
	{
		n = 0;
		ig = NULL;
		jg = NULL;
		di = NULL;
	}

	SparseMatrixPortrait::SparseMatrixPortrait(SparseMatrixPortrait& portarait)
	{
		int* ig = portarait.ig;
		int* jg = portarait.jg;
		int n = portarait.n;

		this->di = new double[n];
		this->n = n;
		this->ig = new int[n + 1];
		arrayspace::copy(this->ig, ig, n + 1);

		int n_jg = ig[n];
		this->jg = new int[n_jg];
		arrayspace::copy(this->jg, jg, n_jg);

	}
	void SparseMatrixPortrait::copy(SparseMatrixPortrait& M)
	{
		n = M.n;
		ig = M.ig;
		jg = M.jg;
	}	
		//virtual void printInFullFormat() = 0;

	SparseMatrixAsym::SparseMatrixAsym() : SparseMatrixPortrait()
	{
		ggl = NULL;
		ggu = NULL;
	}

	void SparseMatrixAsym::allocateMemoryForElems()
	{
		ggl = new double[ig[n]];
		ggu = new double[ig[n]];
	}

	void SparseMatrixAsym::copy(SparseMatrixAsym& M)
	{
		SparseMatrixPortrait::copy(M);
		allocateMemoryForElems();

		arrayspace::copy(ggl, M.ggl, n);
		arrayspace::copy(ggu, M.ggu, n);
	}

	void SparseMatrixAsym::copy(SparseMatrixSym& M)
	{
		SparseMatrixPortrait::copy(M);
		allocateMemoryForElems();

		arrayspace::copy(ggl, M.gg, n);
		arrayspace::copy(ggu, M.gg, n);
	}


	SparseMatrixSym::SparseMatrixSym() : SparseMatrixPortrait()
	{
		gg = NULL;
	}

	void SparseMatrixSym::allocateMemoryForElems()
	{
		gg = new double[ig[n]];
	}

	void SparseMatrixSym::copy(SparseMatrixSym& M)
	{
		SparseMatrixPortrait::copy(M);
		allocateMemoryForElems();
		arrayspace::copy(gg, M.gg, n);
	}

	//SparseMatrixSym::operator SparseMatrixAsym() const
	//{
	//	SparseMatrixAsym M;
	//	M.copy(*this);
	//	return M;
	//}

	void SparseMatrixPortrait::mult(double* v, double* res)
	{
		int i, j, ij;
		int i_beg, i_end;

		for (i = 0; i < n; i++)
		{
			res[i] = 0;
			i_beg = ig[i];
			i_end = ig[i + 1];
			for (ij = i_beg; ij < i_end; ij++)
			{
				j = jg[ij];
				multMatElemAndVecElem(i, j, ij, v, res);
			}
			res[i] += di[i] * v[i];
		}
	}
	double* SparseMatrixPortrait:: operator * (double* v)
	{
		double* res = new double[n];
		mult(v, res);
		return res;
	}

	void SparseMatrixSym::multMatElemAndVecElem(int i, int j, int ij, double* v, double* res) 
	{
		res[i] += gg[ij] * v[j];
		res[j] += gg[ij] * v[i];
	}

	void SparseMatrixAsym::multMatElemAndVecElem(int i, int j, int ij, double* v, double* res)
	{
		res[i] += ggl[ij] * v[j];
		res[j] += ggu[ij] * v[i];
	}

	int SparseMatrixAsym::decomp_mat_LU(SparseMatrixAsym& LU)
	{

		double* LU_di = LU.di;
		double* LU_ggl = LU.ggl;
		double* LU_ggu = LU.ggu;
		int* LU_ig = LU.ig;
		int* LU_jg = LU.jg;

		int i, j, ij, ik, kj;
		int i_beg, i_end, j_beg, j_end;
		for (i = 0; i < n; i++)
		{
			i_beg = ig[i];
			i_end = ig[i + 1];
			LU_di[i] = di[i];
			for (ij = i_beg; ij < i_end; ij++)
			{
				j = jg[ij];
				j_beg = ig[j];
				j_end = ig[j + 1];
				LU_ggl[ij] = ggl[ij];
				LU_ggu[ij] = ggu[ij];
				for (ik = i_beg, kj = j_beg; kj < j_end && ik < ij;)
				{
					if (jg[ik] == jg[kj])
					{
						LU_ggl[ij] -= LU_ggl[ik] * LU_ggu[kj];
						LU_ggu[ij] -= LU_ggl[kj] * LU_ggu[ik];
						ik++;
						kj++;
					}
					else if (jg[ik] > jg[kj]) kj++;
					else ik++;
				}
				LU_ggu[ij] /= LU_di[j];
				LU_di[i] -= LU_ggl[ij] * LU_ggu[ij];
			}
		}
		return 0;
	}

	void SparseMatrixAsym::fillMatrix(double mean)
	{
		fill_vec(ggl, ig[n], mean);
		fill_vec(ggu, ig[n], mean);
	}

	void SparseMatrixSym::fillMatrix(double mean)
	{
		fill_vec(gg, ig[n], mean);
	}
#ifdef MESH_IS_INCLUDE
	void SparseMatrixPortrait::buildPortrait(meshspace::Mesh mesh)
	{
		set<int>* map;
		map = new set<int>[mesh.knotsStorage.kuzlov];
		int kuzlov = mesh.knotsStorage.kuzlov;
		int ktr = mesh.finalElementStorage.ktr;
		FinalElement* finalElements = mesh.finalElementStorage.endElements;

		int indexes[4];
		
		int i, j, k;
		for(k = 0; k < ktr; k++)
		{ 
			indexes[0] = finalElements[k].ver1;
			indexes[1] = finalElements[k].ver2;
			indexes[2] = finalElements[k].ver3;
			indexes[3] = finalElements[k].ver4;

			for (i = 0; i < 4; i++)
			{
				for (j = 0; j < 4; j++)
				{
					if (indexes[i] > indexes[j])
						map[indexes[i]].insert(indexes[j]);
				}
			}

			ig = new int[kuzlov + 1];
			ig[0] = 0;

			for (i = 0; i < ktr; i++)
			{
				ig[i + 1] = ig[i] + map[i].size();
			}
			jg = new int[ig[ktr]];

			int ijCount;
			int* elem;

			for (i = 0; i < ktr; i++)
			{
				j = ig[i];
				
				for (set<int>::iterator elem = map[i].begin(); elem != map[i].end(); elem++, j++)
					jg[j] = *elem;
				//for (auto& item : map[i])
				//{
				//	jg[j] = item;
				//	j++;
				//}
			}
		}
		
	}
#endif
}