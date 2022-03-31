#include "GlobalMatAssembler.h"
void GlobalMatAssembler::assembleGlobalMatrix(SLAE& slae, Mesh& mesh, double* q)
{
	FinalElement* finalElements = mesh.finalElementStorage.endElements;
	int ktr = mesh.finalElementStorage.ktr;
	int materialInd;
	FinalElemBase finalElemBase;
	Material material;
	slae.A.fillMatrix(0);
	for (int i = 0; i < ktr; i++)
	{
		materialInd = mesh.finalElementMaterialStorage.finalElemMaterials[i];
		material = mesh.materialStorage.findMaterialByNum(materialInd);
		finalElemBase.setFinalElement(finalElements[i], mesh.knotsStorage, material, q);
		finalElemBase.buildNewt(localG, localB);
		finalElemBase.buildLocalM(localM, 0);
		addLocalToGlobalAll(slae, localG, localM, localB, finalElemBase.globalInd);
	}
}

void GlobalMatAssembler::addLocalToGlobalAll(SLAE& slae, double G_local[N_KNOTS_2D][N_KNOTS_2D], double M_local[N_KNOTS_2D][N_KNOTS_2D], double localB[N_KNOTS_2D], int L[N_KNOTS_2D])
{
	int m, l;
	int ind;
	int k_left, k_right;
	double* di = slae.A.di;
	double* ggl = slae.A.ggl;
	double* ggu = slae.A.ggu;
	int* ig = slae.A.ig;
	int* jg = slae.A.jg;
	double* b = slae.b;
	for (m = 0; m < N_KNOTS_2D; m++)
	{
		di[L[m]] += G_local[m][m] + M_local[m][m];
		b[L[m]] += localB[m];
	}

	for (m = 1; m < N_KNOTS_2D; m++)
	{
		k_left = ig[L[m]];
		for (l = 0; l < m; l++)
		{
			k_right = ig[L[m] + 1];
			while (jg[k_left] != L[l])
			{
				ind = (k_left + k_right) / 2; // djpvj;yj
				if (jg[ind] <= L[l])
				{
					k_left = ind;
				}
				else
				{
					k_right = ind;
				}
			}

			ggl[k_left] += G_local[m][l] + M_local[m][l];
			ggu[k_left] += G_local[l][m] + M_local[l][m];
			k_left++;
		}
	}
}

void NewtAssembler::assembleLocal(FinalElemBase finalElemBase)
{
	finalElemBase.buildNewt(localG, localB);
	finalElemBase.buildLocalM(localM, 0);
}

void SimpleAssembler::assembleLocal(FinalElemBase finalElemBase)
{
	finalElemBase.buildLocalG(localG);
	finalElemBase.buildLocalM(localM, 0);
	finalElemBase.buildLocalB(localB);
}