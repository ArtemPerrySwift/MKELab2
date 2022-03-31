#include "array.h"
#include <cstring>
#include <fstream>

using namespace std;
namespace arrayspace
{
	template<typename type>
	void copy(type* arReciver, type* arSource, unsigned int n)
	{
		size_t generalArrraySize = n * sizeof(type); // Размер массивов в байтах
		memcpy(arReciver, arSource, generalArrraySize); // Копируем сразу весь массив целиком
	}

	void copy(double* arReciver, double* arSource, unsigned int n)
	{
		size_t generalArrraySize = n * sizeof(double); // Размер массивов в байтах
		memcpy(arReciver, arSource, generalArrraySize); // Копируем сразу весь массив целиком
	}

	void copy(int* arReciver, int* arSource, unsigned int n)
	{
		size_t generalArrraySize = n * sizeof(int); // Размер массивов в байтах
		memcpy(arReciver, arSource, generalArrraySize); // Копируем сразу весь массив целиком
	}

	template<typename type>
	int read_vec(type* vec, int n, ifstream& in)
	{
		for (int i = 0; i < n; i++)
			in >> vec[i];
		return 0;
	}

	int read_vec(double* vec, int n, ifstream& in)
	{
		for (int i = 0; i < n; i++)
			in >> vec[i];
		return 0;
	}

	int read_vec(int* vec, int n, ifstream& in)
	{
		for (int i = 0; i < n; i++)
			in >> vec[i];
		return 0;
	}

	template<typename type>
	int fill_vec(type* vec, int n, type elem)
	{
		for (int i = 0; i < n; i++)
			vec[i] = elem;
		return 0;
	}

	int fill_vec(double* vec, int n, double elem)
	{
		for (int i = 0; i < n; i++)
			vec[i] = elem;
		return 0;
	}

	int fill_vec(int* vec, int n, int elem)
	{
		for (int i = 0; i < n; i++)
			vec[i] = elem;
		return 0;
	}

	double scal(double* v1, double* v2, int n)
	{
		double sum = 0;
		for (int i = 0; i < n; i++)
			sum += v1[i]*v2[i];
		return 0;
	}

	void minus(double* v1, double* v2, double* res, int n)
	{
		for (int i = 0; i < n; i++)
			res[i] = v1[i] - v2[i];

	}
	
	
}