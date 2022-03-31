#pragma once
#include <fstream>

using namespace std;
namespace arrayspace
{
	template<typename type>
	void copy(type* arReciver, type* arSource, unsigned int n);

	void copy(double* arReciver, double* arSource, unsigned int n);

	void copy(int* arReciver, int* arSource, unsigned int n);

	template<typename type>
	int read_vec(type* vec, int n, ifstream& in);

	int read_vec(double* vec, int n, ifstream& in);

	int read_vec(int* vec, int n, ifstream& in);

	template<typename type>
	int fill_vec(type* vec, int n, type elem);

	int fill_vec(double* vec, int n, double elem);

	int fill_vec(int* vec, int n, int elem);

	double scal(double* v1, double* v2, int n);

	void minus(double* v1, double* v2, double* res, int n);

}
