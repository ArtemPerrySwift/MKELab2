#include <iostream>
#include "mesh.h"
#include <string>
#include <locale>

using namespace std;
using namespace meshspace;

int main()
{
	
	string path = "D:\\Метод конечных элементов\\TELMA\\TELMA\\PRIMER\\";
	/*
	KnotsStorage knots;
	knots.read(path);

	FinalElementStorage endElements;
	endElements.read(path);

	KnotsWithFirstConditionStorage knotsWithFirstCondition;
	knotsWithFirstCondition.read(path);
	*/
	Mesh mesh(path);
	

	return 0;
}