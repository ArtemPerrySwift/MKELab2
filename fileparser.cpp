#include <fstream>
#include "fileparser.h"
#include "programlog.h"

using namespace programlog;

namespace fileparser
{
	template <typename T>
	bool fileParser::getParam(string paramName, string fileName, T& param)
	{
		ifstream paramFile;
		paramFile.open(fileName);

		if (!paramFile.is_open())
		{
			writeErr("Ошибка открытия файла - " + fileName);
			return false;
		}


		string buf;
		for (paramFile >> buf; buf.find(paramName) != string::npos && !paramFile.eof(); paramFile >> buf);
		if (buf.find(paramName) == string::npos)
		{
			paramFile >> param;
			return true;
		}

		paramFile.close();
		return false;

	}

	bool fileParser::getParam(string paramName, string fileName, int& param)
	{
		ifstream paramFile;
		paramFile.open(fileName);

		if (!paramFile.is_open())
		{
			writeErr("Ошибка открытия файла - " + fileName);
			return false;
		}


		string buf;
		for (paramFile >> buf; buf.find(paramName) != string::npos && !paramFile.eof(); paramFile >> buf);
		if (buf.find(paramName) == string::npos)
		{
			paramFile >> param;
			return true;
		}

		paramFile.close();
		return false;

	}

	bool fileParser::getParam(string paramName, string fileName, double& param)
	{
		ifstream paramFile;
		paramFile.open(fileName);

		if (!paramFile.is_open())
		{
			writeErr("Ошибка открытия файла - " + fileName);
			return false;
		}


		string buf;
		for (paramFile >> buf; buf.find(paramName) != string::npos && !paramFile.eof(); paramFile >> buf);
		if (buf.find(paramName) == string::npos)
		{
			paramFile >> param;
			return true;
		}

		paramFile.close();
		return false;

	}
	int fileParser::getNumLinesInFile(string fileName)
	{
		int number_of_lines = 0;
		string line;
		ifstream fileForCount(fileName);

		if (!fileForCount.is_open())
		{
			writeErr("Ошибка открытия файла - " + fileName);
			return -1;
		}

		while (getline(fileForCount, line))
			++number_of_lines;

		fileForCount.close();
		return number_of_lines;
	}
}
