#pragma once
#include <string>

using namespace std;

namespace fileparser
{
	const string FILEPARSER_OPEN_FILE_ERR = "не удалось открыть файл";
	const string PARAM_READ_ERR = "ошибка получения параметра";
	const string READ_ERR = "ошибка чтения";
	const string WRITE_ERR = "ошибка записи";

	static class fileParser
	{
	public:
		template <typename T>
		static bool getParam(string paramName, string fileName, T& param);

		static bool getParam(string paramName, string fileName, int& param);
		static bool getParam(string paramName, string fileName, double& param);
		static int getNumLinesInFile(string fileName);
	};
}

