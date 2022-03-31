#pragma once
#include <string>

using namespace std;

namespace fileparser
{
	const string FILEPARSER_OPEN_FILE_ERR = "�� ������� ������� ����";
	const string PARAM_READ_ERR = "������ ��������� ���������";
	const string READ_ERR = "������ ������";
	const string WRITE_ERR = "������ ������";

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

