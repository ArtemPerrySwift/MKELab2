#include <fstream>
#include <iostream>
#include "programlog.h"
#include "timeinfo.h"

const string LogFileName = "ProgramLog.txt";

namespace programlog
{
	void writeErrInLogFile(string error)
	{
		ofstream errFile;
		errFile.open(LogFileName, ios::app);

		errFile << timeinfo::getCurrentTime() << " [Error] - " << error << endl;
	}

	void writeErrInConcole(string error)
	{
		cout << timeinfo::getCurrentTime() << " [Error] - " << error << endl;
	}

	void writeErr(string error)
	{
		writeErrInLogFile(error);
		writeErrInConcole(error);
	}
}
