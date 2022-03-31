#pragma once
#include <string>

using namespace std;

namespace programlog
{
	void writeErrInLogFile(string error);

	void writeErrInConcole(string error);

	void writeErr(string error);
}
