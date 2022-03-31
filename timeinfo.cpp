#define _CRT_SECURE_NO_WARNINGS
#include <chrono>
#include <ctime>
#include "timeinfo.h"

string timeinfo::getCurrentTime()
{
	auto chronoCurrentTime = chrono::system_clock::now();
	time_t CurrentTime = chrono::system_clock::to_time_t(chronoCurrentTime);

	return ctime(&CurrentTime);
}