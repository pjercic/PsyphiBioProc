#ifndef IINPUT_H
#define IINPUT_H

#include <iostream>
#include <string>
using namespace std;

class IInput
{
public:
	// pure virtaul function
	virtual void Open(string fileName, int channel) = 0;
	virtual void Close() = 0;
	virtual double* getSignal(int& noSamples, int& smpFreq) = 0;
	virtual double getNextDataPoint() = 0;
};

#endif