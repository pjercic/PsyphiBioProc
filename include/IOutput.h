#ifndef IOUTPUT_H
#define IOUTPUT_H

#include <iostream>
#include <string>
#include <deque>
using namespace std;

class IOutput
{
public:
	// pure virtaul function
	virtual void saveSignal(double* sig, int noSamp, int smp_freq, string sigLabel, string sigDim) = 0;
    virtual void saveNextDataPoint() = 0;
	virtual void DisplayMessage(string message) = 0;
	virtual void DisplayHistogram(deque<int> histogram) = 0;
};

#endif