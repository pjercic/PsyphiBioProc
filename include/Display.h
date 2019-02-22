#ifndef DISPLAY_H
#define DISPLAY_H

#include <iostream>
#include <string>
#include <deque>
#include "IOutput.h"
using namespace std;

class Display : public IOutput {
public:
	void saveSignal(double* sig, int noSamp, int smp_freq, string sigLabel, string sigDim);
	void saveNextDataPoint();
	void DisplayMessage(string message);
	void DisplayHistogram(deque<int> histogram);
};

#endif