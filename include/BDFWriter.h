#ifndef BDF_WRITER_H
#define BDF_WRITER_H

#include <iostream>
#include <string>
#include <deque>
#include "IOutput.h"
#include <edflib.h>
using namespace std;

class BDFWriter : public IOutput {
public:
	BDFWriter();
	void saveSignal(double* sig, int noSamp, int smp_freq, string sigLabel, string sigDim);
	void saveNextDataPoint();
	void DisplayMessage(string message);
	void DisplayHistogram(deque<int> histogram);

private:
	deque<double*>	signals;
	deque<int>		smpFreqs;
	deque<int>		nosSamples;
	deque<string>		sigLabels;
	deque<string>		sigDim;
};

#endif