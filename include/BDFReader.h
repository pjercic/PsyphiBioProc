#ifndef BDF_READER_H
#define BDF_READER_H

#include <iostream>
#include <string>
#include <ctime>
#include "IInput.h"
#include "edflib.h"
using namespace std;

class BDFReader : public IInput {
public:
	BDFReader();
	double* getSignal(int& noSamples, int& smpFreq);
	double getNextDataPoint();
	void Open(string fileName, int channel);
	void Close();

private:
	struct edf_hdr_struct hdr;			// BDF/EDF header info structure
	struct tm* datetime;						// Initial timestamp
	float	samplingFreq;				// Sampling Frequency of the channel
	double *buf;						// Memory for samples read (converted)
	int		channelNo;					// Number of a channel 
	int		sampleCounter;			// Counter for the individual signal data point

	void PrintHeader (int nChannel);
	void ReadData(int nChannel);
};

#endif