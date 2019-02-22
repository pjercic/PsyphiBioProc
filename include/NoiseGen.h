#ifndef NOISE_GEN_H
#define NOISE_GEN_H

#include <iostream>
#include <IInput.h>

class NoiseGen : public IInput {
public:
	NoiseGen();
	double* getSignal(int& noSamples,  int& smpFreq);
	double getNextDataPoint();
	void Open(string fileName, int channel);
	void Close();

private:
	double*		rndSignalBuf;			// Memory for samples read (converted)
	int			sampleCounter;			// Counter for the individual signal data point
	int			noRndSamples;			// Number of created samples
};

#endif