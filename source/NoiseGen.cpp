#include <ctime>
#include <cmath>
#include <cstdlib>
#include <NoiseGen.h>
using namespace std;

NoiseGen::NoiseGen () {
	sampleCounter = 0;
}

double* NoiseGen::getSignal(int& noSamples,  int& smpFreq) {

	noSamples = noRndSamples;
	smpFreq = noRndSamples;
	return rndSignalBuf;
}

double NoiseGen::getNextDataPoint() {

	return rndSignalBuf[sampleCounter++];
}

void NoiseGen::Open(string fileName, int channel) {

	srand((unsigned)time(0));

	noRndSamples = 204800;

	// Allocate memory for the converted sampes read
	rndSignalBuf=new double[noRndSamples];

	switch(channel) {
	case 0:
		for(int i=0; i<noRndSamples; i++)
			rndSignalBuf[i] = rand();
		break;

	case 1:
		for(int i=0; i<noRndSamples; i++)
			rndSignalBuf[i] = rand() + rand();
		break;

	case 2:
		for(int i=0; i<noRndSamples; i++) {
			rndSignalBuf[i] = 0;
			for(int j=0; j<12; j++)
				rndSignalBuf[i] += rand();
			rndSignalBuf[i] /= 12;
		}
		break;

	// A normally distributed random number,
	case 3:
		for(int i=0; i<noRndSamples; i++)
			rndSignalBuf[i] = powf(-2 * log10((float)rand()), 1 / 2) + cos(2 * 3.141592 * rand());
		break;
	}
	
}


void NoiseGen::Close() {

	delete [] rndSignalBuf;
}
