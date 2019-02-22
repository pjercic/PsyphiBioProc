#include "Display.h"
using namespace std;

void Display::saveSignal(double* sig, int noSamp, int smp_freq, string sigLabel, string sigDim) {
	cout << "Console:saveSignal" << endl;
}

void Display::saveNextDataPoint() {
	cout << "Console:saveNextDataPoint" << endl;
}

void Display::DisplayMessage(string message) {
	cout << message << endl;
}

void Display::DisplayHistogram(deque<int> histogram) {

	int height = 50, width = 30;
	int sumValues = 0;
	float max = 0;

	float* probabilityDensityFcn = new float[histogram.size()];

	for(int i = 0; i < histogram.size(); i++)
		sumValues += histogram[i];

	for(int i = 0; i < histogram.size(); i++) {

		probabilityDensityFcn[i] = float(histogram[i]) / sumValues;

		if(probabilityDensityFcn[i] > max)
			max = probabilityDensityFcn[i];
	}

	for(int i = 0; i < histogram.size(); i += histogram.size() / width) {
		int noOfStars = probabilityDensityFcn[i] / max * height;
		for(int j = 0; j < noOfStars; j++)
			cout << "*";
		cout << endl;
	}
}