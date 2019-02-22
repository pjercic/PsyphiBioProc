#include "ECG.h"
using namespace std;

void ECG::setInput(IInput* input) {
	_input = input;
}

void ECG::setOutput(IOutput* output) {
	_output = output;
}

void ECG::Process() {
	ProcesInput();
	ProcesSignal();
}

ECG::ECG() {
	_input = 0;
	_output = 0;
}

ECG::ECG(IInput* input, IOutput* output) {
	_input = input;
	_output = output;
}

void ECG::ProcesInput() {
	_input->getNextDataPoint();
	_output->saveNextDataPoint();
}

void ECG::ProcesSignal() {
	int noSamples, smpFreq;
	double* sig = _input->getSignal(noSamples, smpFreq);
	_output->saveSignal(sig, noSamples, smpFreq, "a", "b");
}