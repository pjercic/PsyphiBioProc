#include <iostream>
#include "IInput.h"
#include "IOutput.h"
#include "IAnalysis.h"
#include "BDFReader.h"
#include "NoiseGen.h"
#include "ECG.h"
#include "DSPStatistic.h"
#include "Display.h"
#include <BDFWriter.h>
using namespace std;

int main (int argc, char *argv[]) {

	IInput* IObjInput = new BDFReader();
	//IInput* IObjInput = new NoiseGen();
	IObjInput->Open(argv[1], atoi(argv[2]));
	//IOutput* IObjOutput = new Display();
	IOutput* IObjOutput = new BDFWriter();

	// Interface Dependency Insertion
	IAnalysis* IObjAnalysis = new DSPStatistic();
	IObjAnalysis->setInput(IObjInput);
	IObjAnalysis->setOutput(IObjOutput);
	IObjAnalysis->Process();

	delete IObjAnalysis;
	delete IObjOutput;
	delete IObjInput;

	return 0;
}