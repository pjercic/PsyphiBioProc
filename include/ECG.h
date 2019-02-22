#ifndef ECG_H
#define ECG_H

#include <iostream>
#include "IInput.h"
#include "IOutput.h"
#include "IAnalysis.h"

class ECG : public IAnalysis {
public:
	void setInput(IInput* input);

	void setOutput(IOutput* output);

	void Process();

	ECG();

	ECG(IInput*input, IOutput* output);

private:
	IInput* _input;
	IOutput* _output;

	void ProcesInput();

	void ProcesSignal();
};

#endif