#ifndef IANALYSIS_H
#define IANALYSIS_H

#include <iostream>
#include "IInput.h"
#include "IOutput.h"

class IAnalysis
{
public:
	// pure virtaul function
	virtual void setInput(IInput* input) = 0;
    virtual void setOutput(IOutput* input) = 0;
    virtual void Process() = 0;
};

#endif