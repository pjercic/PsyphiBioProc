/**
* @file DSPStatistic.cpp
* @brief This source file defines an Analysis component using methods of DSP on signals.
*
* @author Petar Jerčić
*
* @date 09/12/2012
*/

#include <sstream>
#include <cmath>
#include <ctime>
#include <deque>
#include "DSPStatistic.h"

using namespace std;

void DSPStatistic::setInput(IInput* input) {
	_input = input;
}

void DSPStatistic::setOutput(IOutput* output) {
	_output = output;
}

void DSPStatistic::Process() {
	ProcesSignal();
}

DSPStatistic::DSPStatistic() {
	_input = 0;
	_output = 0;
}

DSPStatistic::DSPStatistic(IInput* input, IOutput* output) {
	_input = input;
	_output = output;
}

void DSPStatistic::ProcesInput() {
	_input->getNextDataPoint();
	_output->saveNextDataPoint();
}

void DSPStatistic::ProcesSignal() {

	ostringstream sstream;

	signal = _input->getSignal(noSamples, smpFreq);

	double avg = 0;
	double var = 0;
	double stdev = 0;
	double snr = 0;
	double cv = 0;

	CalculateDescriptivesStatic(signal, noSamples, avg, stdev, var, snr, cv);

	sstream << "\nStatic\n\nAvg\t";
	sstream << avg;
	sstream << "\nStdev\t";
	sstream << stdev;
	sstream << "\nVar\t";
	sstream << var;
	sstream << "\nSNR\t";
	sstream << snr;
	sstream << "\nCV\t";
	sstream << cv;

	//_output->saveSignal(signal, noSamples, smpFreq, "ECG", "mV");
	//_output->DisplayMessage(sstream.str());

	sstream.str("");
	sstream.clear();

	for(int i = 0; i < noSamples; i++)
		CalculateDescriptivesRunning(_input->getNextDataPoint(), avg, stdev, var, snr, cv);

	sstream << "\nRunning\n\nAvg\t";
	sstream << avg;
	sstream << "\nStdev\t";
	sstream << stdev;
	sstream << "\nVar\t";
	sstream << var;
	sstream << "\nSNR\t";
	sstream << snr;
	sstream << "\nCV\t";
	sstream << cv;

	//_output->DisplayMessage(sstream.str());

	sstream.str("");
	sstream.clear();

	deque<int> histogram = CalculateDescriptivesStaticHistogram(signal, noSamples, avg, stdev, var, snr, cv);

	sstream << "\nHistogram\n\nAvg\t";
	sstream << avg;
	sstream << "\nStdev\t";
	sstream << stdev;
	sstream << "\nVar\t";
	sstream << var;
	sstream << "\nSNR\t";
	sstream << snr;
	sstream << "\nCV\t";
	sstream << cv;

	//_output->DisplayMessage(sstream.str());
	//_output->DisplayHistogram(histogram);

	sstream.str("");
	sstream.clear();

	double* ir = new double[3];
	ir[0] = 0; ir[1] = 1; ir[2] = 0;

	double* ors1 = NULL;
	int orSize1;

	ConvolutionInputSide(signal, noSamples, ir, 3, &ors1, orSize1);

	CalculateDescriptivesStatic(ors1, orSize1, avg, stdev, var, snr, cv);

	sstream << "\nConvolution Input\n\nAvg\t";
	sstream << avg;
	sstream << "\nStdev\t";
	sstream << stdev;
	sstream << "\nVar\t";
	sstream << var;
	sstream << "\nSNR\t";
	sstream << snr;
	sstream << "\nCV\t";
	sstream << cv;

	//output->saveSignal(ors1, orSize1, smpFreq, "ECGconv1", "mV");
	//_output->DisplayMessage(sstream.str());

	sstream.str("");
	sstream.clear();

	double* ors2 = NULL;
	int orSize2;

	ConvolutionOutputSide(signal, noSamples, ir, 3, &ors2, orSize2);

	CalculateDescriptivesStatic(ors2, orSize2, avg, stdev, var, snr, cv);

	sstream << "\nConvolution Output\n\nAvg\t";
	sstream << avg;
	sstream << "\nStdev\t";
	sstream << stdev;
	sstream << "\nVar\t";
	sstream << var;
	sstream << "\nSNR\t";
	sstream << snr;
	sstream << "\nCV\t";
	sstream << cv;

	//_output->saveSignal(ors2, orSize2, smpFreq, "ECGconv2", "mV");
	//_output->DisplayMessage(sstream.str());

	sstream.str("");
	sstream.clear();

	double* ors3 = NULL;
	int orSize3;

	FirstDifference(signal, noSamples, &ors3, orSize3);

	//_output->saveSignal(ors3, orSize3, smpFreq, "FirstDerivate", "dmV");
	//_output->DisplayMessage(sstream.str());

	double* relX = NULL;
	int relXSize;
	double* imgX = NULL;
	int imgXSize;

	DiscreteFourierTransformationTimeSide(signal, noSamples, &relX, relXSize, &imgX, imgXSize);
	//_output->saveSignal(relX, relXSize, smpFreq, "RealX", "freq");
	//_output->saveSignal(imgX, imgXSize, smpFreq, "ImgX", "freq");

	double* ors4 = NULL;
	int orSize4;

	InverseDiscreteFourierTransformationTimeSide(relX, relXSize, imgX, imgXSize, &ors4, orSize4);

	//_output->saveSignal(ors4, orSize4, smpFreq, "Recreated", "mV");
	//_output->DisplayMessage(sstream.str());

	double* mag = NULL;
	int magSize;
	double* phase = NULL;
	int phaseSize;
	double* uphase = NULL;
	int uphaseSize;

	RectangularToPolarConversion(relX, relXSize, imgX, imgXSize, &mag, magSize, &phase, phaseSize);
	PhaseUnwrapping(phase, phaseSize, &uphase, uphaseSize);
	//_output->saveSignal(mag, magSize, smpFreq, "Magnitude", "freq");
	//_output->saveSignal(uphase, uphaseSize, smpFreq, "Phase", "freq");


	/*double* y = NULL;

	y = MovingAverageFilter(signal, noSamples);
	_output->saveSignal(signal, noSamples, smpFreq, "Original", "mV");
	_output->saveSignal(y, noSamples, smpFreq, "MvgAvgFilter", "mV");

	double* y_rec = NULL;

	y_rec = MovingAverageFilterRecursion(signal, noSamples);
	_output->saveSignal(y_rec, noSamples, smpFreq, "MvgAvgFilterRec", "mV");

	_output->DisplayMessage(sstream.str());*/

	double* y = NULL;

	y = SinglePoleLowPassRecursiveFilter(signal, noSamples, 0.05);
	_output->saveSignal(signal, noSamples, smpFreq, "Original", "mV");
	_output->saveSignal(y, noSamples, smpFreq, "LowPass", "mV");

	_output->DisplayMessage(sstream.str());
}

/**
* Descriptives
* @author Petar Jerčić
* @date 26/02/2013
* @TODO N/A
*/
void	DSPStatistic::CalculateDescriptivesStatic(double* sig, int noSamp, double& avg, double& stdev, double& var, double &snr, double &cv) {

	avg = 0;						// Average value
	var = 0;						// Variance value
	stdev = 0;						// Standard deviance value
	snr = 0;						// Signal-to-noise ratio
	cv = 0;							// Coefficient of Variation

	// Calculate Average
	for(int i = 0; i < noSamp; i++)
		avg += sig[i];

	avg /= noSamp;

	// Calculate Variance
	for(int i = 0; i < noSamp; i++)
		var += pow(sig[i] - avg, 2);

	var /= noSamp - 1;

	// Calculate Standard Deviation
	stdev = sqrt(var);

	// Calculate Signal-to-noise ratio
	snr = avg / stdev;

	// Calculate Coefficient of Variation
	cv = stdev / avg * 100;
}

void	DSPStatistic::CalculateDescriptivesRunning(double data, double& avg, double& stdev, double& var, double &snr, double &cv) {

	// Calculates estimates since we are estimating statistic of the underlying signal generation process
	static int		n = 0;
	static double	sum = 0;
	static double	sumSquares = 0;

	n++;
	sum += data;
	sumSquares += pow(data, 2);

	avg = sum/n;
	if(n>1)	stdev = sqrt((sumSquares - pow(sum, 2)/n)/(n-1));
	var = pow(stdev, 2);

	// Calculate Signal-to-noise ratio
	snr = avg / stdev;

	// Calculate Coefficient of Variation
	cv = stdev / avg * 100;
}

deque<int>	DSPStatistic::CalculateDescriptivesStaticHistogram(double* sig, int noSamp, double& avg, double& stdev, double& var, double &snr, double &cv){

	// Since floating point data binning up to a whole integer is used

	deque<int> histogram;	// Histogram
	int zeroPosition = 0;	// Position of 0 value

	avg = 0;						// Average value
	var = 0;						// Variance value
	stdev = 0;						// Standard deviance value
	snr = 0;						// Signal-to-noise ratio
	cv = 0;							// Coefficient of Variation

	// Calculate Histogram
	for(int i = 0; i < noSamp; i++) {
		int nData = sig[i];	// Binning formula - up to the nearest integer
		if(nData >= 0) {
			if(nData > histogram.size() - zeroPosition) {
				for(int j = histogram.size(); j < nData + zeroPosition; j++)
					histogram.push_back(0);
				histogram.push_back(1);
			} else {
				if(nData + zeroPosition >= histogram.size() || nData + zeroPosition < 0)
					histogram.push_back(0);
				histogram[nData + zeroPosition]++;
			}

		} else {
			if(nData < -zeroPosition) {
				for(int j = -zeroPosition; j > nData + 1; j--)
					histogram.push_front(0);
				histogram.push_front(1);
				zeroPosition = abs(nData);
			} else {
				histogram[abs(zeroPosition + nData)]++;
			}
		}
	}

	// Calculate Average
	for(int i = 0; i < histogram.size(); i++)
		avg += (i - zeroPosition) * histogram[i];
	avg /= noSamp + 1;

	// Calculate Variance
	for(int i = 0; i < histogram.size(); i++)
		var += histogram[i] * pow((i - zeroPosition) - avg, 2);
	var /= noSamp;

	// Calculate Standard Deviation
	stdev = sqrt(var);

	// Calculate sig-to-noise ratio
	snr = avg / stdev;

	// Calculate Coefficient of Variation
	cv = stdev / avg * 100;

	return histogram;
}

void	DSPStatistic::ConvolutionInputSide(double* sig, int noSamp, double* impulseResponse, int irSize, double** outputSignal, int& osSize) {
	
	// Output Signal size = Signal + Impulse response size
	osSize = noSamp + irSize;

	// Allocate the output Signal
	double* _outputSignal = new double[osSize];

	// Zero the output array
	for(int i = 0; i < osSize; i++)
		_outputSignal[i] = 0;

	// Loop over each point in impulse resp for each point in Signal
	for(int i = 0; i < noSamp; i++)
		for(int j = 0; j < irSize; j++)
			_outputSignal[i + j] += sig[i] * impulseResponse[j];

	// Store the output signal
	if(*outputSignal != NULL) delete[] *outputSignal;
	*outputSignal = _outputSignal;
}

void	DSPStatistic::ConvolutionOutputSide(double* sig, int noSamp, double* impulseResponse, int irSize, double** outputSignal, int& osSize) {
	
	// Output signal size = signl + Impulse response size
	osSize = noSamp + irSize;

	// Allocate the output Signal
	double* _outputSignal = new double[osSize];

	// Loop for each point of y[] output signal
	for(int i = 0; i < osSize; i++) {
		_outputSignal[i] = 0;

		// Loop for each point of h[] Impulse response
		for(int j = 0; j < irSize; j++) {
			if(i - j >= 0 && i - j <= noSamp)
				_outputSignal[i] +=  impulseResponse[j] * sig[i - j];
		}
	}

	// Loop over each point in impulse resp for each point in signal
	for(int i = 0; i < noSamp; i++)
		for(int j = 0; j < irSize; j++)
			_outputSignal[i + j] += sig[i] * impulseResponse[j];

	// Store the output signal
	if(*outputSignal != NULL) delete[] *outputSignal;
	*outputSignal = _outputSignal;
}

void	DSPStatistic::FirstDifference(double* sig, int noSamp, double** outputSignal, int& osSize) {

	// Output signal size = signl
	osSize = noSamp;

	// Allocate the output Signal
	double* _outputSignal = new double[osSize];

	// Edge holds no data
	_outputSignal[0] = 0;

	// Loop for each point of y[] output signal
	for(int i = 1; i < osSize; i++)
		_outputSignal[i] = sig[i] - sig[i - 1];

	// Store the output signal
	if(*outputSignal != NULL) delete[] *outputSignal;
	*outputSignal = _outputSignal;
}

void	DSPStatistic::RunningSumIntegration(double* sig, int noSamp, double** outputSignal, int& osSize) {

	// Output signal size = signl
	osSize = noSamp;

	// Allocate the output Signal
	double* _outputSignal = new double[osSize];

	// Edge holds no data
	_outputSignal[0] = sig[0];

	// Loop for each point of y[] output signal
	for(int i = 1; i < osSize; i++)
		_outputSignal[i] = _outputSignal[i - 1] + sig[i];

	// Store the output signal
	if(*outputSignal != NULL) delete[] *outputSignal;
	*outputSignal = _outputSignal;
}

/**
* Decompose a time domain signal into sinusoids in frequency domain signals from time side.
* Describes how an individual sample in the frequency domain is affected by all of the samples in the time domain. That is, the program calculates each of the values in the frequency domain in succession, not as a group. 
* @author Petar Jerčić
* @date 09/12/2012
* @TODO Check correct sizes of the pek anh imk arrays.
*/
void	DSPStatistic::DiscreteFourierTransformationTimeSide(double* timeSignal, int timeSigSize, double** realX, int& reXSize, double** imaginaryX, int& imXSize){

	double PI = 4 * atan(1.0);

	/// Alocate array for real and imaginary parts of frequency domain
	reXSize = timeSigSize / 2 + 1;
	imXSize = timeSigSize / 2 + 1;
	double* _realX = new double[reXSize];
	double* _imaginary = new double[imXSize];

	/// Zero out the \a realX and \a imaginaryX so it can be used as accumulators
	for(int i = 0; i < reXSize; i++) {
		_realX[i] = 0;
		_imaginary[i] = 0;
	}

	/// Decomposition method
	/// Correlate \a timeSignal with cosine and sine waves
	for(int i = 0; i < reXSize; i++)
		for(int j = 0; j < timeSigSize; j++) {
			_realX[i] += timeSignal[j] * cos(2 * PI * i * j / timeSigSize);
			_imaginary[i] += - timeSignal[j] * sin(2 * PI * i * j / timeSigSize);
		}

	// Store the output signal
	*realX = _realX;
	*imaginaryX = _imaginary;
}

/**
* Decompose a time domain signal into sinusoids in frequency domain signals from time side.
* The program loops through each sample in the time domain, calculating the contribution of that point to the frequency domain. The overall frequency domain is found by adding the contributions from the individual time domain points.
* @author Petar Jerčić
* @date 09/12/2012
* @TODO Check correct sizes of the pek anh imk arrays.
*/
void	DSPStatistic::DiscreteFourierTransformationFreqSide(double* timeSignal, int timeSigSize, double** realX, int& reXSize, double** imaginaryX, int& imXSize){

	double PI = 4 * atan(1.0);

	/// Alocate array for real and imaginary parts of frequency domain
	reXSize = timeSigSize / 2 + 1;
	imXSize = timeSigSize / 2 + 1;
	double* _realX = new double[reXSize];
	double* _imaginary = new double[imXSize];

	/// Zero out the \a realX and \a imaginaryX so it can be used as accumulators
	for(int i = 0; i < reXSize; i++) {
		_realX[i] = 0;
		_imaginary[i] = 0;
	}

	/// Decomposition method
	/// loops through each sample in the time domain \a timeSignal, calculating the contribution of that point to the frequency domain.
	for(int i = 0; i < reXSize; i++)
		for(int j = 0; j < timeSigSize; j++) {
		
			_realX[i] += timeSignal[j] * cos(2 * PI * i * j / timeSigSize);
			_imaginary[i] += - timeSignal[j] * sin(2 * PI * i * j / timeSigSize);
		}

	// Store the output signals
	*realX = _realX;
	*imaginaryX = _imaginary;
}

///Synthesis of time domain signal from frequency domain signals from frequency side.
///Each of the scaled sinusoids are generated one at a time and added to an accumulation array, which ends up becoming the time domain signal.
///@author Petar Jerčić
///@date 09/12/2012
///@TODO Check correct sizes of the pek anh imk arrays.
void	DSPStatistic::InverseDiscreteFourierTransformationFreqSide(double* realX, int reXSize, double* imaginaryX, int imXSize, double** timeSignal, int& timeSigSize) {

	double PI = 4 * atan(1.0);

	/// Alocate array for time domain signal (minus two extra data points)
	timeSigSize = 2 * reXSize - 2;
	double* _timeSignal = new double[timeSigSize];

	/// Alocate array for cosine and sine wave amplitudes
	double* realXAmplitudes = new double[reXSize];
	double* imaginaryXAmplitudes = new double[imXSize];

	/// Find the cosine and sine wave amplitudes
	for(int i = 0; i < reXSize; i++) {

		realXAmplitudes[i] = realX[i] / (timeSigSize / 2);
		imaginaryXAmplitudes[i] = - imaginaryX[i] / (timeSigSize / 2);
	}

	/// Cover two special cases
	realXAmplitudes[0] = realX[0] / 2;
	realXAmplitudes[reXSize - 1] = realX[reXSize - 1] / 2;

	/// Zero out the \a timeSignal so it can be used as accumulator
	for(int i = 0; i < timeSigSize; i++)
		_timeSignal[i] = 0;

	/// Synthesis method
	/// Loop through each frequency generating the entire length of the sine and cosine waves, and add them to the accumulator signal \a timeSignal
	for(int i = 0; i < reXSize; i++)
		for(int j = 0; j < timeSigSize; j++) {

			_timeSignal[j] += realXAmplitudes[i] * cos(2 * PI * i * j / timeSigSize);
			_timeSignal[j] += imaginaryXAmplitudes[i] * sin(2 * PI * i * j / timeSigSize);
		}

	// Store the output signal
	*timeSignal = _timeSignal;
}

/**
* Synthesis of time domain signal from frequency domain signals from time side.
* Each sample in the time domain signal is calculated one at a time, as the sum of all the corresponding samples in the cosine and sine waves.
* @author Petar Jerčić
* @date 09/12/2012
* @TODO N/A.
*/
void	DSPStatistic::InverseDiscreteFourierTransformationTimeSide(double* realX, int reXSize, double* imaginaryX, int imXSize, double** timeSignal, int& timeSigSize) {

	double PI = 4 * atan(1.0);

	/// Alocate array for time domain signal (minus two extra data points)
	timeSigSize = 2 * reXSize - 2;
	double* _timeSignal = new double[timeSigSize];

	/// Alocate array for cosine and sine wave amplitudes
	double* realXAmplitudes = new double[reXSize];
	double* imaginaryXAmplitudes = new double[imXSize];

	/// Find the cosine and sine wave amplitudes
	for(int i = 0; i < reXSize; i++) {

		realXAmplitudes[i] = realX[i] / (timeSigSize / 2);
		imaginaryXAmplitudes[i] = - imaginaryX[i] / (timeSigSize / 2);
	}

	/// Cover two special cases
	realXAmplitudes[0] = realX[0] / 2;
	realXAmplitudes[reXSize - 1] = realX[reXSize - 1] / 2;

	/// Zero out the \a timeSignal so it can be used as accumulator
	for(int i = 0; i < timeSigSize; i++)
		_timeSignal[i] = 0;

	/// Synthesis method
	/// Loop through each sample in the time domain and sum the corresponding samples from each cosine and sine waves.
	for(int i = 0; i < reXSize; i++)
		for(int j = 0; j < timeSigSize; j++) {
			_timeSignal[j] += realXAmplitudes[i] * cos(2 * PI * i * j / timeSigSize);
			_timeSignal[j] += imaginaryXAmplitudes[i] * sin(2 * PI * i * j / timeSigSize);
		}

	// Store the output signal
	*timeSignal = _timeSignal;
}

/**
* Rectangular to polar Conversion.
* From real and imginary to magnitude and phase.
* @author Petar Jerčić
* @date 09/12/2012
* @TODO N/A.
*/
void	DSPStatistic::RectangularToPolarConversion(double* realX, int reXSize, double* imaginaryX, int imXSize, double** mag, int& magSize, double** phase, int& phaseSize) {

	double PI = 4 * atan(1.0);

	/// Alocate array for magnitude and phase parts of frequency domain
	magSize = reXSize;
	phaseSize = imXSize;
	double* _mag = new double[reXSize];
	double* _phase = new double[imXSize];

	/// Rectangular to polar conversion
	for(int i = 0; i < reXSize; i++) {

		_mag[i] = sqrt(pow(realX[i], 2) + pow(imaginaryX[i], 2));

		/// Prevent divide by zero
		if(realX[i] = 0)	realX[i] = pow(1.0, -20);

		_phase[i] = atan(imaginaryX[i] / realX[i]);

		/// Correct arcus tangens ambiguity
		if(realX[i] < 0 && imaginaryX[i] < 0)	_phase[i] -= PI;
		if(realX[i] < 0 && imaginaryX[i] >= 0)	_phase[i] += PI;
	}

	// Store the output signals
	if(*mag != NULL) delete[] *mag;
	*mag = _mag;
	if(*phase != NULL) delete[] *phase;
	*phase = _phase;
}

/**
* Polar To Rectangular Conversion
* From magnitude and phase to real and imginary
* @author Petar Jerčić
* @date 09/12/2012
* @TODO N/A.
*/
void	DSPStatistic::PolarToRectangularConversion(double* mag, int magSize, double* phase, int phaseSize, double** realX, int& reXSize, double** imaginaryX, int& imXSize) {

	double PI = 4 * atan(1.0);

	/// Alocate array for magnitude and phase parts of frequency domain
	reXSize = magSize;
	imXSize = phaseSize;
	double* _realX = new double[reXSize];
	double* _imaginaryX = new double[imXSize];

	/// Polar To Rectangular Conversion
	for(int i = 0; i < magSize; i++) {

		_realX[i] = mag[i] * cos(phase[i]);
		_imaginaryX[i] = mag[i] * sin(phase[i]);
	}

	// Store the output signals
	if(*realX != NULL) delete[] *realX;
	*realX = _realX;
	if(*imaginaryX != NULL) delete[] *imaginaryX;
	*imaginaryX = _imaginaryX;
}

/**
* Unwrapping the phase. It is often easier to understand the phase if it does not have these discontinuities, even if it means that the phase extends above π, or below -π.
* a multiple of 2π is added or subtracted from each value of the phase.
* @author Petar Jerčić
* @date 09/12/2012
* @TODO N/A.
*/
void	DSPStatistic::PhaseUnwrapping(double* phase, int phaseSize, double** uwPhase, int& uwPhaseSize) {

	double PI = 4 * atan(1.0);
	
	/// Alocate array for magnitude and phase parts of frequency domain
	uwPhaseSize = phaseSize;
	double* _uwPhase = new double[uwPhaseSize];

	/// First point of all phase signals is 0
	_uwPhase[0] = 0;

	/// Go through the unwrappnig algorithm
	for(int i = 1; i < phaseSize; i++) {
		int c = (_uwPhase[i - 1] - phase[i]) / (2 * PI);
		_uwPhase[i] = phase[i] + c * 2 * PI;
	}

	// Store the output signal
	if(*uwPhase != NULL) delete[] *uwPhase;
	*uwPhase = _uwPhase;
}

/**
* This subroutine creates the complex frequency domain from the real frequency domain, generating negative frequencies between samples N/2 + 1 and N - 1.
* Upon entering to this method \a realFreqSize contains number of points in the signal, and \a reX and \a imX contain the real frequency domain in samples 0 to N/2.
* @author Petar Jerčić
* @date 16/12/2012
* @todo N/A.
*/
void	DSPStatistic::NegativeFrequencyGeneration(int realFreqSize, double* reX, double* imX) {

	/// Mirror the real part and negative mirror the imaginary against \a realFreqSize/2
	for(int i = realFreqSize / 2; i < realFreqSize; i++) {

		reX[i] = reX[realFreqSize - i - 1];
		imX[i] = - imX[realFreqSize - i - 1];
	}
}

/**
* This subroutine calculates the complex DFT by correlating the time domain signal with sine and cosine waves. Data are passed to this FFT subroutine in the arrays: \a reX and \a imX, each running from sample 0 to N.
* Upon entering to this method \a signalsSize contains number of points in the signal, and \a xR and \a xI contain the real and imaginary parts of time domain. All samples run from 0 to N/2.
* Upon return from the subroutine,\a reX and \a imX are overwritten with the frequency domain data.
* @author Petar Jerčić
* @date 17/12/2012
* @todo N/A.
*/
void	DSPStatistic::ComplexDFTByCorrelationTime(int signalsSize, double* xR, double* xI, double* reX, double* imX) {

	double PI = 4 * atan(1.0);

	/// Zero out the \a reX and \a imX so they can be used as accumulators during the correlation
	for(int i = 0; i < signalsSize; i++) {

		reX[i] = 0;
		imX[i] = 0;
	}

	/// Loop through each value in frequency domain
	for(int i = 0; i < signalsSize; i++)
		/// Correlate with the complex sinusoide \a sinReal & \a sinImag
		for(int j = 0; j < signalsSize; j++) {

			/// Calculate complex sinusoide
			double sinReal = cos(2 * PI * i * j / signalsSize);
			double sinImag = - sin(2 * PI * i * j / signalsSize);

			reX[i] += xR[j] * sinReal - xI[j] * sinImag; 
			imX[i] += xR[j] * sinImag - xI[j] * sinReal;
		}
}

/**
* This subroutine calculates the Fast Fourier Transform. The bit reversal sorting is done after the three nested loops. 
* Data are passed to this FFT subroutine in the arrays: \a reX and \a imX, each running from sample 0 to N.
* Upon entering to this method \a signalsSize contains number of points in the signal. All samples run from 0 to N.
* Upon return from the subroutine,\a reX and \a imX are overwritten with the frequency domain data.
* @author Petar Jerčić
* @date 17/12/2012
* @todo N/A.
*/
void	DSPStatistic::FastFourierTransform(int signalsSize, double* reX, double* imX) {

	double PI = 4 * atan(1.0);

	/// Set constants
	int nM1 = signalsSize - 1;
	int nD2 = signalsSize - 2;
	int m = log((float)signalsSize) / log(2.0);
	int k = nD2;
	double tR = 0;
	double tI = 0;

	/// Bit reversal sorting
	for(int i = 1; i < signalsSize - 1; i++) {

		if(i < k) {

			tR = reX[k];
			tI = imX[k];
			reX[k] = reX[i];
			imX[k] = imX[i];
			reX[i] = tR;
			imX[i] = tI;
		}

		int j = nD2;

		while(j <= k) {

			k -= j;
			j /= 2;
		}

		k += j; 
	}

	/// Loop for each stage
	for(int i = 0; i < m; i++) {

		int le = pow(2.0, i);
		int le2 = le / 2;

		double uR = 1;
		double uI = 0;

		/// Calculate sinusoide and cosinusoide values
		double sinReal = cos(PI / le2);
		double sinImag = - sin(PI / le2);

		/// Loop for each sub DFT
		for(int j = 0; j < le2; j++) {

			int jm1 = j - 1;

			/// Loop for each butterfly
			for(int l = 0; l < le2; l++) {

				int ip = l + le2;

				/// Butterfly calculation
				tR = reX[ip] * uR - imX[ip] * uI;
				tI = reX[ip] * uI - imX[ip] * uR;
				reX[ip] = reX[l] - tR;
				imX[ip] = imX[l] - tI;
				reX[l] += tR;
				imX[l] += tI;
			}

			tR = uR;
			uR = tR * sinReal - uI * sinImag;
			uI = tR * sinImag - uI * sinReal;
		}
	}
}

/**
* This subroutine calculates the Inverse Fast Fourier Transform calculating a Forward FFT, and then adjusting the data.
* Upon entering to this method \a signalsSize contains number of points in the signal for Inverse Fast Fourier Transform. All samples run from 0 to N.
* Data are passed to this FFT subroutine in the arrays: \a reX and \a imX, each holding real and imaginary part of the complex frequency domain.
* Upon return from the subroutine,\a reX and \a imX are overwritten with the complex time domain data.
* @author Petar Jerčić
* @date 17/12/2012
* @todo N/A.
*/
void	DSPStatistic::InverseFastFourierTransform(int signalsSize, double* reX, double* imX) {

	/// Change the sign of \a imX
	for(int i = 0; i < signalsSize; i++)
		imX[i] = - imX[i];

	/// Calculate forward FFT
	FastFourierTransform(signalsSize, reX, imX);

	/// Divido the time domain by \a signalsSize and change the sign of \a imX
	for(int i = 0; i < signalsSize; i++) {

		reX[i] /= signalsSize;
		imX[i] = - imX[i] / signalsSize;
	}
}

/**
* This subroutine calculates the Fast Fourier Transform for real numbers.
* Upon entering to this method \a signalsSize contains number of points in the \a signal for Fast Fourier Transform. All samples run from 0 to N.
* Upon return from the subroutine,\a reX and \a imX are overwritten with the real and imaginary data.
* @author Petar Jerčić
* @date 13/03/2013
* @todo Implement.
*/
void	DSPStatistic::FastFourierTransformReal(double* signal, int signalSize, double** reX, double** imX) {

	/// Alocate array for the real part of input signal
	double*	_reX = new double[signalSize];

	/// Alocate array for the imaginary part of input signal
	double*	_imX = new double[signalSize];

	/// Separate even and odd points
	int nh = signalSize / 2;
	for (int i = 0; i < nh; i++) {

		_reX[i] = signal[2 * i];
		_imX[i] = signal[2 * i + 1];
	}

	/// Calculate the \a signalSize / 2 point FFT (Table 12-3)
	FastFourierTransform(signalSize / 2, _reX, _imX);

	/// Even/Odd frequency domain decomposition
	int nm1 = signalSize;
	int nd2 = signalSize / 2;
	int n4 = signalSize / 2;
	for (int i = 1; i < n4; i++) {

		int im = nd2 - i;
		int ip2 = i + nd2;
		int ipm = im + nd2;
		_reX[ip2] = (_imX[i] + _imX[im]) / 2;
		_reX[ipm] = _reX[ip2];
		_imX[ip2] = -(_reX[i] - _reX[im]) / 2;
		_imX[ipm] = -_imX[ip2];
		_reX[i] = (_reX[i] + _reX[im]) / 2;
		_reX[im] = _reX[i];
		_imX[i] = (_imX[i] - _imX[im]) / 2;
		_imX[im] = -_imX[i];
	}

	_reX[signalSize * 3 / 4] = _imX[signalSize / 4];
	_reX[nd2] = _imX[0];
	_imX[signalSize * 3 / 4] = 0;
	_imX[nd2] = 0;
	_imX[signalSize / 4] = 0;
	_imX[0] = 0;

	/// Complete the last FFT stage
	double PI = 4 * atan(1.0);
	int l = log((double)signalSize) / log(2.0);
	int le = pow(2.0, l);
	int le2 = le / 2;
	double ur = 1;
	double ui = 0;
	double sr = cos(PI / le2);
	double si = -sin(PI / le2);
	double tr = 0;
	double ti = 0;
	for (int j = 1; j <= le2; j++) {

		int jm1 = j - 1;
		for (int i = jm1; i < nm1; i += le) {

			int ip = i + le2;
			tr = _reX[ip] * ur - _imX[ip] * ui;
			ti = _reX[ip] * ui + _imX[ip] * ur;
			_reX[ip] = _reX[i] - tr;
			_imX[ip] = _imX[i] - ti;
			_reX[i] = _reX[i] + tr;
			_imX[i] = _imX[i] + ti;
		}

		tr = ur;
		ur = tr * sr - ui * si;
		ui = tr * si + ui * sr;
	}

	// Store the output signal
	*reX = _reX;
	*imX = _imX;
}

/**
* This subroutine calculates the Inverse Fast Fourier Transform for real numbers.
* Upon entering to this method \a signalsSize contains number of points in the signal for Inverse Fast Fourier Transform. All samples run from 0 to N.
* Data are passed to this FFT subroutine in the arrays: \a reX and \a imX, each holding real and imaginary part of the frequency domain from 0 to n / 2.
* @author Petar Jerčić
* @date 16/01/2013
* @todo N/A
*/
double*	DSPStatistic::InverseFastFourierTransformReal(int signalsSize, double* reX, double* imX) { 

	/// Alocate array for real time domain
	double*	x = new double[signalsSize];

	/// Make the frequency domain symetrical
	for (int i = signalsSize / 2; i < signalsSize; i++) {

		reX[i] = reX[signalsSize - i];
		imX[i] = -imX[signalsSize - i];
	}

	/// Add real and imaginary parts together
	for (int i = 0; i < signalsSize; i++)
		x[i] = reX[i] + imX[i];

	/// Calculate the forward real DFT
	double* _reX;
	double* _imX;
	FastFourierTransformReal(x, signalsSize, &_reX, &_imX);

	/// Add the real and imaginary parts together and divide the time domain
	for (int i = 0; i < signalsSize; i++)
		x[i] = (_reX[i] + _imX[i]) / signalsSize;

	return x;
}

/**
* This method filters 5000 samples with a 101 point moving average filter, resulting in 4900 samples of filtered data.
* Upon entering to this method \a x contains the input signal and \a xSize its size.
* @author Petar Jerčić
* @date 10/01/2013
* @todo N/A.
*/
double*	DSPStatistic::MovingAverageFilter(double* x, int xSize) {

	/// Alocate array for the output signal
	double*	y = new double[xSize];

	/// Loop for each point in the output signal
	for (int i = 50; i < xSize - 50; i++) {

		/// Zero the output so it can be used as an accumulator
		y[i] = 0;

		/// Calculate the summation
		for (int j = -50; j <= 50; j++)
			y[i] += x[i + j];

		/// Complete the average by dividing
		y[i] /= 101;
	}

	return y;
}

/**
* Notice that this equation use two sources of data to calculate each point in the output: points from the input and previously calculated points from the output. 
* This is called a recursive equation, meaning that the result of one calculation is used in future calculations. 
* Upon entering to this method \a x contains the input signal and \a xSize its size.
* @author Petar Jerčić
* @date 10/01/2013
* @todo N/A.
*/
double*	DSPStatistic::MovingAverageFilterRecursion(double* x, int xSize) {

	/// Alocate array for the output signal
	double*	y = new double[xSize];

	/// Define a double precision accumulator
	double acc;

	acc = 0;
	/// Find y[50] by averaging points from x[0] to x[100]
	for (int i = 0; i <= 100; i++)
		acc += x[i];

	y[50] = acc / 101;

	/// Recursive moving average filter
	for (int i = 51; i < xSize - 50; i++) {

		acc += x[i + 50] - x[i - 51];
		y[i] = acc / 101;
	}

	return y;
}

/**
* This program will filter the iput data with a 101 point window-sinc filter. In this example, we will assume that the EEG signal has been amplified by analog electronics, and then digitized at a sampling rate of 100 samples per second. 
* Our goal is to separate the alpha from the beta rhythms. To do this, we will design a digital low-pass filter with a cutoff frequency of 14 hertz, or 0.14 of the sampling rate. 
* The transition bandwidth will be set at 4 hertz, or 0.04 of the sampling rate. From Eq. 16-3, the filter kernel needs to be about 101 points long, and we will arbitrarily choose to use a Hamming window.
* For the filter to have unity gain at DC, the constant K must be chosen such that the sum of all the samples is equal to one. 
* In practice, ignore K during the calculation of the filter kernel, and then normalize all of the samples as needed. Also notice how the calculation is handled at the center of the sinc, i = M/2, which involves a division by zero.
* Upon entering to this method \a x contains the input signal and \a xSize its size.
* @author Petar Jerčić
* @date 13/01/2013
* @todo N/A.
*/
double*	DSPStatistic::LowPassWindowdSincFilter(double* x, int xSize) {

	double PI = 4 * atan(1.0);

	/// Alocate array for the output signal
	double*	y = new double[xSize];

	/// Set the cutoff frequency (between 0 and 0.5)
	double fc = 0.14;

	/// Set the filter length (+ 1 = 101 point)
	int M = 100;

	/// Alocate array for filter kernel
	double*	h = new double[M + 1];

	/// Calculate the filter kernel  via Eq. 16-4
	for (int i = 0; i < M + 1; i++) {

		if (i - M / 2 == 0)		h[i] = 2 * PI * fc;
		else					h[i] = sin(2 * PI * fc * (i - M / 2)) / (i - M / 2);

		h[i] *= (0.54 - 0.46 * cos(2 * PI * i / M));
	}

	/// Normalize the low pass filter kernel for unity gain at DC
	double sum = 0;
	for (int i = 0; i < M + 1; i++)
		sum += h[i];

	for (int i = 0; i < M + 1; i++)
		h[i] /= sum;

	/// Convolve the input signal and the filter kernel
	for (int i = M; i < xSize; i++) {

		y[i] = 0;

		for (int j = 0; j < M + 1; j++)
			y[i] += x[i - j] * h[j];
	}

	return y;
}

/**
* This program will calculate a 801 point band-pass filter kernel. In a second example, we will design a band-pass filter to isolate a signaling tone in an audio signal, such as when a button on a telephone is pressed. 
* We will assume that the signal has been digitized at 10 kHz, and the goal is to isolate an 80 hertz band of frequencies centered on 2 kHz. In terms of the sampling rate, we want to block all frequencies below 0.196 and above 0.204 (corresponding to 1960 hertz and 2040 hertz, respectively). 
* To achieve a transition bandwidth of 50 hertz (0.005 of the sampling rate), we will make the filter kernel 801 points long, and use a Blackman window. The design involves several steps. 
* First, two low-pass filters are designed, one with a cutoff at 0.196, and the other with a cutoff at 0.204. This second filter is then spectrally inverted, making it a high-pass filter (see Chapter 14, Fig. 14-6). 
* Next, the two filter kernels are added, resulting in a band-reject filter (see Fig. 14-8). Finally, another spectral inversion makes this into the desired band-pass filter.
* Upon entering to this method \a fiterKernel Will hold the band pass filter kernel and \a fiterKernelSize Will hold the band pass filter kernel size.
* @author Petar Jerčić
* @date 16/01/2013
* @todo N/A.
*/
void	DSPStatistic::BandPassWindowdSincFilter(double** fiterKernel, int& fiterKernelSize) {

	double PI = 4 * atan(1.0);

	/// Set the filter length (+ 1 = 801 point)
	int M = 800;

	/// Alocate array for the lower cutoff
	double*	a = new double[M + 1];

	/// Alocate array for the upper cutoff
	double*	b = new double[M + 1];

	/// Alocate array for filter kernel
	double*	h = new double[M + 1];

	/// Calculte the first low-pass filter kernel via Eq. 16-4 with a cutoff frequency of 0.196, store in \a a
	double fc = 0.196;
	for (int i = 0; i < M + 1; i++) {

		if (i - M / 2 == 0)		a[i] = 2 * PI * fc;
		else					a[i] = sin(2 * PI * fc * (i - M / 2)) / (i - M / 2);

		a[i] *= (0.42 - 0.5 * cos(2 * PI * i / M) + 0.08 * cos(4 * PI * i / M));
	}

	/// Normalize the first low pass filter kernel for unity gain at DC
	double sum = 0;
	for (int i = 0; i < M + 1; i++)
		sum += a[i];

	for (int i = 0; i < M + 1; i++)
		a[i] /= sum;

	/// Calculte the second low-pass filter kernel via Eq. 16-4 with a cutoff frequency of 0.204, store in \a b
	fc = 0.204;
	for (int i = 0; i < M + 1; i++) {

		if (i - M / 2 == 0)		b[i] = 2 * PI * fc;
		else					b[i] = sin(2 * PI * fc * (i - M / 2)) / (i - M / 2);

		b[i] *= (0.42 - 0.5 * cos(2 * PI * i / M) + 0.08 * cos(4 * PI * i / M));
	}

	/// Normalize the second low pass filter kernel for unity gain at DC
	sum = 0;
	for (int i = 0; i < M + 1; i++)
		sum += b[i];

	for (int i = 0; i < M + 1; i++)
		b[i] /= sum;

	/// Change the low-pass filter in \a b into a high-pass filter kernel using spectral inversion
	for (int i = 0; i < M + 1; i++)
		b[i] = -b[i];
	b[400] += 1;

	/// Add the low-pass filter kernel in \a a to a high-pass filter kernel in \a b, to form a band-reject filter kernel stored in \a h
	for (int i = 0; i < M + 1; i++)
		h[i] = a[i] + b[i];

	/// Change the band-reject filter kernel into a band-pass filter kernel by using spectral inversion
	for (int i = 0; i < M + 1; i++)
		h[i] = -h[i];
	h[400] += 1;

	*fiterKernel	= h;
	fiterKernelSize = M + 1;
}

/**
* This program will calculate a arbitrary size point band-pass filter kernel. In a second example, we will design a band-pass filter to isolate a signaling tone in an audio signal, such as when a button on a telephone is pressed. 
* We will assume that the signal has been digitized at 10 kHz, and the goal is to isolate an 80 hertz band of frequencies centered on 2 kHz. In terms of the sampling rate, we want to block all frequencies below 0.196 and above 0.204 (corresponding to 1960 hertz and 2040 hertz, respectively). 
* To achieve a transition bandwidth of 50 hertz (0.005 of the sampling rate), we will make the filter kernel arbitrarily points long, and use a Blackman window. The design involves several steps. 
* First, two low-pass filters are designed, one with a cutoff at 0.196, and the other with a cutoff at 0.204. This second filter is then spectrally inverted, making it a high-pass filter (see Chapter 14, Fig. 14-6). 
* Next, the two filter kernels are added, resulting in a band-reject filter (see Fig. 14-8). Finally, another spectral inversion makes this into the desired band-pass filter.
* Upon entering to this method \a fiterKernelSize holds the band pass filter kernel size.
* @author Petar Jerčić
* @date 16/01/2013
* @todo N/A.
*/
double*	DSPStatistic::BandPassWindowdSincFilterCustom(int fiterKernelSize) {

	double PI = 4 * atan(1.0);

	/// Set the filter length (fiterKernelSize - 1 points)
	int M;
	if (fiterKernelSize % 2 != 0)	M = fiterKernelSize - 1;
	else							M = fiterKernelSize;

	/// Alocate array for the lower cutoff
	double*	a = new double[M + 1];

	/// Alocate array for the upper cutoff
	double*	b = new double[M + 1];

	/// Alocate array for filter kernel
	double*	h = new double[M + 1];

	/// Calculte the first low-pass filter kernel via Eq. 16-4 with a cutoff frequency of 0.196, store in \a a
	double fc = 0.196;
	for (int i = 0; i < M + 1; i++) {

		if (i - M / 2 == 0)		a[i] = 2 * PI * fc;
		else					a[i] = sin(2 * PI * fc * (i - M / 2)) / (i - M / 2);

		a[i] *= (0.42 - 0.5 * cos(2 * PI * i / M) + 0.08 * cos(4 * PI * i / M));
	}

	/// Normalize the first low pass filter kernel for unity gain at DC
	double sum = 0;
	for (int i = 0; i < M + 1; i++)
		sum += a[i];

	for (int i = 0; i < M + 1; i++)
		a[i] /= sum;

	/// Calculte the second low-pass filter kernel via Eq. 16-4 with a cutoff frequency of 0.204, store in \a b
	fc = 0.204;
	for (int i = 0; i < M + 1; i++) {

		if (i - M / 2 == 0)		b[i] = 2 * PI * fc;
		else					b[i] = sin(2 * PI * fc * (i - M / 2)) / (i - M / 2);

		b[i] *= (0.42 - 0.5 * cos(2 * PI * i / M) + 0.08 * cos(4 * PI * i / M));
	}

	/// Normalize the second low pass filter kernel for unity gain at DC
	sum = 0;
	for (int i = 0; i < M + 1; i++)
		sum += b[i];

	for (int i = 0; i < M + 1; i++)
		b[i] /= sum;

	/// Change the low-pass filter in \a b into a high-pass filter kernel using spectral inversion
	for (int i = 0; i < M + 1; i++)
		b[i] = -b[i];
	b[M / 2] += 1;

	/// Add the low-pass filter kernel in \a a to a high-pass filter kernel in \a b, to form a band-reject filter kernel stored in \a h
	for (int i = 0; i < M + 1; i++)
		h[i] = a[i] + b[i];

	/// Change the band-reject filter kernel into a band-pass filter kernel by using spectral inversion
	for (int i = 0; i < M + 1; i++)
		h[i] = -h[i];
	h[M / 2] += 1;

	return h;
}

/**
* This program converts an aliased 1024 points impulse response into an M + 1 point filter kernel. 
* The frequency response we want the filter to produce is shown in (a). To say the least, it is very irregular and would be virtually impossible to obtain with analog electronics. This ideal frequency response is defined by an array of numbers that have been selected, not some mathematical equation. In this example, there are 513 samples spread between 0 and 0.5 of the sampling rate. More points could be used to better represent the desired frequency response, while a smaller number may be needed to reduce the computation time during the filter design. However, these concerns are usually small, and 513 is a good length for most applications.
* The next step is to take the Inverse DFT to move the filter into the time domain. The quickest way to do this is to convert the frequency domain to rectangular form, and then use the Inverse FFT. This results in a 1024 sample signal running from 0 to 1023, as shown in (b). This is the impulse response that corresponds to the frequency response we want; however, it is not suitable for use as a filter kernel (more about this shortly). Just as in the last chapter, it needs to be shifted, truncated, and windowed. In this example, we will design the filter kernel with M = 40, i.e., 41 points running from sample 0 to sample 40.
* Upon entering to this method \a reX contains the signal to convert and \a reXSize its size.
* @author Petar Jerčić
* @date 13/01/2013
* @todo Besides the desired magnitude array shown in (a), there must be a corresponding phase array of the same length. In this example, the phase of the desired frequency response is entirely zero (this array is not shown in Fig. 17-1). Just as with the magnitude array, the phase array can be loaded with any arbitrary curve you would like the filter to produce. However, remember that the first and last samples (i.e., 0 and 512) of the phase array must have a value of zero (or a multiple of 2π, which is the same thing). The frequency response can also be specified in rectangular form by defining the array entries for the real and imaginary parts, instead of using the magnitude and phase.
*/
double*	DSPStatistic::CustomFilterDesign(double* reX, int reXSize) {

	double PI = 4 * atan(1.0);

	/// Set the filter kernel length (+ 1 = 41 point)
	int M = 40;

	/// Alocate the temporary storage buffer
	double* temp = new double[reXSize];

	/// Shift (rotate) the signal M / 2 points to the right
	for (int i = 0; i < reXSize; i++) {

		int index = i + M / 2;

		if (index >= reXSize) index -= reXSize;

		temp[index] = reX[i];
	}

	for (int i = 0; i < reXSize; i++)
		reX[i] = temp[i];

	/// Truncate the window the signal
	for (int i = 0; i < reXSize; i++) {
		if (i <= M) reX[i] *= (0.54 - 0.46 * cos(2 * PI * i / M));
		if (i > M)	reX[i] = 0;
	}

	/// The filter kernel now resides in reX[0] to reX[40]
	return reX;
}

/**
* his program filters a 10 million point signal by convolving it with a 400 point filter kernel. This is done by breaking the input signal into 16000 segments, with each segment having 625 points. When each of these segments is convolved with the filter kernel, an output segment of 625 + 400 - 1 = 1024 points is produced. Thus, 1024 point FFTs are used. After defining and initializing all the arrays (lines 130 to 230), the first step is to calculate and store the frequency response of the filter (lines 250 to 310). Line 260 calls a mythical subroutine that loads the filter kernel into XX[0] through XX[399], and sets XX[400] through XX[1023] to a value of zero. The subroutine in line 270 is the FFT, transforming the 1024 samples held in XX[ ] into the 513 samples held in REX[ ] & IMX[ ], the real and imaginary parts of the frequency response. These values are transferred into the arrays REFR[ ] & IMFR[ ] (for: REal and IMaginary Frequency Response), to be used later in the program.
* The FOR-NEXT loop between lines 340 and 580 controls how the 16000 segments are processed. In line 360, a subroutine loads the next segment to be processed into XX[0] through XX[624], and sets XX[625] through XX[1023] to a value of zero. In line 370, the FFT subroutine is used to find this segment's frequency spectrum, with the real part being placed in the 513 points of REX[ ], and the imaginary part being placed in the 513 points of IMX[ ]. Lines 390 to 430 show the multiplication of the segment's frequency spectrum, held in REX[ ] & IMX[ ], by the filter's frequency response, held in REFR[ ] and IMFR[ ]. The result of the multiplication is stored in REX[ ] & IMX[ ], overwriting the data previously there. Since this is now the frequency spectrum of the output segment, the IFFT can be used to find the output segment. This is done by the mythical IFFT subroutine in line 450, which transforms the 513 points held in REX[ ] & IMX[ ] into the 1024 points held in XX[ ], the output segment.
* Lines 470 to 550 handle the overlapping of the segments. Each output segment is divided into two sections. The first 625 points (0 to 624) need to be combined with the overlap from the previous output segment, and then written to the output signal. The last 399 points (625 to 1023) need to be saved so that they can overlap with the next output segment.
* To understand this, look back at Fig 18-1. Samples 100 to 199 in (g) need to be combined with the overlap from the previous output segment, (f), and can then be moved to the output signal (i). In comparison, samples 200 to 299 in (g) need to be saved so that they can be combined with the next output segment, (h).
* Now back to the program. The array OLAP[ ] is used to hold the 399 samples that overlap from one segment to the next. In lines 470 to 490 the 399 values in this array (from the previous output segment) are added to the output segment currently being worked on, held in XX[ ]. The mythical subroutine in line 550 then outputs the 625 samples in XX[0] to XX[624] to the file holding the output signal. The 399 samples of the current output segment that need to be held over to the next output segment are then stored in OLAP[ ] in lines 510 to 530.
* After all 0 to 15999 segments have been processed, the array, OLAP[ ], will contain the 399 samples from segment 15999 that should overlap segment 16000. Since segment 16000 doesn't exist (or can be viewed as containing all zeros), the 399 samples are written to the output signal in line 600. This makes the length of the output signal points. This matches the length of input signal, plus the length of the filter kernel, minus 1.
* Upon entering to this method \a x contains the input signal and \a xSize its size.
* @author Petar Jerčić
* @date 15/01/2013
* @todo N/A
*/
double*	DSPStatistic::FFTConvolution(double* x, int xSize) {

	/// Initialize the arrays

	/// Initialize the output
	double* xConvOut = new double[xSize];

	/// Initialize the time domain signal array (for the FFT)
	double* xX = new double[1024];

	/// Initialize the real part of the frequency domain (for the FFT)
	double* reX;

	/// Initialize the imaginary part of the frequency domain (for the FFT)
	double* imX;

	/// Initialize the real part of the filtes's frequency response
	double* reFr = new double[512];

	/// Initialize the imaginary part of the filtes's frequency response
	double* imFr = new double[512];

	/// Initialize the overlaping samples from segment to segment
	double* olap = new double[400];

	/// Zero the array holding the overlaping samples
	for (int i = 0; i < 400; i++)
		olap[i] = 0;

	/// Find and store the filters frequency response

	/// Load the filter kernel into \a xX
	double* tempFiltrKer = BandPassWindowdSincFilterCustom(401);
	for (int i = 0; i < 1024; i++) {
		
		if (i < 400)	xX[i] = tempFiltrKer[i];
		else			xX[i] = 0;
	}

	/// FFT the loaded filter kernel into \a xX
	FastFourierTransformReal(xX, 1024, &reX, &imX);

	/// Save the frequency response to \a reFr and \a imFr
	for (int i = 0; i < 512; i++) {
		
		reFr[i] = reX[i];
		imFr[i] = imX[i];
	}

	/// Process each of 16000 segments
	for (int i = 0; i < 16000; i++) {

		/// Load the next input segment into \a xX
		xX = &x[i * (xSize / 16000)];

		/// FFT the loaded input segment into \a xX
		FastFourierTransformReal(xX, xSize / 16000, &reX, &imX);

		/// Multiply the frequency spectrum by the frequency response
		for (int j = 0; j < 512; j++) {

			double temp = reX[j] * reFr[j] - imX[j] * imFr[j];
			imX[j] = reX[j] * imFr[j] + imX[j] * reFr[j];
			reX[j] = temp;
		}

		/// IFFT on the real and imaginary part of the frequency domain
		xX = InverseFastFourierTransformReal(xSize, reX, imX);	///@TODO Chech xSize here

		/// Add the last's segment overlap to this segment
		for (int j = 0; j < 400; j++)
			xX[j] += olap[j];

		/// Save the samples that will overlap next segment
		for (int j = 625; j < 1024; j++)
			olap[j - 625] = xX[j];

		/// Store the signal
		for (int j = i * 625; j < (i + 1) * 625; j++)
			xConvOut[j] = olap[j - i * 625];	///@TODO xConvOut size
	}

	return xConvOut;
}

/**
* The "a" and "b" values that define the filter are called the recursion coefficients. In actual practice, no more than about a dozen recursion coefficients can be used or the filter becomes unstable 
* (i.e., the output continually increases or oscillates).
* Upon entering to this method \a x contains the input signal and \a xSize its size, \a a contains the recursion coefficients and \a aSize its array size, 
* \a b contains the recursion coefficients, b[0] is ignored, and \a bSize its array size.
* @author Petar Jerčić
* @date 17/01/2013
* @todo zero out the phase by bidirectional filtering.
*/
double*	DSPStatistic::RecursiveFilter(double* x, int xSize, double* a, int aSize, double* b, int bSize) {

	/// Allocate the array for the output signal
	double* y = new double[xSize];

	for (int i = 0; i < max(aSize, bSize); i++)
		y[i] = 0;

	for (int i = max(aSize, bSize); i < xSize; i++) {
		y[i] = 0;
		for (int j = 0; j < aSize; j++)
			y[i] += a[j] * x[i - j];
		for (int j = 1; j < bSize; j++)
			y[i] += b[j] * y[i - j];
	}

	return y;
}

/**
* The "a" and "b" values that define the filter are called the recursion coefficients. In actual practice, no more than about a dozen recursion coefficients can be used or the filter becomes unstable 
* (i.e., the output continually increases or oscillates).
* Upon entering to this method \a x contains the input signal and \a xSize its size, \a a contains the recursion coefficients and \a aSize its array size, 
* \a b contains the recursion coefficients, b[0] is ignored, and \a bSize its array size. Standard a0 = 0.15, b0 = 0.85
* @author Petar Jerčić
* @date 17/01/2013
* @todo N/A
*/
double*	DSPStatistic::SinglePoleLowPassRecursiveFilter(double* signal, int signalSize, double cutoffFreq) {

	double PI = 4 * atan(1.0);

	/// Alocate the array for the \a a recursion coefficients
	double* a = new double[1];

	/// Alocate the array for the \a b recursion coefficients
	double* b = new double[2];

	/// Calculate the \a x, \a a, \a b coefficients
	double X = exp(-2 * PI * cutoffFreq);
	a[0] = 1 - X;
	b[1] = X;

	return RecursiveFilter(signal, signalSize, a, 1, b, 2);
}

/**
* The "a" and "b" values that define the filter are called the recursion coefficients. In actual practice, no more than about a dozen recursion coefficients can be used or the filter becomes unstable 
* (i.e., the output continually increases or oscillates).
* Upon entering to this method \a x contains the input signal and \a xSize its size, \a a contains the recursion coefficients and \a aSize its array size, 
* \a b contains the recursion coefficients, b[0] is ignored, and \a bSize its array size. Standard a0 = 0.93, a1 = -0.93, b0 = 0.86.
* @author Petar Jerčić
* @date 17/01/2013
* @todo N/A
*/
double*	DSPStatistic::SinglePoleHighPassRecursiveFilter(double* signal, int signalSize, double cutoffFreq) {

	double PI = 4 * atan(1.0);

	/// Alocate the array for the \a a recursion coefficients
	double* a = new double[2];

	/// Alocate the array for the \a b recursion coefficients
	double* b = new double[2];

	/// Calculate the \a x, \a a, \a b coefficients
	double X = exp(-2 * PI * cutoffFreq);
	a[0] = (1 + X) / 2;
	a[1] = -(1 + X) / 2;
	b[1] = X;

	return RecursiveFilter(signal, signalSize, a, 2, b, 2);
}

/**
* The "a" and "b" values that define the filter are called the recursion coefficients. In actual practice, no more than about a dozen recursion coefficients can be used or the filter becomes unstable 
* (i.e., the output continually increases or oscillates).
* Upon entering to this method \a x contains the input signal and \a xSize its size, \a a contains the recursion coefficients and \a aSize its array size, 
* \a b contains the recursion coefficients, b[0] is ignored, and \a bSize its array size. Standard a0 = 0.93, a1 = -0.93, b0 = 0.86.
* @author Petar Jerčić
* @date 17/01/2013
* @todo N/A
*/
double*	DSPStatistic::FourStageLowPassRecursiveFilter(double* signal, int signalSize, double cutoffFreq) {

	double PI = 4 * atan(1.0);

	/// Alocate the array for the \a a recursion coefficients
	double* a = new double[1];

	/// Alocate the array for the \a b recursion coefficients
	double* b = new double[5];

	/// Calculate the \a x, \a a, \a b coefficients, but with changed eq.
	double X = exp(-14.445 * cutoffFreq);
	a[0] = pow(1 - X, 4);
	b[1] = 4 * X;
	b[2] = -6 * pow(X, 2);
	b[3] = 4 * pow(X, 3);
	b[4] = -pow(X, 4);

	return RecursiveFilter(signal, signalSize, a, 1, b, 5);
}

/**
* The "a" and "b" values that define the filter are called the recursion coefficients. In actual practice, no more than about a dozen recursion coefficients can be used or the filter becomes unstable 
* (i.e., the output continually increases or oscillates).
* The narrowest bandwidth that can be obtain with single precision is about 0.0003 of the sampling frequency. When pushed beyond this limit, the attenuation of the notch will degrade.
* Upon entering to this method \a x contains the input signal and \a xSize its size, \a a contains the recursion coefficients and \a aSize its array size, 
* \a b contains the recursion coefficients, b[0] is ignored, and \a bSize its array size. Standard a0 = 0.93, a1 = -0.93, b0 = 0.86.
* @author Petar Jerčić
* @date 17/01/2013
* @todo N/A
*/
double*	DSPStatistic::BandPassRecursiveFilter(double* signal, int signalSize, double centerFreq, double bandwidth) {

	double PI = 4 * atan(1.0);

	/// Alocate the array for the \a a recursion coefficients
	double* a = new double[3];

	/// Alocate the array for the \a b recursion coefficients
	double* b = new double[3];

	/// Calculate the \a R, \a K, \a a, \a b coefficients.
	double R = 1 - 3 * bandwidth;
	double K = (1 - 2 * cos(2 * PI * centerFreq) + pow(R, 2)) / (2 - 2 * cos(2 * PI * centerFreq));
	a[0] = 1 - K;
	a[1] = 2 * (K - R) * cos(2 * PI * centerFreq);
	a[2] = pow(R, 2) - K;
	b[1] = 2 * R * cos(2 * PI * centerFreq);
	a[2] = -pow(R, 2);

	return RecursiveFilter(signal, signalSize, a, 3, b, 3);
}

/**
* Two parameters must be selected before using these equations: f, the center frequency, and BW, the bandwidth (measured at an amplitude of 0.707). Both of these are expressed as a fraction of the sampling frequency, and therefore must be between 0 and 0.5. From these two specified values, calculate the intermediate variables: R and K, and then the recursion coefficients. 
* The narrowest bandwidth that can be obtain with single precision is about 0.0003 of the sampling frequency. When pushed beyond this limit, the attenuation of the notch will degrade.
* Upon entering to this method \a x contains the input signal and \a xSize its size, \a a contains the recursion coefficients and \a aSize its array size, 
* \a b contains the recursion coefficients, b[0] is ignored, and \a bSize its array size. Standard a0 = 0.93, a1 = -0.93, b0 = 0.86.
* @author Petar Jerčić
* @date 17/01/2013
* @todo N/A
*/
double*	DSPStatistic::BandRejectRecursiveFilter(double* signal, int signalSize, double centerFreq, double bandwidth) {

	double PI = 4 * atan(1.0);

	/// Alocate the array for the \a a recursion coefficients
	double* a = new double[3];

	/// Alocate the array for the \a b recursion coefficients
	double* b = new double[3];

	/// Calculate the \a R, \a K, \a a, \a b coefficients.
	double R = 1 - 3 * bandwidth;
	double K = (1 - 2 * cos(2 * PI * centerFreq) + pow(R, 2)) / (2 - 2 * cos(2 * PI * centerFreq));
	a[0] = K;
	a[1] = -2 * K * cos(2 * PI * centerFreq);
	a[2] = K;
	b[1] = 2 * R * cos(2 * PI * centerFreq);
	a[2] = -pow(R, 2);

	return RecursiveFilter(signal, signalSize, a, 3, b, 3);
}

/**
* Two parameters must be selected before using these equations: f, the center frequency, and BW, the bandwidth (measured at an amplitude of 0.707). Both of these are expressed as a fraction of the sampling frequency, and therefore must be between 0 and 0.5. From these two specified values, calculate the intermediate variables: R and K, and then the recursion coefficients. 
* The narrowest bandwidth that can be obtain with single precision is about 0.0003 of the sampling frequency. When pushed beyond this limit, the attenuation of the notch will degrade.
* Upon entering to this method \a x contains the input signal and \a xSize its size, \a a contains the recursion coefficients and \a aSize its array size, 
* \a b contains the recursion coefficients, b[0] is ignored, and \a bSize its array size. Standard a0 = 0.93, a1 = -0.93, b0 = 0.86.
* @author Petar Jerčić
* @date 17/01/2013
* @todo N/A
*/
void	DSPStatistic::ChebyshevFilterRecursionCoefficients(double cutoffFreq, bool lowPass, int percRipple, int noPoles, double** _a, double** _b) {

	double PI = 4 * atan(1.0);

	/// Alocate the array for the \a a recursion coefficients
	double* a = new double[22];

	/// Alocate the array for the \a b recursion coefficients
	double* b = new double[22];

	/// Alocate temporary array for the \a a recursion coefficients combining stages
	double* tempA = new double[22];

	/// Alocate temporary array for the \a b recursion coefficients combining stages
	double* tempB = new double[22];

	/// Zero out the arrays for recursion coefficients
	for (int i = 0; i < 22; i++) {

		a[i] = 0;
		b[i] = 0;
	}

	a[2] = 1;
	b[2] = 1;

	/// Loop for each pole pair
	for (int i = 0; i < noPoles / 2; i++) {

		double* calcA;
		double* calcB;
		ChebyshevFilterParametars(cutoffFreq, lowPass, percRipple, noPoles, i, &calcA, &calcB);

		/// Add coefficients to the cascade
		for (int j = 0; j < 22; j++) {

			tempA[j] = a[j];
			tempB[j] = b[j];
		}

		for (int j = 2; j < 22; j++) {

			a[j] = calcA[0] * tempA[j] + calcA[1] * tempA[j - 1] + calcA[2] * tempA[j - 2];
			b[j] =tempB[j] + calcB[1] * tempB[j - 1] + calcB[2] * tempB[j - 2];
		}
	}

	/// Finish combining coefficients
	b[2] = 0;
	for (int i = 0; i < 20; i++) {

		a[i] = a[i + 2];
		b[i] = - b[i + 2];
	}

	/// Normalize the gain
	double sa = 0;
	double sb = 0;

	for (int i = 0; i < 20; i++) {
		if (lowPass) {
			sa += a[i];
			sb += b[i];
		}
		else {
			sa += a[i] * pow(-1.0, i);
			sb += b[i] * pow(-1.0, i);
		}
	}

	double gain = sa / (1 - sb);

	for (int i = 0; i < 20; i++)
		a[i] /= gain;

	*_a = a;
	*_b = b;
}

/**
* Two parameters must be selected before using these equations: f, the center frequency, and BW, the bandwidth (measured at an amplitude of 0.707). Both of these are expressed as a fraction of the sampling frequency, and therefore must be between 0 and 0.5. From these two specified values, calculate the intermediate variables: R and K, and then the recursion coefficients. 
* The narrowest bandwidth that can be obtain with single precision is about 0.0003 of the sampling frequency. When pushed beyond this limit, the attenuation of the notch will degrade.
* Upon entering to this method \a x contains the input signal and \a xSize its size, \a a contains the recursion coefficients and \a aSize its array size, 
* \a b contains the recursion coefficients, b[0] is ignored, and \a bSize its array size. Standard a0 = 0.93, a1 = -0.93, b0 = 0.86.
* @author Petar Jerčić
* @date 17/01/2013
* @todo N/A
*/
void	DSPStatistic::ChebyshevFilterParametars(double cutoffFreq, bool lowPass, int percRipple, int noPoles, int count, double ** _a, double **_b) {

	double PI = 4 * atan(1.0);

	/// Alocate the array for the \a a recursion coefficients
	double* a = new double[22];

	/// Alocate the array for the \a b recursion coefficients
	double* b = new double[22];

	/// Calculate the pole location on the unit circle
	double rp = -cos(PI / (noPoles * 2) + (count - 1) * PI / noPoles);
	double ip = sin(PI / (noPoles * 2) + (count - 1) * PI / noPoles);


	/// Warp from the circle to ellipse
	if (percRipple > 0) {

		double es = sqrt(pow(100.0 / (100 - percRipple), 2) - 1);
		double vx = (1 / noPoles) * log((1 / es) + sqrt((1 / pow(es, 2)) + 1));
		double kx = (1 / noPoles) * log((1 / es) + sqrt((1 / pow(es, 2)) - 1));
		kx = (exp(kx) + exp(-kx)) / 2;
		rp *= ((exp(vx) - exp(-vx)) / 2) / kx;
		ip *= ((exp(vx) + exp(-vx)) / 2) / kx;
	}

	/// Alocate the temporory arrays
	double* x = new double[3];
	double* y = new double[3];

	/// s-damain to z-domain conversion
	double t = 2 * tan(1.0 / 2);
	double w = 2 * PI * cutoffFreq;
	double m = pow(rp, 2) + pow(ip, 2);
	double d = 4 - 4 * rp * t + m * pow(t, 2);
	x[0] = pow(t, 2) / d;
	x[1] = 2 * pow(t, 2) / d;
	x[2] = pow(t, 2) / d;
	y[1] = (8 - 2 * m * pow(t, 2)) / d;
	y[2] = (-4 -4 * rp * t - m * pow(t, 2)) / d;

	/// Low pass to low pass, or low pass to high pass
	double k;
	if (lowPass)	k = sin(1 / 2 - w / 2) / sin(1 / 2 + w / 2);
	else			k = -cos(w / 2 + 1 / 2) / cos(w / 2 - 1 / 2);

	d = 1 + y[1] * k - y[2] * pow(k, 2);

	a[0] = (x[0] - x[1] * k + x[2] * pow(k, 2)) / d;
	a[1] = (-2 * x[0] * k + x[1] + x[1] * pow(k, 2) - 2 * x[2] * k) /d;
	a[2] = (x[0] * pow(k, 2) - x[1] * k + x[2]) /d;
	b[1] = (2 * k + y[1] + y[1] * pow(k, 2) - 2 * y[2] * k) / d;
	b[2] = (-pow(k, 2) - y[1] * k + y[2]) / d;

	if (!lowPass) {
		a[1] = -a[1];
		b[1] = -b[1];
	}

	*_a = a;
	*_b = b;
}

/**
* The neural network used in this example is the traditional three-layer, fully interconnected architecture, as shown in Figs. 26-5 and 26-6. There are 101 nodes in the input layer (100 pixel values plus a bias node), 10 nodes in the hidden layer, and 1 node in the output layer. When a 100 pixel image is applied to the input of the network, we want the output value to be close to one if a vowel is present, and near zero if a vowel is not present. Don't be worried that the input signal was acquired as a two-dimensional array (10×10), while the input to the neural network is a one-dimensional array. This is your understanding of how the pixel values are interrelated; the neural network will find relationships of its own. 
* @author Petar Jerčić
* @date 22/01/2013
* @todo check \a esum for conversion
*/
void	DSPStatistic::NeuralNetworkTraining() {

	/// Initialize
	double MU = 0.000005;

	/// Alocate the array for the input layer signal + bias term
	/// The array elements: X1[0] through X1[99], hold the input layer values. In addition, X1[100] always holds a value of 1, providing the input to the bias node.
	double* x1 = new double[101];

	/// Alocate the array for the hidden layer signal
	double* x2 = new double[10];

	/// Alocate the array for the hidden layer weights
	double** wh = new double*[10];
	for(int i = 0; i < 10; i++)
		wh[i] = new double[101];

	/// Alocate the array for the output layer weights
	double* wo = new double[10];

	/// initialize random seed
	srand (time(NULL));

	/// Set weights to random values
	for (int i = 0; i < 10; i++ ) {

		/// Output layer weights, from -0.5 to 0.5
		wo[i] = rand() / RAND_MAX - 0.5;

		/// Hidden layer weights
		for (int j = 0; j < 101; j++)
			wh[i][j] = (rand() / RAND_MAX - 0.5) / 1000;
	}

	/// Iteration loop
	/// Loop for 800 iterations
	for (int iter = 0; iter < 800; iter++) {

		/// Clear the error accumulator \a esum
		double esum = 0;

		/// Loop for each data point in the training set
		for (int dataPointCount = 0; dataPointCount < 260; dataPointCount++) {

			/// Load \a x1 with the training set
			bool target = LoadNeuralNetworkTrainingSet(x1, 101, dataPointCount % 2);

			/// Find the error for this data point \a edata
			double x3;
			double edata = CalculateNeuralNetworkError(x1, 101, x2, 10, &x3, wh, wo, target);

			/// Accumulate error for this iteration
			esum += pow(edata, 2);

			/// Find the new weights
			FindNeuralNetworkWeights(x1, 101, x2, 10, x3, wh, wo, edata, MU);
		}

		// check \a esum for conversion
	}

	/// Save the weights
	// sub
}

/**
* Implementation
* @author Petar Jerčić
* @date 22/01/2013
* @todo N/A
*/
bool	DSPStatistic::LoadNeuralNetworkTrainingSet(double* x1, int x1Size, int dataSet) {

	srand (time(NULL));

	/// Fill the training data data
	for (int i = 0; i < x1Size -1; i++) {

		if(dataSet == 0)
			x1[i] = 0;
		else
			x1[i] = rand() / RAND_MAX;
	}

	/// Set the bias value
	x1[x1Size -1] = 1;

	if(dataSet == 0)
		return false;
	else
		return true;
}

/**
* Implementation
* @author Petar Jerčić
* @date 22/01/2013
* @todo N/A
*/
double	DSPStatistic::CalculateNeuralNetworkError(double* x1, int x1Size, double* x2, int x2Size, double* x3, double** wh, double* wo, bool target) {

	/// Find the hidden node values \x2

	/// Loop for each hidden nodes
	for (int i = 0; i < x2Size; i++) {

		/// Clear the accumulator
		double acc = 0;

		/// Weight and sum each input node
		for (int j = 0; j < x1Size; j++)
			acc += x1[j] * wh[i][j];

		/// Pass summed value through sigmoid
		x2[i] = 1 / (1 + exp(-acc));
	}

	/// Find the output value \a x3

	/// Clear the accumulator
	double acc = 0;

	/// Weight and sum each hidden node
	for (int i = 0; i < x2Size; i++)
		acc += x2[i] * wo[i];

	/// Pass summed value through sigmoid
	double x3 = 1 / (1 + exp(-acc));

	/// Find error for this data point

	/// Find the error
	double edata = target - x3;

	/// Give extra weights to targets
	if (target) edata *= 5;

	return edata;
}

/**
* Implementation
* @author Petar Jerčić
* @date 22/01/2013
* @todo N/A
*/
void	DSPStatistic::FindNeuralNetworkWeights(double* x1, int x1Size, double* x2, int x2Size, double x3, double** wh, double* wo, double edata, double mu) {

	/// Find new weights for hidden layer
	for (int i = 0; i < x2Size; i++)
		for (int j = 0; j < x1Size; j++) {

			double slopeo = x3 * (1 - x3);
			double slopeh = x2[i] * (1 - x2[i]);
			double dx3dw = x1[j] * slopeh * wo[i] * slopeo;
			wh[i][j] += + dx3dw * edata * mu;
		}

	/// Find new weights for output layer
	for (int i = 0; i < x2Size; i++) {

		double slopeo = x3 * (1 - x3);
		double dx3dw = x2[i] * slopeo;
		wo[i] += dx3dw * edata * mu;
	}
}

/**
* With other weights, the outputs might classify the objects as: metal or non-metal, biological or nonbiological, enemy or ally, etc. No algorithms, no rules, no procedures; only a relationship between the input and output dictated by the values of the weights selected.
* @author Petar Jerčić
* @date 22/01/2013
* @todo N/A
*/
void	DSPStatistic::NeuralNetwork() {

	/// Alocate the array for the input values
	double* x1 = new double[15];

	/// Alocate the array for the values exiting hidden layer
	double* x2 = new double[4];

	/// Alocate the array for the values exiting output layer
	double* x3 = new double[3];

	/// Alocate the array for the hidden layer weights
	double** wh = new double*[4];
	for(int i = 0; i < 4; i++)
		wh[i] = new double[5];

	/// Alocate the array for the output layer weights
	double** wo = new double*[2];
	for(int i = 0; i < 2; i++)
		wo[i] = new double[4];

	/// Load the data into \a x1
	// sub
	///Calculate the weights \a wh and \a wo
	// sub

	/// Find the hidden node values, \a x2
	/// Loop for each hidden layer node
	for (int i = 0; i < 4; i++) {

		/// Clear the accumulator
		double acc = 0;

		/// Weight and sum each input node
		for (int j = 0; j < 15; j++)
			acc += x1[j] * wh[i][j];

		/// Pass summed value through sigmoid
		x2[i] = 1 / (1 + exp(-acc));
	}

	/// Find the output value \a x3

	/// Loop for each hidden nodes
	for (int i = 0; i < 2; i++) {

		/// Clear the accumulator
		double acc = 0;

		/// Weight and sum each hidden node
		for (int j = 0; j < 4; j++)
			acc += x2[j] * wo[i][j];

		/// Pass summed value through sigmoid
		double x3 = 1 / (1 + exp(-acc));
	}
}

/**
* The array, T[ ], holds the desired frequency response, some kind of curve that we have manually designed. Since this program is based around the FFT, the lengths of the signals must be a power of two. As written, this program uses an FFT length of 256, as defined by the variable, N%, in line 130. This means that T[0] to T[128] correspond to the frequencies between 0 and 0.5 of the sampling rate. Only the magnitude is contained in this array; the phase is not controlled in this design, and becomes whatever it becomes.
* The recursion coefficients are set to their initial values in lines 270-310, typically selected to be the identity system. Don't use random numbers here, or the initial filter will be unstable. The recursion coefficients are held in the arrays, A[ ] and B[ ]. The variable, NP%, sets the number of poles in the designed filter. For example, if NP% is 5, the "a" coefficients run from A[0] to A[5], while the "b" coefficients run from B[1] to B[5]. 
* @author Petar Jerčić
* @date 23/01/2013
* @todo Implement saving coefficients
*/
void	DSPStatistic::IterativeRecursiveFilter() {

	/// Initialize

	/// Number of points in FFT
	int n = 256;

	/// Number of poles in the filter
	int np = 8;

	/// Perturbation increment
	double delta = 0.00001;

	/// Iteration step size
	double mu = 0.2;

	/// Alocate the array for the real part of the signal during FFT
	double* reX = new double[n];

	/// Alocate the array for the imaginary part of the signal during FFT
	double* imX = new double[n];

	/// Alocate the array for the desired frequency response
	double* t = new double[n / 2];

	/// Alocate the array for the \a a recursion coefficients
	double* a = new double[np];

	/// Alocate the array for the \a b recursion coefficients
	double* b = new double[np];

	/// Alocate the array for the slope of \a a recursion coefficients
	double* sa = new double[np];

	/// Alocate the array for the slope of \a b recursion coefficients
	double* sb = new double[np];

	/// Load \a t
	//sub

	/// Initialize coefficients to the identify system
	for (int i = 0; i < np; i++) {

		a[i] = 0;
		b[i] = 0;
	}
	a[0] = 1;

	/// Iteration loop

	/// Loop for desired number of iterations
	for (int i = 0; i < 100; i++) {

		double eold, enew;

		/// Calculate new coefficients
		CalculateRecursionCofficients(a, b, sa, sb, t, reX, imX, np, n, delta, mu, eold, enew);

		if (enew > eold)	mu /= 2;
	}

	/// Implement saving coefficients
}

/**
* Subroutine 2000 updates the recursion coefficients according to the steepest decent method: calculate the slope for each coefficient, and then change the coefficient an amount proportional to its slope. Lines 2080-2130 calculate the slopes for the "a" coefficients, storing them in the array, SA[ ]. Likewise, lines 2150-2200 calculate the slopes for the "b" coefficients, storing them in the array, SB[ ]. Lines 2220-2250 then modify each of the recursion coefficients by an amount proportional to these slopes. In this program, the proportionality constant is simply the step size, MU. No error term is required in the proportionality constant because there is only one example to be matched: the desired frequency response.
* The last issue is how the program calculates the slopes of the recursion coefficients. In the neural network example, an equation for the slope was derived. This procedure cannot be used here because it would require taking the derivative across the DFT. Instead, a brute force method is applied: actually change the recursion coefficient by a small increment, and then directly calculate the new value of ER. The slope is then found as the change in ER divided by the amount of the increment. Specifically, the current value of ER is found in lines 2040-2050, and stored in the variable, EOLD. The loop in lines 2080-2130 runs through each of the "a" coefficients. The first action inside this loop is to add a small increment, DELTA, to the recursion coefficient being worked on (line 2090). Subroutine 3000 is invoked in line 2100 to find the value of ER with the modified coefficient. Line 2110 then calculates the slope of this coefficient as: (ER - EOLD)/DELTA. Line 2120 then restores the modified coefficient by subtracting the value of DELTA.
* @author Petar Jerčić
* @date 23/01/2013
* @todo 
*/
void	DSPStatistic::CalculateRecursionCofficients(double* a, double* b, double* sa, double* sb, double* t, double* reX, double* imX, int noPoles, int noFFTPoints, double delta, double mu, double& eold, double& enew) {

	/// Find the current error
	//sub
	
	/// Store the current error in the variable \a eold
	eold = CalculateFrequencyDomainError(a, b, t, reX, imX, noPoles, noFFTPoints);

	/// Find the error slopes

	/// Loop through each \a a coefficient
	for (int i = 0; i < noPoles; i++) {

		/// Add the small increment to the coefficient
		a[i] += delta;

		/// Find the error with a chage
		double er = CalculateFrequencyDomainError(a, b, t, reX, imX, noPoles, noFFTPoints);

		/// Calculate the error slope, store in \a sa
		sa[i] = (er - eold) / delta;

		/// Return coefficients to oreginal value
		a[i] -= delta;
	}

	/// Repeat process for each \a b coefficient
	for (int i = 0; i < noPoles; i++) {

		/// Add the small increment to the coefficient
		b[i] += delta;

		/// Find the error with a chage
		double er = CalculateFrequencyDomainError(a, b, t, reX, imX, noPoles, noFFTPoints);

		/// Calculate the error slope, store in \a sb
		sb[i] = (er - eold) / delta;

		/// Return coefficients to oreginal value
		b[i] -= delta;
	}

	/// Calculate the new coefficients

	/// Loop through each coefficient
	for (int i = 0; i < noPoles; i++) {

		/// Adjust coefficients to move downhill
		a[i] -= sa[i] * mu;
		b[i] -= sb[i] * mu;
	}

	/// Find the new error and store in variable \a enew
	enew = CalculateFrequencyDomainError(a, b, t, reX, imX, noPoles, noFFTPoints);
}

/**
* The array, T[ ], holds the desired frequency response, some kind of curve that we have manually designed. Since this program is based around the FFT, the lengths of the signals must be a power of two.
* As previously mentioned, the iterative procedure requires a single value that describes how well the current system is functioning. This is provided by the variable, ER (for error), and is calculated in subroutine 3000. Lines  3040 to 3080 load an impulse in the array, IMX[ ]. Next, lines 3100-3150 use this impulse as an input signal to the recursive filter defined by the current values of A[ ] and B[ ]. The output of this filter is thus the impulse response of the current system, and is stored in the array, REX[ ]. The system's frequency response is then found by taking the FFT of the impulse response, as shown in line 3170. Subroutine 1000 is the FFT program listed in Table 12-4 in Chapter 12. This FFT subroutine returns the frequency response in rectangular form, overwriting the arrays REX[ ] and IMX[ ].
* Lines 3200-3250 then calculate ER, the mean squared error between the magnitude of the current frequency response, and the desired frequency response. Pay particular attention to how this error is found. The iterative action of this program optimizes this error, making the way it is defined very important. The FOR-NEXT loop runs through each frequency in the frequency response. For each frequency, line 3220 calculates the magnitude of the current frequency response from the rectangular data. In line 3230, the error at this frequency is found by subtracting the desired magnitude, T[ ], from the current magnitude, MAG. This error is then squared, and added to the accumulator variable, ER. After looping through each frequency, line 3250 completes the calculation to make ER the mean squared error of the entire frequency response.
* Lines 340 to 380 control the iteration loop of the program. Subroutine 2000 is where the changes to the recursion coefficients are made. The first action in this subroutine is to determine the current value of ER, and store it in another variable, EOLD (lines 2040 & 2050). After the subroutine updates the coefficients, the value of ER is again determined, and assigned to the variable, ENEW (lines 2270 and 2280).
* The variable, MU, controls the iteration step size, just as in the previous neural network program. An advanced feature is used in this program: an automated adjustment to the value of MU. This is the reason for having the two variables, EOLD and ENEW. When the program starts, MU is set to the relatively high value of 0.2 (line 160). This allows the convergence to proceed rapidly, but will limit how close the filter can come to an optimal solution. As the iterations proceed, points will be reached where no progress is being made, identified by ENEW being higher than EOLD. Each time this occurs, line 370 reduces the value of MU.
* @author Petar Jerčić
* @date 23/01/2013
* @todo 
*/
double	DSPStatistic::CalculateFrequencyDomainError(double* a, double* b, double* t, double* reX, double* imX, int noPoles, int noFFTPoints) {

	/// Load shifted impulse intu \a imX
	for (int i = 0; i < noFFTPoints; i++) {

		reX[i] = 0;
		imX[i] = 0;
	}
	imX[12] = 1;

	/// Calculate impulse response
	for (int i = 12; i < noFFTPoints; i++)
		for(int j = 0; j < noPoles; j++)
			reX[i] += a[j] * imX[i - j] + b[j] * reX[i - j];
	imX[12] = 1;

	/// Calculate the FFT
	FastFourierTransform(noFFTPoints, reX, imX);

	/// Find the frequency domain error

	/// Zero ER to use as an accumulator
	double er = 0;

	/// Loop through each positive frequency
	for (int i = 0; i < noFFTPoints / 2; i++) {

		/// Rectangular polar conversion
		double mag = sqrt(pow(reX[i], 2) + pow(imX[i], 2));

		/// Calculate and accumulate squared error
		er += pow(mag - t[i], 2);
	}

	/// Finish calculation of error \a er
	er = sqrt(er / (noFFTPoints / 2 + 1));

	return er;
}

/**
* This program changes the recursive coeff from each of the individual stages into transfer function of Eq. 31-3. Each stage is associated with coeff in arrays and we can manipulate polynomials just by manipulating their coeff's. When two poly's are multiplied, their coeff's are convolved. Do this stepwise for cascaded and parallel.
* @author Petar Jerčić
* @date 26/02/2013
* @todo N/A
*/
void	DSPStatistic::CombineRecursionCoeffCascadeParallelStages (double* a1, int a1Size, double* b1, int b1Size, double* a2, int a2Size,  double* b2, int b2Size, double** a3, int& a3Size, double** b3, int& b3Size, bool bParallel) {

	/// INITIALIZE VARIABLES

	/// Alocate the array for \a a and \a b coeff for sys 3, the combined system
	a3Size = a1Size + a2Size;
	b3Size = b1Size + b2Size;
	double* a3_ = new double[a3Size];
	double* b3_ = new double[b3Size];

	/// Convert the recursion coeff into transfer function
	for (int i = 0; i < 8; i++) {

		b2[i] = b2[i];
		b1[i] = - b1[i];
	}
	b1[0] = 1;
	b2[0] = 1;

	/// Multiply the polynomials by convolving
	for (int i = 0; i < 16; i++) {

		a3_[i] = 0;
		b3_[i] = 0;

		for (int j = 0; j < 8; j++) {

			if (i - j >= 0 || i - j < 8) {

				if (!bParallel)	a3_[i] += a1[j] * a2[i - j];
				if (bParallel)	a3_[i] += a1[j] * b2[i - j] + a2[j] * b1[i - j];
				b3_[i] += b1[j] * b2[i - j];
			}
		}
	}

	/// Convert the transfer function int recursion coeff.
	for (int i = 0; i < 16; i++)
		b3_[i] = -b3_[i];
	b3_[0] = 0;

	/// The recursion coeff of the combined system now reside in \a a3_[] and \a b3_[]
	*a3 = a3_;
	*b3 = b3_;
}

/**
* This is a leaky integrator
* @author Petar Jerčić
* @date 23/03/2013
* @todo N/A
*/
double	DSPStatistic::LeakyIntegrator(double x) {

	static const double lambda = 0.9;
	static double y = 0;

	y = lambda * y + (1 - lambda) * x;

	return y;
}

/**
* This is an integrator
* @author Petar Jerčić
* @date 26/02/2013
* @todo N/A
*/
double*	DSPStatistic::Integrator(double* x, int size) {

	double* y = new double[size];

	for(int i = 0; i < size; i++)
		y[i] = LeakyIntegrator(x[i]);

	return y;
}

/**
* This is a moving average
* @author Petar Jerčić
* @date 23/03/2013
* @todo N/A
*/
double	DSPStatistic::MovingAvegare(double x){

	static const int M = 10;
	static double z[M];
	static int ix = -1;

	int n;
	double avg = 0;

	if (ix == -1) {
		for(n =0; n < M; n++)
			z[n] = 0;
		ix = 0;
	}

	z[ix] = x;
	ix = (ix + 1) % M;

	for (n = 0; n < M; n++)
		avg += z[n];

	return avg / M;
}

/**
* This is an average
* @author Petar Jerčić
* @date 23/03/2013
* @todo N/A
*/
double*	DSPStatistic::Average(double* x, int size) {

	double* y = new double[size];

	for(int i = 0; i < size; i++)
		y[i] = MovingAvegare(x[i]);

	return y;
}

/**
* This is an average
* @author Petar Jerčić
* @date 24/03/2013
* @todo N/A
*/
void		DSPStatistic::Callback(const void* input, void* output, unsigned long samples) {

	/// Cast pointers to convert them to the right data type
	float*	pIn = (float *) input;
	float*	pOut = (float *) output;

	/// Call processing function for every data point in the input buffer
	for (int n = 0; n < samples; n++)
		*pOut++ = LeakyIntegrator(*pIn++);
}

/**
* This is an echo
* @author Petar Jerčić
* @date 24/03/2013
* @todo N/A
*/
float	DSPStatistic::Echo(double* signal, int size, int sampling_freq, int index) {

	/// Create a buffer mask
	//enum{BUFF_MASK = size -1};
	enum{BUFF_MASK = 0x10000 -1};

	/// Coaefs for the three replicas of the signal
	static float a = 1;
	static float b = 0.7f;
	static float c = 0.5f;

	/// Normalizing coeffs
	static float norm = 1.0f / (a + b + c);
	/// Introduce a delay based on sampling frequency
	static int N = (int)(0.3 * sampling_freq);

	/// Weight sum of the three delayed replicas of the input
	return norm * ( a * signal[index]
					+ b * signal[(index + size - N)		& BUFF_MASK]
					+ c * signal[(index + size - 2 * N)	& BUFF_MASK]);
}

/**
* This is an echo
* @author Petar Jerčić
* @date 24/03/2013
* @todo N/A
*/
float	DSPStatistic::NaturalEcho(double* signal, double* processed_signal, int size, int sampling_freq, int index_x, int index_y) {

	/// Create a buffer mask
	//enum{BUFF_MASK = size -1};
	enum{BUFF_MASK = 0x10000 -1};

	/// Decaying factor
	static float a = 0.7;

	/// Factor for the leaky integrator (loww pass)
	static float L = 0.6;
	static float c = 0.5f;

	/// Normalizing coeffs
	static float norm = 1.0f / (1 + a);

	/// Introduce a delay based on sampling frequency
	static int N = (int)(0.3 * sampling_freq);

	/// Weight sum of the three delayed replicas of the input
	return norm * ( signal[index_x]
					- L * signal[(index_x + size - 1)		& BUFF_MASK]
					+ L * processed_signal[(index_y + size - 1)		& BUFF_MASK]
					+ a * (1 - L) * processed_signal[(index_y + size - N)	& BUFF_MASK]);
}

/**
* This is an reverb
* @author Petar Jerčić
* @date 24/03/2013
* @todo N/A
*/
float	DSPStatistic::Reverb(double* signal, double* processed_signal, int size, int sampling_freq, int index_x, int index_y) {

	/// Create a buffer mask
	//enum{BUFF_MASK = size -1};
	enum{BUFF_MASK = 0x10000 -1};

	/// Decaying factor
	static float a = 0.8;

	/// Introduce a delay based on sampling frequency
	static int N = (int)(0.02 * sampling_freq);

	/// Weight sum of the three delayed replicas of the input
	return			( - a * signal[index_x]
					- signal[(index_x + size - N)		& BUFF_MASK]
					+ a * processed_signal[(index_y + size - N)	& BUFF_MASK]);
}

/**
* This is an reverb
* @author Petar Jerčić
* @date 24/03/2013
* @todo N/A
*/
float	DSPStatistic::Fuzz(double* signal, int index_x) {


	static float limit = 0.005;
	static float G = 5;

	float y = signal[index_x];

	if (y > limit)
		y = limit;
	if (y < -limit)
		y = - limit;

	return G * y;
}

/**
* This is an tremolo
* @author Petar Jerčić
* @date 24/03/2013
* @todo N/A
*/
float	DSPStatistic::Tremolo(double* signal, int sampling_freq, int index_x) {

	double PI = 4 * atan(1.0);

	/// Speed of the modulated sinusoide
	static double phi = 5 * 2 * PI / sampling_freq;	/// 5Hz LFO
	static double omega = 0;

	omega += phi;

	return (1 + cos(omega) / 4) * signal[index_x];
}

/**
* This is an flan
* @author Petar Jerčić
* @date 24/03/2013
* @todo N/A
*/
float	DSPStatistic::Flanger(double* signal, int size, int sampling_freq, int index_x) {

	double PI = 4 * atan(1.0);

	/// Create a buffer mask
	//enum{BUFF_MASK = size -1};
	enum{BUFF_MASK = 0x10000 -1};

	/// wax delay
	static int N = (int)(0.002 * sampling_freq);

	/// Speed of the modulated sinusoide
	static double phi = 0.1 * 2 * PI / sampling_freq;	/// 0.1Hz LFO
	static double omega = 0;

	int D = (int)(N * (1 + cos(omega)));
	omega += phi;

	return 0.5f + (signal[index_x] + signal[(index_x + size - D)		& BUFF_MASK]);
}
