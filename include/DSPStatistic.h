/**
* @file DSPStatistic.h
* @brief This header file defines an Analysis component using methods of DSP on signals.
*
* @author Petar Jerčić
*
* @date 09/12/2012
*/

#ifndef DSP_STASTIC
#define DSP_STASTIC

#include <iostream>
#include "IInput.h"
#include "IOutput.h"
#include "IAnalysis.h"

class DSPStatistic : public IAnalysis {
public:

	/**
	* This method will be used to create random noise signal for other components to use.
	* @author Petar Jerčić
	* @param fileName	N/A
	* @param channel	N/A
	* @return	N/A
	* @sa	something
	* @date 08/12/2012
	*/
	void setInput(IInput* input);

	void setOutput(IOutput* output);

	void Process();

	DSPStatistic();

	DSPStatistic(IInput*input, IOutput* output);

public:
	IInput* _input;
	IOutput* _output;
	double* signal;
	int		noSamples;
	int		smpFreq;

	void ProcesInput();

	void ProcesSignal();

	void	CalculateDescriptivesStatic(double* sig, int noSamp, double& avg, double& stdev, double& var, double &snr, double &cv);
	void	CalculateDescriptivesRunning(double data, double& avg, double& stdev, double& var, double &snr, double &cv);
	deque<int>	CalculateDescriptivesStaticHistogram(double* sig, int noSamp, double& avg, double& stdev, double& var, double &snr, double &cv);
	void	ConvolutionInputSide(double* sig, int noSamp, double* impulseResponse, int irSize, double** outputSignal, int& osSize);
	void	ConvolutionOutputSide(double* sig, int noSamp, double* impulseResponse, int irSize, double** outputSignal, int& osSize);
	void	FirstDifference(double* sig, int noSamp, double** outputSignal, int& osSize);
	void	RunningSumIntegration(double* sig, int noSamp, double** outputSignal, int& osSize);

	/**
	* Decompose a time domain signal into sinusoids in frequency domain signals from time side.
	* Frequency domain signals, held in \a realX and \a imaginaryX, are calculated from the time domain signal, held in \a timeSignal.
	* @param [in]	timeSignal		Will hold the time domain signal
	* @param [in]	timeSigSize		Will hold the time domain signal size
	* @param [out]	realX			Holds the real part of the frequency domain
	* @param [out]	reXSize			Holds the real part of the frequency domain size
	* @param [out]	imaginaryX		Holds the imaginary part of the frequency domain
	* @param [out]	imXSize			Holds the imaginary part of the frequency domain size
	*/
	void	DiscreteFourierTransformationTimeSide(double* timeSignal, int timeSigSize, double** realX, int& reXSize, double** imaginaryX, int& imXSize);

	/**
	* Decompose a time domain signal into sinusoids in frequency domain signals from frequency side.
	* Frequency domain signals, held in \a realX and \a imaginaryX, are calculated from the time domain signal, held in \a timeSignal.
	* @param [in]	timeSignal		Will hold the time domain signal
	* @param [in]	timeSigSize		Will hold the time domain signal size
	* @param [out]	realX			Holds the real part of the frequency domain
	* @param [out]	reXSize			Holds the real part of the frequency domain size
	* @param [out]	imaginaryX		Holds the imaginary part of the frequency domain
	* @param [out]	imXSize			Holds the imaginary part of the frequency domain size
	*/
	void	DiscreteFourierTransformationFreqSide(double* timeSignal, int timeSigSize, double** realX, int& reXSize, double** imaginaryX, int& imXSize);

	///Synthesis of time domain signal from frequency domain signals from frequency side.
	///Time domain signal, held in \a timeSignal, is calculated from the frequency domain signal, held in \a realX and \a imaginaryX
	///@param [in]	realX			Holds the real part of the frequency domain
	///@param [in]	reXSize			Holds the real part of the frequency domain size
	///@param [in]	imaginaryX		Holds the imaginary part of the frequency domain
	///@param [in]	imXSize			Holds the imaginary part of the frequency domain size
	///@param [out]	timeSignal		Will hold the time domain signal
	///@param [out]	timeSigSize		Will hold the time domain signal size
	void	InverseDiscreteFourierTransformationFreqSide(double* realX, int reXSize, double* imaginaryX, int imXSize, double** timeSignal, int& timeSigSize);

	/**
	* Synthesis of time domain signal from frequency domain signals from time side.
	* Time domain signal, held in \a timeSignal, is calculated from the frequency domain signal, held in \a realX and \a imaginaryX
	* @param [in]	realX			Holds the real part of the frequency domain
	* @param [in]	reXSize			Holds the real part of the frequency domain size
	* @param [in]	imaginaryX		Holds the imaginary part of the frequency domain
	* @param [in]	imXSize			Holds the imaginary part of the frequency domain size
	* @param [out]	timeSignal		Will hold the time domain signal
	* @param [out]	timeSigSize		Will hold the time domain signal size
	*/
	void	InverseDiscreteFourierTransformationTimeSide(double* realX, int reXSize, double* imaginaryX, int imXSize, double** timeSignal, int& timeSigSize);

	/**
	* Rectangular to polar Conversion
	* From real and imginary to magnitude and phase
	* @param [in]	realX			Holds the real part of the frequency domain
	* @param [in]	reXSize			Holds the real part of the frequency domain size
	* @param [in]	imaginaryX		Holds the imaginary part of the frequency domain
	* @param [in]	imXSize			Holds the imaginary part of the frequency domain size
	* @param [out]	mag				Holds the magnitude
	* @param [out]	magSize			Holds the magnitude size
	* @param [out]	phase			Holds the phase
	* @param [out]	phaseSize		Holds the phase size
	*/
	void	RectangularToPolarConversion(double* realX, int reXSize, double* imaginaryX, int imXSize, double** mag, int& magSize, double** phase, int& phaseSize);

	/**
	* Polar To Rectangular Conversion
	* From magnitude and phase to real and imginary
	* @param [in]	mag				Holds the magnitude
	* @param [in]	magSize			Holds the magnitude size
	* @param [in]	phase			Holds the phase
	* @param [in]	phaseSize		Holds the phase size
	* @param [out]	realX			Holds the real part of the frequency domain
	* @param [out]	reXSize			Holds the real part of the frequency domain size
	* @param [out]	imaginaryX		Holds the imaginary part of the frequency domain
	* @param [out]	imXSize			Holds the imaginary part of the frequency domain size
	*/
	void	PolarToRectangularConversion(double* mag, int magSize, double* phase, int phaseSize, double** realX, int& reXSize, double** imaginaryX, int& imXSize);

	/**
	* Unwrapping the phase
	* It is often easier to understand the phase if it does not have these discontinuities, even if it means that the phase extends above π, or below -π.
	* @param [in]	phase			Holds the phase
	* @param [in]	phaseSize		Holds the phase size
	* @param [out]	uwPhase			Holds the unwrapped phase
	* @param [out]	uwPhaseSize		Holds the unwrapped phase size
	*/
	void	PhaseUnwrapping(double* phase, int phaseSize, double** uwPhase, int& uwPhaseSize);

	/**
	* @brief This subroutine creates the complex frequency domain from the real frequency domain.
	*
	* Helper method generating negative frequencies between samples N/2 + 1 and N - 1.
	* @param [in]	realFreqSize	Holds the size of the signals array
	* @param [out]	reX				Holds the real frequency domain
	* @param [out]	imX				Holds the imaginary frequency domain
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 12 > Table 12-1
	*/
	void	NegativeFrequencyGeneration(int realFreqSize, double* reX, double* imX);

	/**
	* @brief This subroutine calculates the complex DFT. 
	*
	* It does so by correlating the time domain signal with sine and cosine waves. The algorithm is called decimation in time.
	* @param [in]	signalsSize		Holds the size of the signals array
	* @param [in]	xR				Holds the real part of the time domain
	* @param [in]	xI				Holds the imaginary part of the time domain
	* @param [out]	reX				Holds the real frequency domain
	* @param [out]	imX				Holds the imaginary frequency domain
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 12 > Table 12-2
	*/
	void	ComplexDFTByCorrelationTime(int signalsSize, double* xR, double* xI, double* reX, double* imX);

	/**
	* @brief This subroutine calculates the Fast Fourier Transform.
	*
	* It does so by correlating the time domain signal with sine and cosine waves. The algorithm is called decimation in frequency.
	* @param [in]	signalsSize		Holds the size of the signals array
	* @param [out]	reX				Holds the real part of the input
	* @param [out]	imX				Holds the imaginary part of the input
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 12 > Table 12-4
	*/
	void	FastFourierTransform(int signalsSize, double* reX, double* imX);

	/**
	* @brief This subroutine calculates the Inverse Fast Fourier Transform.
	*
	* The easiest way to calculate an Inverse FFT is to calculate a Forward FFT, and then adjust the data.
	* @param [in]	signalsSize		Holds the size of the signals array for Inverse Fast Fourier Transform
	* @param [out]	reX				Holds the real part of the complex frequency domain
	* @param [out]	imX				Holds the imaginary part of the complex frequency domain
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 12 > Table 12-5
	*/
	void	InverseFastFourierTransform(int signalsSize, double* reX, double* imX);

	/**
	* @brief This subroutine calculates the Fast Fourier Transform for real numbers.
	*
	* There are two small disadvantages in using the real FFT. First, the code is about twice as long. While your computer doesn't care, you must take the time to convert someone else's program to run on your computer. Second, debugging these programs is slightly harder because you cannot use symmetry as a check for proper operation. These algorithms force the imaginary part of the time domain to be zero, and the frequency domain to have left-right symmetry. For debugging, check that these programs produce the same output as the conventional FFT algorithms. 
	* @param [in]	signal			Holds the signals array
	* @param [in]	signalSize		Holds the size of the signals array, has to be even number
	* @param [out]	reX				Will Hold the real part of the input
	* @param [out]	imX				Will Hold the imaginary part of the input
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 12 > Table 12-7
	*/
	void	FastFourierTransformReal(double* signal, int signalSize, double** reX, double** imX);

	/**
	* @brief This subroutine calculates the Inverse Fast Fourier Transform for real numbers.
	*
	* The easiest way to calculate an Inverse FFT is to calculate a Forward FFT, and then adjust the data.
	* @param [in]	signalsSize		Holds the size of the signals array for Inverse Fast Fourier Transform
	* @param [in]	reX				Holds the real part of the complex frequency domain
	* @param [in]	imX				Holds the imaginary part of the complex frequency domain
	* @return A double pointer holding the real time domain.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 12 > Table 12-6
	*/
	double*	InverseFastFourierTransformReal(int signalsSize, double* reX, double* imX);

	/**
	* @brief  The moving average filter operates by averaging a number of points from the input signal to produce each point in the output signal.
	*
	* You should recognize that the moving average filter is a convolution using a very simple filter kernel. 
	* For example, a 5 point filter has the filter kernel: …0, 0, 1/5, 1/5, 1/5, 1/5, 1/5, 0, 0…. 
	* That is, the moving average filter is a convolution of the input signal with a rectangular pulse having an area of one.
	* @param [in]	x		Holds the input signal
	* @param [in]	xSize	Holds the input signal size
	* @return A double pointer holding the output signal.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 15 > Table 15-1
	*/
	double*	MovingAverageFilter(double* x, int xSize);

	/**
	* @brief  The moving average filter operates by averaging a number of points from the input signal to produce each point in the output signal.
	*
	* This algorithm is faster than other digital filters for several reasons. 
	* First, there are only two computations per point, regardless of the length of the filter kernel. 
	* Second, addition and subtraction are the only math operations needed, while most digital filters require time-consuming multiplication. 
	* Third, the indexing scheme is very simple.
	* @param [in]	x		Holds the input signal
	* @param [in]	xSize	Holds the input signal size
	* @return A double pointer holding the output signal.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 15 > Table 15-2
	*/
	double*	MovingAverageFilterRecursion(double* x, int xSize);

	/**
	* @brief  Windowed-sinc filters are used to separate one band of frequencies from another. Get high-pass filter with simple spectral inversion of the low-pass filter kernel (check BandPassWindowdSincFilter).
	*
	* To design a windowed-sinc, two parameters must be selected: the cutoff frequency, fc, and the length of the filter kernel, M. 
	* The cutoff frequency is expressed as a fraction of the sampling rate, and therefore must be between 0 and 0.5.
	* @param [in]	x		Holds the input signal
	* @param [in]	xSize	Holds the input signal size
	* @return A double pointer holding the output signal.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 16 > Table 16-1
	*/
	double*	LowPassWindowdSincFilter(double* x, int xSize);

	/**
	* @brief   Our goal is to separate the alpha from the beta rhythms. To do this, we will design a digital low-pass filter with a cutoff frequency of 14 hertz, or 0.14 of the sampling rate.
	*
	* If you close your eyes and relax, the predominant EEG pattern will be a slow oscillation between about 7 and 12 hertz. This waveform is called the alpha rhythm, and is associated with contentment and a decreased level of attention. 
	* In this example, we will assume that the EEG signal has been amplified by analog electronics, and then digitized at a sampling rate of 100 samples per second.
	* @param [out]	fiterKernel		Will hold the band pass filter kernel
	* @param [out]	fiterKernelSize	Will hold the band pass filter kernel size
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 16 > Table 16-2
	*/
	void	BandPassWindowdSincFilter(double** fiterKernel, int& fiterKernelSize);

	/**
	* @brief   Create a band-pass filter kernel of arbitrary size
	*
	* If you close your eyes and relax, the predominant EEG pattern will be a slow oscillation between about 7 and 12 hertz. This waveform is called the alpha rhythm, and is associated with contentment and a decreased level of attention. 
	* In this example, we will assume that the EEG signal has been amplified by analog electronics, and then digitized at a sampling rate of 100 samples per second.
	* @param [in]	fiterKernelSize	Holds the band-pass filter kernel size, has to be odd number
	* @return A double pointer holding the band pass filter kernel
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 16 > Table 16-2
	*/
	double*	BandPassWindowdSincFilterCustom(int fiterKernelSize);

	/**
	* @brief   Convert the aliased impulse response into a filter kernel
	*
	* @param [in]	reX		Holds the signal being converted
	* @param [in]	reXSize	Holds the signal size
	* @return A double pointer holding the output signal.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 17 > Table 17-1
	*/
	double*	CustomFilterDesign(double* reX, int reXSize);

	/**
	* @brief   Method to carry out FFT convolution.
	*
	* FFT convolution uses the principle that multiplication in the frequency domain corresponds to convolution in the time domain. The input signal is transformed into the frequency domain using the DFT, 
	* multiplied by the frequency response of the filter, and then transformed back into the time domain using the Inverse DFT.
	* Keep FFT convolution tucked away for when you have a large amount of data to process and need an extremely long filter kernel. Think in terms of a million sample signal and a thousand point filter kernel.
	* @param [in]	x		Holds the input signal
	* @param [in]	xSize	Holds the input signal size
	* @return A double pointer holding the output signal.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 18 > Table 18-1
	*/
	double*	FFTConvolution(double* x, int xSize);

	/**
	* @brief   Recursive filters are an efficient way of achieving a long impulse response, without having to perform a long convolution.
	*
	* They execute very rapidly, but have less performance and flexibility than other digital filters. Recursive filters are also called Infinite Impulse Response (IIR) filters, since their impulse responses are composed of decaying exponentials.
	* @param [in]	x		Holds the input signal
	* @param [in]	xSize	Holds the input signal size
	* @param [in]	a		Holds the a recursion coefficients
	* @param [in]	aSize	Holds the a recursion coefficients array size
	* @param [in]	b		Holds the b recursion coefficients
	* @param [in]	bSize	Holds the b recursion coefficients array size
	* @return A double pointer holding the output signal.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 19 > Table 19-1
	*/
	double*	RecursiveFilter(double* x, int xSize, double* a, int aSize, double* b, int bSize);

	/**
	* @brief   SinglePoleLowPass Recursive filters are an efficient way of achieving a long impulse response, without having to perform a long convolution.
	*
	* These single pole recursive filters are definitely something you want to keep in your DSP toolbox. You can use them to process digital signals just as you would use RC networks to process analog electronic signals. This includes everything you would expect: DC removal, high-frequency noise suppression, wave shaping, smoothing, etc.
	* @param [in]	signal		Holds the input signal
	* @param [in]	signalSize	Holds the input signal size
	* @param [in]	cutoffFreq	Holds the cutoff frequency of the filter
	* @return A double pointer holding the output signal.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 19 > Table 19-1
	*/
	double*	SinglePoleLowPassRecursiveFilter(double* signal, int signalSize, double cutoffFreq);

	/**
	* @brief   SinglePoleHighPass Recursive filters are an efficient way of achieving a long impulse response, without having to perform a long convolution.
	*
	* These single pole recursive filters are definitely something you want to keep in your DSP toolbox. You can use them to process digital signals just as you would use RC networks to process analog electronic signals. This includes everything you would expect: DC removal, high-frequency noise suppression, wave shaping, smoothing, etc.
	* @param [in]	signal		Holds the input signal
	* @param [in]	signalSize	Holds the input signal size
	* @param [in]	cutoffFreq	Holds the cutoff frequency of the filter
	* @return A double pointer holding the output signal.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 19 > Table 19-1
	*/
	double*	SinglePoleHighPassRecursiveFilter(double* signal, int signalSize, double cutoffFreq);

	/**
	* @brief   The four stage low-pass filter is comparable to the Blackman and Gaussian filters (relatives of the moving average, Chapter 15), but with a much faster execution speed.
	*
	* These single pole recursive filters are definitely something you want to keep in your DSP toolbox. You can use them to process digital signals just as you would use RC networks to process analog electronic signals. This includes everything you would expect: DC removal, high-frequency noise suppression, wave shaping, smoothing, etc.
	* @param [in]	signal		Holds the input signal
	* @param [in]	signalSize	Holds the input signal size
	* @param [in]	cutoffFreq	Holds the cutoff frequency of the filter
	* @return A double pointer holding the output signal.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 19 > Table 19-1
	*/
	double*	FourStageLowPassRecursiveFilter(double* signal, int signalSize, double cutoffFreq);

	/**
	* @brief   Band pass filter A common need in electronics and DSP is to isolate a narrow band of frequencies from a wider bandwidth signal.
	*
	* These single pole recursive filters are definitely something you want to keep in your DSP toolbox. You can use them to process digital signals just as you would use RC networks to process analog electronic signals. This includes everything you would expect: DC removal, high-frequency noise suppression, wave shaping, smoothing, etc.
	* @param [in]	signal		Holds the input signal
	* @param [in]	signalSize	Holds the input signal size
	* @param [in]	centerFreq	Holds the center frequency as fraction of a sampling rate 0 - 0.5
	* @param [in]	bandwidth	Holds the bandwidth of a filter as fraction of a sampling rate 0 - 0.5
	* @return A double pointer holding the output signal.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 19 > Eq 19-7
	*/
	double*	BandPassRecursiveFilter(double* signal, int signalSize, double centerFreq, double bandwidth);

	/**
	* @brief   Band reject filter A common need in electronics and DSP is to isolate a narrow band of frequencies from a wider bandwidth signal.
	*
	* These single pole recursive filters are definitely something you want to keep in your DSP toolbox. You can use them to process digital signals just as you would use RC networks to process analog electronic signals. This includes everything you would expect: DC removal, high-frequency noise suppression, wave shaping, smoothing, etc.
	* @param [in]	signal		Holds the input signal
	* @param [in]	signalSize	Holds the input signal size
	* @param [in]	centerFreq	Holds the center frequency as fraction of a sampling rate 0 - 0.5
	* @param [in]	bandwidth	Holds the bandwidth of a filter as fraction of a sampling rate 0 - 0.5
	* @return A double pointer holding the output signal.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 19 > Eq 19-8
	*/
	double*	BandRejectRecursiveFilter(double* signal, int signalSize, double centerFreq, double bandwidth);

	/**
	* @brief   The Chebyshev filter Recursion Coefficients
	*
	* The Chebyshev response is a mathematical strategy for achieving a faster roll-off by allowing ripple in the frequency response.
	* As the ripple increases (bad), the roll-off becomes sharper (good). The Chebyshev response is an optimal trade-off between these two parameters. When the ripple is set to 0%, the filter is called a maximally flat or Butterworth filter (after S. Butterworth, a British engineer who described this response in 1930). A ripple of 0.5% is a often good choice for digital filters. This matches the typical precision and accuracy of the analog electronics that the signal has passed through.
	* ou must select four parameters to design a Chebyshev filter: (1) a high-pass or low-pass response, (2) the cutoff frequency, (3) the percent ripple in the passband, and (4) the number of poles. Just what is a pole?
	* Number of poles haave to be even.aa
	* @param [in]	signal		Holds the input signal
	* @param [in]	signalSize	Holds the input signal size
	* @param [in]	centerFreq	Holds the center frequency as fraction of a sampling rate 0 - 0.5
	* @param [in]	bandwidth	Holds the bandwidth of a filter as fraction of a sampling rate 0 - 0.5
	* @return A double pointer holding the output signal.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 19 > Eq 19-8
	*/
	void	ChebyshevFilterRecursionCoefficients(double cutoffFreq, bool lowPass, int percRipple, int noPoles, double** _a, double** _b);

	/**
	* @brief   The Chebyshev filter parametars Coefficients
	*
	* These single pole recursive filters are definitely something you want to keep in your DSP toolbox. You can use them to process digital signals just as you would use RC networks to process analog electronic signals. This includes everything you would expect: DC removal, high-frequency noise suppression, wave shaping, smoothing, etc.
	* @param [in]	signal		Holds the input signal
	* @param [in]	signalSize	Holds the input signal size
	* @param [in]	centerFreq	Holds the center frequency as fraction of a sampling rate 0 - 0.5
	* @param [in]	bandwidth	Holds the bandwidth of a filter as fraction of a sampling rate 0 - 0.5
	* @return A double pointer holding the output signal.
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 19 > Eq 19-8
	*/
	void	ChebyshevFilterParametars(double cutoffFreq, bool lowPass, int percRipple, int noPoles, int count, double ** _a, double **_b);

	/**
	* @brief   Determination of weights
	*
	* And save them
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 26 > Eq 26-2
	*/
	void	NeuralNetworkTraining();

	/**
	* @brief   The Chebyshev filter parametars Coefficients
	*
	* These single pole recursive filters are definitely something you want to keep in your DSP toolbox. You can use them to process digital signals just as you would use RC networks to process analog electronic signals. This includes everything you would expect: DC removal, high-frequency noise suppression, wave shaping, smoothing, etc.
	* @param [out]	x1		Holds the input signal, the bias node \a x1[last] asways has a value of one
	* @param [in]	x1Size	Holds the input signal size
	* @param [in]	dataSet	Holds the identifier for the data set we want to train
	* @return Bool value specifying if the training set is target or not-target data
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 26 > Eq 26-3
	*/
	bool	LoadNeuralNetworkTrainingSet(double* x1, int x1Size, int dataSet);

	/**
	* @brief   Calculate the error with the current weights
	*
	* subroutine 2000 passes the data through the current neural network to produce the output node value, X3. In other words, subroutine 2000 is the same as the program in Table 26-1, except for a different number of nodes in each layer. This subroutine also calculates how well the current network identifies the letter as a target or a nontarget.
	* @param [in]	x1		Holds the input signal, the bias node \a x1[last] asways has a value of one
	* @param [in]	x1Size	Holds the input signal size
	* @param [out]	x2		Holds the input signal, the bias node \a x1[last] asways has a value of one
	* @param [in]	x2Size	Holds the input signal size
	* @param [out]	x3		double
	* @param [in]	wh		Holds the input signal, the bias node \a x1[last] asways has a value of one
	* @param [in]	wo		Holds the input signal, the bias node \a x1[last] asways has a value of one
	* @param [in]	target	Bool value specifying if the training set is target or not-target data
	* @return Double value how well the current network identifies the letter as a target or a nontarget. This makes return a value between -1 and 1
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 26 > Eq 26-3
	*/
	double	CalculateNeuralNetworkError(double* x1, int x1Size, double* x2, int x2Size, double* x3, double** wh, double* wo, bool target);

	/**
	* @brief   Find new weights
	*
	* Subroutine 3000 is the heart of the neural network strategy, the algorithm for changing the weights on each iteration. We will use an analogy to explain the underlying mathematics. Consider the predicament of a military paratrooper dropped behind enemy lines. He parachutes to the ground in unfamiliar territory, only to find it is so dark he can't see more than a few feet away. His orders are to proceed to the bottom of the nearest valley to begin the remainder of his mission. The problem is, without being able to see more than a few feet, how does he make his way to the valley floor? Put another way, he needs an algorithm to adjust his x and y position on the earth's surface in order to minimize his elevation. This is analogous to the problem of adjusting the neural network weights, such that the network's error, ESUM, is minimized.
	* @param [in]	x1		Holds the input signal, the bias node \a x1[last] asways has a value of one
	* @param [in]	x1Size	Holds the input signal size
	* @param [in]	x2		Holds the input signal, the bias node \a x1[last] asways has a value of one
	* @param [in]	x2Size	Holds the input signal size
	* @param [out]	wh		Holds the input signal, the bias node \a x1[last] asways has a value of one
	* @param [out]	wo		Holds the input signal, the bias node \a x1[last] asways has a value of one
	* @param [in]	edata	Double value how well the current network identifies the letter as a target or a nontarget. This makes \a edata a value between -1 and 1
	* @param [in]	mu		Double value
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 26 > Eq 26-3
	*/
	void	FindNeuralNetworkWeights(double* x1, int x1Size, double* x2, int x2Size, double x3, double** wh, double* wo, double edata, double mu);

	/**
	* @brief   Neural Network
	*
	* Table 26-1 is a program to carry out the flow diagram of Fig. 26-5. The key point is that this architecture is very simple and very generalized. This same flow diagram can be used for many problems, regardless of their particular quirks. The ability of the neural network to provide useful data manipulation lies in the proper selection of the weights. This is a dramatic departure from conventional information processing where solutions are described in step-by-step procedures.
	* As an example, imagine a neural network for recognizing objects in a sonar signal. Suppose that 1000 samples from the signal are stored in a computer. How does the computer determine if these data represent a submarine, whale, undersea mountain, or nothing at all? Conventional DSP would approach this problem with mathematics and algorithms, such as correlation and frequency spectrum analysis. With a neural network, the 1000 samples are simply fed into the input layer, resulting in values popping from the output layer. By selecting the proper weights, the output can be configured to report a wide range of information. For instance, there might be outputs for: submarine (yes/no), whale (yes/no), undersea mountain (yes/no), etc. 
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 26 > Eq 26-1
	*/
	void	NeuralNetwork();

	/**
	* @brief   Iterative Recursive Filter
	*
	* Start with a generic set of recursion coefficients, and use iteration to slowly mold them into what you want. This technique is important for two reasons. First, it allows custom recursive filters to be designed without having to hassle with the mathematics of the z-transform. Second, it shows that the ideas from conventional DSP and neural networks can be combined to form superb algorithms. 
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 26 > Eq 26-4
	*/
	void	IterativeRecursiveFilter();

	/**
	* @brief   Calculate the new recursion cofficients
	*
	* general
	* @param [out]	a			Holds the \a a recursion coefficients
	* @param [out]	b			Holds the \a b recursion coefficients
	* @param [in]	sa			Holds the slope of \a a recursion coefficients
	* @param [in]	sb			Holds the slope of \a b recursion coefficients
	* @param [in]	t			Holds the array for the desired frequency response
	* @param [in]	reX			Holds the real part of the signal during FFT
	* @param [in]	imX			Holds the imaginary part of the signal during FFT
	* @param [in]	noPoles		Holds the Number of poles in the filter
	* @param [in]	noFFTPoints	Holds the Number of points in FFT
	* @param [in]	delta		Holds the Perturbation increment
	* @param [in]	mu			Holds the Iteration step size
	* @param [out]	eold		Holds the old error coeff
	* @param [out]	enew		Holds the new error coeff
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 26 > Eq 26-5
	*/
	void	CalculateRecursionCofficients(double* a, double* b, double* sa, double* sb, double* t, double* reX, double* imX, int noPoles, int noFFTPoints, double delta, double mu, double& eold, double& enew);

	/**
	* @brief   Calculate the frequency domain error
	*
	* general
	* @param [in]	a			Holds the \a a recursion coefficients
	* @param [in]	b			Holds the \a b recursion coefficients
	* @param [in]	t			Holds the desired frequency response
	* @param [in]	reX			Holds the real part of the signal during FFT
	* @param [in]	imX			Holds the imaginary part of the signal during FFT
	* @param [in]	noPoles		Holds the Number of poles in the filter
	* @param [in]	noFFTPoints	Holds the Number of points in FFT
	* @return Double value holding the current frequency domain error
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 26 > Eq 26-5
	*/
	double	CalculateFrequencyDomainError(double* a, double* b, double* t, double* reX, double* imX, int noPoles, int noFFTPoints);

	/**
	* @brief   Combining Recursion Coefficients of Cascade and Parallel Stages
	*
	* Frequency responses of systems in a cascade are combined by multiplication, while in parallel they are combined by addition. These same rules follow z-domain transfer functions, so we can do these operation in z-domain.
	* @param [in]	a1			array for \a a coeff for sys 1, one of the stages
	* @param [in]	a1Size		Its size
	* @param [in]	b1			array for \a b coeff for sys 1, one of the stages
	* @param [in]	b1Size		Its size
	* @param [in]	a2			array for \a a coeff for sys 2, one of the stages
	* @param [in]	a2Size		Its size
	* @param [in]	b2			array for \a b coeff for sys 2, one of the stages
	* @param [in]	b2Size		Its size
	* @param [in]	a3			array for \a a coeff for sys 3, the combined system
	* @param [in]	a3Size		Its size
	* @param [in]	b3			array for \a b coeff for sys 3, the combined system
	* @param [in]	b3Size		Its size
	* @param [in]	bParallel	Indicate parallel or cascade combination
	* @return N/A
	* @sa "The Scientist and Engineer's Guide to Digital Signal Processing" > Ch 31 > Table 31-1
	*/
	void	CombineRecursionCoeffCascadeParallelStages(double* a1, int a1Size, double* b1, int b1Size, double* a2, int a2Size,  double* b2, int b2Size, double** a3, int& a3Size, double** b3, int& b3Size, bool bParallel);

	/**
	* @brief This is a leaky integrator
	*
	* Integrates
	* @param [in]	signal		Holds the signals point
	* @return Double integrated signal
	* @sa Digital Signal Processing Coursera class > 5.8
	*/
	double	LeakyIntegrator(double x);

	/**
	* @brief This is an integrator
	*
	* Integrates
	* @param [in]	signal		Holds the signals array
	* @param [in]	a1Size		Its size
	* @return Double integrated point
	* @sa Digital Signal Processing Coursera class > 5.8
	*/
	double*	Integrator(double* x, int size);

	/**
	* @brief This is a moving average
	* Averages
	* @param [in]	signal		Holds the signals point
	* @return Double averaged signal
	* @sa Digital Signal Processing Coursera class > 5.8
	*/
	double	MovingAvegare(double x);

	/**
	* @brief This is an average
	*
	* Averages
	* @param [in]	signal		Holds the signals array
	* @param [in]	a1Size		Its size
	* @return Double integrated point
	* @sa Digital Signal Processing Coursera class > 5.8
	*/
	double*	Average(double* x, int size);

	/**
	* @brief This is a callback for producer/consumer read/write buffer
	*
	* Callback to process ready data
	* @param [in]	input		Holds the input data buffer
	* @param [in]	output		Holds the output data buffer
	* @param [in]	samples		Number of samples to read
	* @return Double integrated point
	* @sa Digital Signal Processing Coursera class > 5.11
	*/
	void		Callback(const void* input, void* output, unsigned long samples);

	/**
	* @brief This is an Three point echo
	*
	* No natural low pass filtering after a bounce repetition
	* @param [in]	signal			Holds the input data buffer
	* @param [in]	size			Size of the signal buffer
	* @param [in]	sampling_freq	Sampling frequency of the signal
	* @param [in]	index			Index of the data point you want to process
	* @return Float processed point
	* @sa Digital Signal Processing Coursera class > 5.11
	*/
	float		Echo(double* signal, int size, int sampling_freq, int index);

	/**
	* @brief This is an Three point natural echo
	*
	* No natural low pass filtering after a bounce repetition
	* @param [in]	signal			Holds the input data buffer
	* @param [in]	processed_signal	Holds the output data buffer
	* @param [in]	size			Size of the signal buffer
	* @param [in]	sampling_freq	Sampling frequency of the signal
	* @param [in]	index_x			Index of the data point you want to process
	* @param [in]	index_y			Index of the data point processed
	* @return Float processed point
	* @sa Digital Signal Processing Coursera class > 5.11
	*/
	float	NaturalEcho(double* signal, double* processed_signal, int size, int sampling_freq, int index_x, int index_y);

	/**
	* @brief This is an reverb
	*
	* Reverb
	* @param [in]	signal			Holds the input data buffer
	* @param [in]	processed_signal	Holds the output data buffer
	* @param [in]	size			Size of the signal buffer
	* @param [in]	sampling_freq	Sampling frequency of the signal
	* @param [in]	index_x			Index of the data point you want to process
	* @param [in]	index_y			Index of the data point processed
	* @return Float processed point
	* @sa Digital Signal Processing Coursera class > 5.11
	*/
	float	Reverb(double* signal, double* processed_signal, int size, int sampling_freq, int index_x, int index_y);

	/**
	* @brief This is an fuzz
	*
	* Reverb
	* @param [in]	signal			Holds the input data buffer
	* @param [in]	index_x			Index of the data point you want to process
	* @return Float processed point
	* @sa Digital Signal Processing Coursera class > 5.11
	*/
	float	Fuzz(double* signal, int index_x);

	/**
	* @brief This is an tremolo
	*
	* Reverb
	* @param [in]	signal			Holds the input data buffer
	* @param [in]	sampling_freq	Sampling frequency of the signal
	* @param [in]	index_x			Index of the data point you want to process
	* @return Float processed point
	* @sa Digital Signal Processing Coursera class > 5.11
	*/
	float	Tremolo(double* signal, int sampling_freq, int index_x);

	/**
	* @brief This is an flanger
	*
	* Reverb
	* @param [in]	signal			Holds the input data buffer
	* @param [in]	size			Size of the signal buffer
	* @param [in]	sampling_freq	Sampling frequency of the signal
	* @param [in]	index_x			Index of the data point you want to process
	* @return Float processed point
	* @sa Digital Signal Processing Coursera class > 5.11
	*/
	float	Flanger(double* signal, int size, int sampling_freq, int index_x);
};

#endif
