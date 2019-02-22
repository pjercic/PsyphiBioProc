#include <iostream>
#include <gtest.h>
#include <IInput.h>
#include <BDFReader.h>
#include <DSPStatistic.h>
using namespace std;

// Tests add function
TEST(BDFReader, ComapareSignalData) {
	IInput* bdfReader;
	IInput* bdfReader_a;

	bdfReader = new BDFReader();
	bdfReader_a = new BDFReader();

	bdfReader->Open("test.bdf", 0);
	bdfReader_a->Open("test.bdf", 0);
	
	int noSamples = 0;
	int sampFreq = 0;
	double* signal = bdfReader->getSignal(noSamples, sampFreq);

	double signalDataPoint = 0;
	signalDataPoint = bdfReader->getNextDataPoint();
	EXPECT_EQ(signal[0], signalDataPoint) << "both numbers a and b are the same";

	signalDataPoint = bdfReader->getNextDataPoint();
	EXPECT_EQ(signal[1], signalDataPoint) << "both numbers a and b are the same";

	signalDataPoint = bdfReader->getNextDataPoint();
	EXPECT_EQ(signal[2], signalDataPoint) << "both numbers a and b are the same";

	bdfReader_a->Close();
	bdfReader->Close();

	delete bdfReader_a;
	delete bdfReader;
}

namespace {

	// The fixture for testing class Foo.
	class InputTest : public ::testing::Test {

	protected:

		IInput* bdfReader;
		IInput* bdfReader_a;

		// You can remove any or all of the following functions if its body
		// is empty.

		InputTest() {
			// You can do set-up work for each test here.
			bdfReader = new BDFReader();
			bdfReader_a = new BDFReader();
		}

		virtual ~InputTest() {
			// You can do clean-up work that doesn't throw exceptions here.
			delete bdfReader_a;
			delete bdfReader;
		}

		// If the constructor and destructor are not enough for setting up
		// and cleaning up each test, you can define the following methods:

		virtual void SetUp() {
			// Code here will be called immediately after the constructor (right
			// before each test).
			bdfReader->Open("test.bdf", 0);
			bdfReader_a->Open("test.bdf", 0);
		}

		virtual void TearDown() {
			// Code here will be called immediately after each test (right
			// before the destructor).
			bdfReader_a->Close();
			bdfReader->Close();
		}

		// Objects declared here can be used by all tests in the test case for Foo.
	};

	// Tests that the Foo::Bar() method does Abc.
	TEST_F(InputTest, GetSignal) {
		int noSamples = 0;
		int sampFreq = 0;
		double* signal = bdfReader->getSignal(noSamples, sampFreq);

		ASSERT_EQ(20480, noSamples) << "Number of samples in bdf test file is ";
		ASSERT_EQ(4495.0229443325370, signal[0]) << "both numbers a and b are the same";
		ASSERT_EQ(4500.3979344009122, signal[1]) << "both numbers a and b are the same";
		ASSERT_EQ(4497.1791903483381, signal[2]) << "both numbers a and b are the same";
	}

	// Tests that Foo does Xyz.
	TEST_F(InputTest, GetNextDataPoint) {
		// Exercises the Xyz feature of Foo.
		double signalDataPoint = 0;
		signalDataPoint = bdfReader->getNextDataPoint();
		EXPECT_EQ(4495.0229443325370, signalDataPoint) << "both numbers a and b are the same";

		signalDataPoint = bdfReader->getNextDataPoint();
		EXPECT_EQ(4500.3979344009122, signalDataPoint) << "both numbers a and b are the same";

		signalDataPoint = bdfReader->getNextDataPoint();
		EXPECT_EQ(4497.1791903483381, signalDataPoint) << "both numbers a and b are the same";
	}

}  // namespace

namespace {

	// The fixture for testing descriptives
	class DescriptivesTest : public ::testing::Test {

	protected:

		DSPStatistic* dspStatistic;
		double* signal;

		// You can remove any or all of the following functions if its body
		// is empty.

		DescriptivesTest() {
			// You can do set-up work for each test here.
			dspStatistic = new DSPStatistic();
			
		}

		virtual ~DescriptivesTest() {
			// You can do clean-up work that doesn't throw exceptions here.
			delete dspStatistic;
		}

		// If the constructor and destructor are not enough for setting up
		// and cleaning up each test, you can define the following methods:

		virtual void SetUp() {
			// Code here will be called immediately after the constructor (right
			// before each test).
			;
		}

		virtual void TearDown() {
			// Code here will be called immediately after each test (right
			// before the destructor).
			;
		}

		// Objects declared here can be used by all tests in the test case for Foo.
	};

	// Tests that the Foo::Bar() method does Abc.
	TEST_F(DescriptivesTest, DescriptivesStatic) {

		signal = new double[1000];

		for (int i = 0; i < 1000; i++)
				signal[i] = i + 1;

		double avg, stdev, var, snr, cv;
		
		dspStatistic->CalculateDescriptivesStatic(signal, 1000, avg, stdev, var, snr, cv);

		ASSERT_EQ(500.5, avg) << "Number of samples in bdf test file is ";
		ASSERT_EQ(288.81943609574938, stdev) << "both numbers a and b are the same";
		ASSERT_EQ(83416.666666666672, var) << "both numbers a and b are the same";
		ASSERT_EQ(1.7329166165744965, snr) << "both numbers a and b are the same";
		ASSERT_EQ(57.706181038111758, cv) << "both numbers a and b are the same";

		delete signal;
	}



	// Tests that the Foo::Bar() method does Abc.
	TEST_F(DescriptivesTest, DescriptivesRunning) {

		signal = new double[1000];

		for (int i = 0; i < 1000; i++)
				signal[i] = i + 1;

		double avg, stdev, var, snr, cv;
		
		for (int i = 0; i < 1000; i++)
			dspStatistic->CalculateDescriptivesRunning(signal[i], avg, stdev, var, snr, cv);

		ASSERT_EQ(500.5, avg) << "Number of samples in bdf test file is ";
		ASSERT_EQ(288.81943609574938, stdev) << "both numbers a and b are the same";
		ASSERT_EQ(83416.666666666657, var) << "both numbers a and b are the same";
		ASSERT_EQ(1.7329166165744965, snr) << "both numbers a and b are the same";
		ASSERT_EQ(57.706181038111758, cv) << "both numbers a and b are the same";

		delete signal;
	}

	// Tests that the Foo::Bar() method does Abc.
	TEST_F(DescriptivesTest, DiscreteFourierTransformationTime) {

		double PI = 4 * atan(1.0);

		signal = new double[16];

		for (int i = 0; i < 16; i++) {
			double a = sin((2 * PI) * ((double)i / 8));
			signal[i] = sin((2 * PI) * ((double)i / 8));
		}

		double* realX, *imaginaryX;
		int reXSize, imXSize;
		double* reconstructedTime;
		int reconstructedTimeSize;
		
			dspStatistic->DiscreteFourierTransformationTimeSide(signal, 16, &realX, reXSize, &imaginaryX, imXSize);

			dspStatistic->InverseDiscreteFourierTransformationTimeSide(realX, reXSize, imaginaryX, imXSize, &reconstructedTime, reconstructedTimeSize);

		ASSERT_EQ(signal[1], reconstructedTime[1]) << "Number of samples in bdf test file is ";
	}

	// Tests that the Foo::Bar() method does Abc.
	TEST_F(DescriptivesTest, DiscreteFourierTransformationFreq) {

		double PI = 4 * atan(1.0);

		signal = new double[16];

		for (int i = 0; i < 16; i++) {
			double a = sin((2 * PI) * ((double)i / 8));
			signal[i] = sin((2 * PI) * ((double)i / 8));
		}

		double* realX, *imaginaryX;
		int reXSize, imXSize;
		double* reconstructedTime;
		int reconstructedTimeSize;
		
			dspStatistic->DiscreteFourierTransformationFreqSide(signal, 16, &realX, reXSize, &imaginaryX, imXSize);

			dspStatistic->InverseDiscreteFourierTransformationFreqSide(realX, reXSize, imaginaryX, imXSize, &reconstructedTime, reconstructedTimeSize);

		ASSERT_EQ(signal[1], reconstructedTime[1]) << "Number of samples in bdf test file is ";
	}

	// Tests that the Foo::Bar() method does Abc.
	TEST_F(DescriptivesTest, FastFourierTransformationInteger) {

		double PI = 4 * atan(1.0);

		signal = new double[16];

		for (int i = 0; i < 16; i++) {
			double a = sin((2 * PI) * ((double)i / 8));
			signal[i] = sin((2 * PI) * ((double)i / 8));
		}

		double* realX, *imaginaryX;
		int reXSize, imXSize;
		double* reconstructedTime;
		int reconstructedTimeSize;
		
		dspStatistic->FastFourierTransformReal(signal, 16, &realX, &imaginaryX);
		reconstructedTime = dspStatistic->InverseFastFourierTransformReal(16, realX, imaginaryX);

		ASSERT_EQ(signal[1], reconstructedTime[1]) << "Number of samples in bdf test file is ";
	}
	
	// Tests that the Foo::Bar() method does Abc.
	TEST_F(DescriptivesTest, IntegratorTest) {

		double* signal = new double[20];
		double* integrated_signal;
		signal[0] = 1;
		for(int i = 1; i < 20; i++)
			signal[i] = 0;

		integrated_signal = dspStatistic->Integrator(signal, 20);

		// Should be a linearly decaying signal
		ASSERT_EQ(signal[0], integrated_signal[0]) << "Number of samples in bdf test file is ";
	}

}  // namespace