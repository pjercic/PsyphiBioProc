#include <BDFWriter.h>
using namespace std;

BDFWriter::BDFWriter () {
	;
}

void BDFWriter::DisplayMessage(string message) {

	cout << message << endl;

	int hdl, chns;

	chns = signals.size();

	//  Opens an new file for writing. Warning, an already existing file with the same name will be silently overwritten without advance warning!
	hdl = edfopen_file_writeonly("sinus.bdf", EDFLIB_FILETYPE_BDFPLUS, chns);

	if(hdl<0) {
		printf("error: edfopen_file_writeonly()\n");
		return;
	}

	// Sets the samplefrequency of signal edfsignal.
	// This function is required for every signal and can be called only after opening a file in writemode and before the first sample write action
	for(int i=0; i < chns; i++) {
		if(edf_set_samplefrequency(hdl, i, smpFreqs[i])) {
			printf("error: edf_set_samplefrequency()\n");
			return;
		}
	}	

	for(int i=0; i < chns; i++) {

		double phy_max = 0;
		for(int j=0; j<nosSamples[i]; j++)
			if(signals[i][j] > phy_max) phy_max = signals[i][j];

		// Sets the maximum physical value of signal edfsignal.
		// This function is required for every signal and can be called only after opening a file in writemode and before the first sample write action.
		// Sets the minimum physical value of signal edfsignal.
		// Usually this will be (-(phys_max))
		// This function is required for every signal and can be called only after opening a file in writemode and before the first sample write action.
		if(edf_set_physical_maximum(hdl, i, phy_max)) {
			printf("error: edf_set_physical_maximum()\n");
			return;
		}

		if(edf_set_physical_minimum(hdl, i, -phy_max)) {
			printf("error: edf_set_physical_minimum()\n");
			return;
		}

		// Sets the maximum digital value of signal edfsignal. Usually, the value 32767 is used for EDF+ and 8388607 for BDF+.
		// Sets the minimum digital value of signal edfsignal. Usually, the value -32768 is used for EDF+ and -8388608 for BDF+.
		// Usually this will be (-(dig_max + 1))
		// This function is required for every signal and can be called only after opening a file in writemode and before the first sample write action.
		if(edf_set_digital_maximum(hdl, i, 8388607)) {
			printf("error: edf_set_digital_maximum()\n");
			return;
		}

		if(edf_set_digital_minimum(hdl, i, -8388608)) {
			printf("error: edf_set_digital_minimum()\n");
			return;
		}

		// Set labels
		if(edf_set_label(hdl, i, sigLabels[i].c_str())) {
			printf("error: edf_set_label()\n");
			return;
		}

		// Sets the physical dimension of signal edfsignal ("uV", "BPM", "mA", "Degr.", etc.).
		if(edf_set_physical_dimension(hdl, i, sigDim[i].c_str())) {
			printf("error: edf_set_physical_dimension()\n");
			return;
		}
	}


	// Writes n physical samples (uV, mA, Ohm) from *buf belonging to one signal where n is the samplefrequency of the signal.
	// The physical samples will be converted to digital samples using the values of physical maximum, physical minimum, digital maximum and digital minimum.
	// The number of samples written is equal to the samplefrequency of the signal.
	// Size of buf should be equal to or bigger than sizeof(double[samplefrequency]).
	// Call this function for every signal in the file. The order is important!
	// When there are 4 signals in the file, the order of calling this function must be: signal 0, signal 1, signal 2, signal 3, signal 0, signal 1, signal 2, etc.
	
	bool bStop = false;
	for(int i=0; !bStop; i++) {
		bStop = true;

		for(int j=0; j < signals.size(); j++) {

			if(smpFreqs[j] * (i + 1) <= nosSamples[j]) {
				bStop = false;
				if(edfwrite_physical_samples(hdl, signals[j] + smpFreqs[j] * i)) {
					printf("error: edfwrite_physical_samples()\n");
					return;
				}
			}
		}
	}

	// Writes an annotation/event to the file.
	// onset is relative to the starttime and startdate of the file.
	// onset and duration are in units of 100 microSeconds! resolution is 0.0001 second!
	// For example: 34.071 seconds must be written as 340710.
	edfwrite_annotation_latin1(hdl, 0LL, -1LL, "Recording starts");

	//edfwrite_annotation_latin1(hdl, noSamples/smp_freq * 10000LL, -1LL, "Recording ends");

	edfclose_file(hdl);
}

void BDFWriter::saveSignal(double* sig, int noSamp, int smp_freq, string sigLabel, string sigDimension) {
	signals.push_back(sig);
	smpFreqs.push_back(smp_freq);
	nosSamples.push_back(noSamp);
	sigLabels.push_back(sigLabel);
	sigDim.push_back(sigDimension);
}

void BDFWriter::saveNextDataPoint() {
	cout << "Console:saveNextDataPoint" << endl;
}

void BDFWriter::DisplayHistogram(deque<int> histogram) {
	cout << "Console:DisplayHistogram" << endl;
}