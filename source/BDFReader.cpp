#include "BDFReader.h"
using namespace std;

BDFReader::BDFReader () {
	sampleCounter = 0;
}

double* BDFReader::getSignal(int& noSamples, int& smpFreq) {

	noSamples = hdr.signalparam[channelNo].smp_in_file;
	smpFreq = samplingFreq;
	return buf;
}

double BDFReader::getNextDataPoint() {

	return buf[sampleCounter++];
}

void BDFReader::Open(string fileName, int channel) {

	if(channel < 0) {
		cout << "signal nr must be positive" << endl;
		return;
	}

	// Open the file for reading
	if(edfopen_file_readonly(fileName.c_str(), &hdr, EDFLIB_READ_ALL_ANNOTATIONS)) {

		// Error reporting
		switch(hdr.filetype) {
			case EDFLIB_MALLOC_ERROR                : printf("\nmalloc error\n\n");
				break;
			case EDFLIB_NO_SUCH_FILE_OR_DIRECTORY   : printf("\ncan not open file, no such file or directory\n\n");
				break;
			case EDFLIB_FILE_CONTAINS_FORMAT_ERRORS : printf("\nthe file is not EDF(+) or BDF(+) compliant\n (it contains format errors)\n\n");
				break;
			case EDFLIB_MAXFILES_REACHED            : printf("\nto many files opened\n\n");
				break;
			case EDFLIB_FILE_READ_ERROR             : printf("\na read error occurred\n\n");
				break;
			case EDFLIB_FILE_ALREADY_OPENED         : printf("\nfile has already been opened\n\n");
				break;
			default                                 : printf("\nunknown error\n\n");
				break;
		}

		return;
	}

	cout << "File: " << fileName << " opened sucesfully." << endl;

	// Check the correct specification of a channel
	if(channel > hdr.edfsignals) {
		printf("\nerror: file has %i signals and you selected signal %i\n\n", hdr.edfsignals, channel);
		edfclose_file(hdr.handle);
		return;
	}

	channelNo = channel;

	cout << "Channel " << channel << ":" << hdr.signalparam[channel].label << " read sucesfully." << endl;
	PrintHeader(channel);
	ReadData(channel);

	time_t rawtime;

  time ( &rawtime );

  datetime = gmtime ( &rawtime );

	datetime->tm_hour = hdr.starttime_hour;
	datetime->tm_mday = hdr.startdate_day;
	datetime->tm_min = hdr.starttime_minute;
	datetime->tm_mon = hdr.startdate_month;
	datetime->tm_sec = hdr.starttime_second;
	datetime->tm_year = hdr.startdate_year;

	cout << asctime (datetime);

	samplingFreq = ((double)hdr.signalparam[channel].smp_in_datarecord / (double)hdr.datarecord_duration) * EDFLIB_TIME_DIMENSION;

	edfclose_file(hdr.handle);
}

void BDFReader::Close() {

	edfclose_file(hdr.handle);
}

void BDFReader::PrintHeader (int nChannel) {

	cout << "-- HEADER -- ****************************" << endl;
	printf("\nlibrary version: %i.%02i\n", edflib_version() / 100, edflib_version() % 100);

	printf("\ngeneral header:\n\n");

	printf("filetype: %i\n", hdr.filetype);
	printf("edfsignals: %i\n", hdr.edfsignals);
#ifdef WIN32
	printf("file duration: %I64d seconds\n", hdr.file_duration / EDFLIB_TIME_DIMENSION);
#else
	printf("file duration: %lli seconds\n", hdr.file_duration / EDFLIB_TIME_DIMENSION);
#endif
	printf("startdate: %i-%i-%i\n", hdr.startdate_day, hdr.startdate_month, hdr.startdate_year);
	printf("starttime: %i:%02i:%02i\n", hdr.starttime_hour, hdr.starttime_minute, hdr.starttime_second);
	printf("patient: %s\n", hdr.patient);
	printf("recording: %s\n", hdr.recording);
	printf("patientcode: %s\n", hdr.patientcode);
	printf("gender: %s\n", hdr.gender);
	printf("birthdate: %s\n", hdr.birthdate);
	printf("patient_name: %s\n", hdr.patient_name);
	printf("patient_additional: %s\n", hdr.patient_additional);
	printf("admincode: %s\n", hdr.admincode);
	printf("technician: %s\n", hdr.technician);
	printf("equipment: %s\n", hdr.equipment);
	printf("recording_additional: %s\n", hdr.recording_additional);
	printf("datarecord duration: %f seconds\n", ((double)hdr.datarecord_duration) / EDFLIB_TIME_DIMENSION);
#ifdef WIN32
	printf("number of datarecords in the file: %I64d\n", hdr.datarecords_in_file);
	printf("number of annotations in the file: %I64d\n", hdr.annotations_in_file);
#else
	printf("number of datarecords in the file: %lli\n", hdr.datarecords_in_file);
	printf("number of annotations in the file: %lli\n", hdr.annotations_in_file);
#endif

	printf("\nsignal parameters:\n\n");

	printf("label: %s\n", hdr.signalparam[nChannel].label);
#ifdef WIN32
	printf("samples in file: %I64d\n", hdr.signalparam[nChannel].smp_in_file);
#else
	printf("samples in file: %lli\n", hdr.signalparam[nChannel].smp_in_file);
#endif
	printf("samples in datarecord: %i\n", hdr.signalparam[nChannel].smp_in_datarecord);
	printf("physical maximum: %f\n", hdr.signalparam[nChannel].phys_max);
	printf("physical minimum: %f\n", hdr.signalparam[nChannel].phys_min);
	printf("digital maximum: %i\n", hdr.signalparam[nChannel].dig_max);
	printf("digital minimum: %i\n", hdr.signalparam[nChannel].dig_min);
	printf("physical dimension: %s\n", hdr.signalparam[nChannel].physdimension);
	printf("prefilter: %s\n", hdr.signalparam[nChannel].prefilter);
	printf("transducer: %s\n", hdr.signalparam[nChannel].transducer);
	printf("samplefrequency: %f\n", ((double)hdr.signalparam[nChannel].smp_in_datarecord / (double)hdr.datarecord_duration) * EDFLIB_TIME_DIMENSION);

	cout << "****************************" << endl;
}

void BDFReader::ReadData(int nChannel) {

	// Allocate memory for the converted sampes read
	buf=new double[hdr.signalparam[nChannel].smp_in_file];
	if(buf==NULL) {
		printf("\allocation error\n");
		edfclose_file(hdr.handle);
		return;
	}

	// Read physical converted samples
	int nSamplesRead;
	nSamplesRead = edfread_physical_samples(hdr.handle, nChannel, hdr.signalparam[nChannel].smp_in_file, buf);

	if(nSamplesRead < 0) {
		printf("\nerror: edf_read_physical_samples()\n");
		edfclose_file(hdr.handle);
		delete [] buf;
		return;
	}
}