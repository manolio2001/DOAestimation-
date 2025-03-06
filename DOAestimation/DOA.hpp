#include <iostream>
#include <vector>
#include "math.h"
#include <matio.h>
#include <complex>

using namespace std;

extern int framesize;
extern double time_step;
extern int fftsize;
extern int Fs;
extern int UpCutOff;
extern int freqzonewidth;
extern int freqzonebins;
extern double zoneoverlap;
extern double zone_step;
extern int numzones;
extern double MMCC;
extern int c;
extern int NumEstPerZone;
extern double Ra;
extern int ignore_bins;

// Function to compute DOAs for the current frame
// Parameters:
// - current_data: Input data for the current frame (time domain or frequency domain data)
// - fft_size: FFT size used in the processing
// - ignore_bins: Number of bins to ignore during processing
// - Nk: Number of frequency zones or bins
// - MMC: Magnitude-squared coherence threshold
// - zone_step: Step size between frequency zones (in bins)
// - K: Number of zones
// - c: Speed of sound in the medium
// - MyFs: Sampling frequency
// - NopEst: Number of DOA estimates per zone
// - radius: Radius of the microphone array
vector<double> DOAsCurrentFrame(vector<vector<double>> current_data, int fft_size, int ignore_bins, int Nk, double MMC, double zone_step, int K, double c, double MyFs, int NopEst, double radius);

// Function to locate zones based on a threshold
// Parameters:
// - r: 2D vector input representing the zones
// - threshold: Threshold value for identifying zones
vector<vector<double>> locate_zones(vector<double> &r, int size, double threshold);

void DOA_Histogram(const vector<double> &Blockestimates,double windowlength,vector<double> &PlainHist,vector<double> &SmoothedHist);