#include "DOA.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

vector<double> convid_same(vector<double> &signal, vector<double> &kernel)
{
    vector<double> result(signal.size(), 0.0);
    int kernelCenter = kernel.size() / 2;

    for (int k = 0; k < signal.size(); k++)
    {
        for (int l = 0; l < kernel.size(); l++)
        {
            int signalIdx = k - l + kernelCenter;

            if (signalIdx >= 0 && signalIdx < signal.size())
            {
                result[k] += signal[signalIdx] * kernel[l];
            }
        }
    }
    return result;
}

void DOA_Histogram(const vector<double> &Blockestimates, double windowlength, vector<double> &PlainHist, vector<double> &SmoothedHist)
{

    vector<int> histbinning(360);
    iota(histbinning.begin(), histbinning.end(), 0);


    vector<double> estangle = Blockestimates;

    // Initialization of histogram with 0
    PlainHist.resize(histbinning.size(), 0);

    // Calculation of bin width
    double binWidth = (histbinning.size() > 1) ? (histbinning[1] - histbinning[0]) : 1.0;

    for (int k = 0; k < estangle.size(); k++)
    {
        int binIndex = static_cast<int>(floor(estangle[k] - histbinning[0]) / binWidth);

        if (binIndex >= 0 && binIndex <  static_cast<int>(histbinning.size()))
        {
            PlainHist[binIndex]++;
        }
    }
    // for (int n = 0; n < 15; n++)
    // {
    //     cout << "n:" << n << " " << PlainHist[n] << endl;
    // }

    //Definition of the Gaussian window parameters

    int gaussWidth = 5; // Standard deviation of the gaussian window
    int windowSize = 2 * gaussWidth + 1;

    vector<double> gaussWindow(windowSize);
    for (int i = 0; i < windowSize; i++)
    {
        double x = i - gaussWidth;
        gaussWindow[i] = exp(-x * x / (2*gaussWidth * gaussWidth));
    }

    double sum = accumulate(gaussWindow.begin(), gaussWindow.end(), 0.0);

    for (int j = 0; j < windowSize; j++)
    {
        gaussWindow[j] /= sum;
    }

    SmoothedHist = convid_same(PlainHist, gaussWindow);
}