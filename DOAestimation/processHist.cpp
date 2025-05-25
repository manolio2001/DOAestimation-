/**
 * processHist.cpp
 * ----------------
 * Implements processHist function: given a histogram of angle counts,
 * subtracts Gaussian-weighted contributions iteratively to estimate
 * NSources peak angles.
 *
 * Functions:
 * - linspace: generate a linearly spaced vector.
 * - circshift: circularly shift a vector.
 * - processHist: estimate source angles from histogram.
 */

#include "DOA.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iterator>

using std::vector;
using std::exp;
using std::max_element;
using std::distance;

/**
 * Generate a linearly spaced vector.
 * @param start  Starting value of the sequence.
 * @param end    Ending value of the sequence.
 * @param num    Number of elements to generate.
 * @return       Vector of length num with values from start to end.
 */
vector<double> linspace(double start, double end, int num) {
    vector<double> result;
    if (num <= 0) {
        // No elements to generate
        return result;
    }
    if (num == 1) {
        // Single element equals start
        result.push_back(start);
        return result;
    }
    // Compute uniform step size
    double step = (end - start) / (num - 1);
    result.reserve(num);
    for (int i = 0; i < num; ++i) {
        result.push_back(start + step * i);
    }
    return result;
}

/**
 * Circularly shift a vector by 'shift' positions.
 * Positive shift moves elements to the right; negative shifts to the left.
 * @param v      Input vector to shift.
 * @param shift  Number of positions to shift.
 * @return       New vector after circular shift.
 */
vector<double> circshift(const vector<double>& v, int shift) {
    int n = static_cast<int>(v.size());
    vector<double> result(n);
    // Normalize shift to range [0, n)
    shift = ((shift % n) + n) % n;
    for (int i = 0; i < n; ++i) {
        int newIndex = (i + shift) % n;
        result[newIndex] = v[i];
    }
    return result;
}

/**
 * Estimate NSources angles from a histogram by iterative Gaussian subtraction.
 * @param Histogram   Input histogram (e.g., angle counts per degree).
 * @param NSources    Number of source peaks to estimate.
 * @param gaussWidth  Standard deviation of Gaussian window (in histogram bins).
 * @return            Vector of estimated angles (in degrees).
 */
vector<double> processHist(vector<double>& Histogram, int NSources, double gaussWidth) {
    int H_length = static_cast<int>(Histogram.size());
    // Reserve space for NSources angle estimates
    vector<double> Estimates;
    Estimates.reserve(NSources);

    // Prepare angle axis [0, 1, ..., 359]
    vector<double> angleAxis(360);
    for (int i = 0; i < 360; ++i) {
        angleAxis[i] = i;
    }

    // Create symmetric x-axis range centered at 0 for Gaussian
    int half = H_length / 2;
    double start = -static_cast<double>(half);
    double end   =  static_cast<double>(half);
    vector<double> x = linspace(start, end, H_length);

    // Compute Gaussian window values using gaussWidth
    vector<double> gaussWindow(H_length);
    for (int i = 0; i < H_length; ++i) {
        // Gaussian formula: exp(-x^2 / (2 * sigma^2))
        gaussWindow[i] = exp(-(x[i] * x[i]) / (2 * gaussWidth * gaussWidth));
    }

    // Shift Gaussian so it's centered at the middle of histogram indices
    vector<double> Gauss_shifted = circshift(gaussWindow, -half);

    // Copy histogram to working residual
    vector<double> yn = Histogram;

    // Iteratively find and subtract NSources peaks
    for (int iS = 0; iS < NSources; ++iS) {
        // Find the peak index in current residual
        auto max_it = max_element(yn.begin(), yn.end());
        int cur_max_index = static_cast<int>(distance(yn.begin(), max_it));

        // Map peak index to angle and store estimate
        double estimateAngle = angleAxis[cur_max_index % angleAxis.size()];
        Estimates.push_back(estimateAngle);

        // Generate Gaussian contribution at the peak location
        vector<double> Gauss_contr = circshift(Gauss_shifted, cur_max_index);

        // Subtract weighted Gaussian from histogram and clip negatives
        vector<double> yn_new(H_length);
        for (int i = 0; i < H_length; ++i) {
            // Weighted subtraction of Gaussian
            yn_new[i] = yn[i] - (Gauss_contr[i] * yn[i]);
            // Ensure no negative values remain
            if (yn_new[i] < 0) {
                yn_new[i] = 0;
            }
        }
        // Update residual for next iteration
        yn = yn_new;
    }

    return Estimates;
}
