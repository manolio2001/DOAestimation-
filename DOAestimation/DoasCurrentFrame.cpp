#include <iostream>
#include <complex>
#include <fftw3.h>
#include <algorithm>
#include <functional>
#include <cmath>
#include <numeric>
#include <fstream>
#include "DOA.hpp"

using namespace std;

double td_adjacent(double l, double phi, int i, double a, double c)
{
    return l * sin(M_PI + a / 2 - phi + (i - 1) * a) / c;
}

vector<double> hamming(int frame_size)
{
    vector<double> window(frame_size);
    for (int i = 0; i < frame_size; i++)
    {
        window[i] = 0.54 - 0.46 * cos(2 * M_PI * i / (frame_size - 1));
    }
    return window;
}

vector<complex<double>> FFT(vector<double> &data, int fftsize)
{
    vector<complex<double>> result(fftsize);
    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize);
    
    if (!in || !out) {
        if (in) fftw_free(in);
        if (out) fftw_free(out);
        throw runtime_error("Failed to allocate FFTW memory");
    }
    
    fftw_plan p = fftw_plan_dft_1d(fftsize, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    for (int i = 0; i < data.size(); i++)
    {
        if (i < data.size())
        {
            in[i][0] = data[i];
            in[i][1] = 0.0;
        }
        else
        {
            in[i][0] = 0.0;
            in[i][1] = 0.0;
        }
    }
    fftw_execute(p);
    for (int l = 0; l < fftsize; l++)
    {
        result[l] = complex<double>(out[l][0], out[l][1]);
    }
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return result;
}

vector<double> Karabasi_estimation_AllbestBinPairs(
    int best_zone, int NoEst, vector<vector<double>> &current_data,
    int frame_size, int fftsize, double zone_step, int K, int M,
    int ignore_bins, double c, int Fs, double radius)
{
    cout << "Karabasi_estimation_AllbestBinPairs called with NoEst: " << NoEst << endl;
    double a = 2 * M_PI / M;               // Angle between microphones
    double l = 2 * radius * sin(M_PI / M); // Distance between adjacent microphones
    vector<double> estimates;
    vector<vector<complex<double>>> FreqCorr(M);
    vector<vector<double>> FreqCorr_abs(M, vector<double>(fftsize, 0.0));
    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize);
    
    if (!in || !out) {
        if (in) fftw_free(in);
        if (out) fftw_free(out);
        throw runtime_error("Failed to allocate FFTW memory");
    }
    
    fftw_plan p = fftw_plan_dft_1d(fftsize, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    
    for (int i = 0; i < M; i++)
    {
        FreqCorr[i].resize(fftsize);
        vector<double> spart1(current_data.size());
        vector<double> spart2(current_data.size());
        for (int j = 0; j < current_data.size(); j++)
        {
            spart1[j] = current_data[j][i];

            if (i + 1 == M)
            {
                spart2[j] = current_data[j][0];
            }
            else
            {
                spart2[j] = current_data[j][i + 1];
            }
        }

        vector<double> s_wndowed1(frame_size);
        vector<double> s_wndowed2(frame_size);
        vector<double> window_data = hamming(frame_size);
        for (int m = 0; m < frame_size; m++)
        {
            s_wndowed1[m] = spart1[m] * window_data[m];
            s_wndowed2[m] = spart2[m] * window_data[m];
        }

        int beg = ignore_bins + 1 + zone_step * (best_zone - 1);
        int fin = ignore_bins + zone_step * (best_zone - 1) + K;

        vector<complex<double>> X = FFT(s_wndowed1, fftsize);
        vector<complex<double>> Y = FFT(s_wndowed2, fftsize);

        vector<complex<double>> Xcomponent(X.begin() + (beg - 1), X.begin() + fin);
        vector<complex<double>> Ycomponent(Y.begin() + (beg - 1), Y.begin() + fin);

        vector<complex<double>> Xhalfull(fftsize);
        vector<complex<double>> Yhalfull(fftsize);

        if (beg >= 1 && beg - 1 + Xcomponent.size() <= fftsize)
        {
            copy(Xcomponent.begin(), Xcomponent.end(), Xhalfull.begin() + (beg - 1));
            copy(Ycomponent.begin(), Ycomponent.end(), Yhalfull.begin() + (beg - 1));
        }

        for (int l = 0; l < fftsize; l++)
        {
            FreqCorr[i][l] = Xhalfull[l] * conj(Yhalfull[l]);
        }
    }

    for (int p = 0; p < M; p++)
    {
        for (int l = 0; l < fftsize; l++)
        {
            FreqCorr_abs[p][l] = abs(FreqCorr[p][l]);
        }
    }

    vector<int> max_frq_cols;
    double max_val = 0.0;
    int max_col = 0, max_row = 0;
    int cols = FreqCorr_abs[0].size();
    for (int n = 0; n < NoEst; n++)
    {
        for (int i = 0; i < M; i++)
        {
            max_val = 0.0;

            for (int j = 0; j < cols; j++)
            {
                if (FreqCorr_abs[i][j] > max_val)
                {
                    max_val = FreqCorr_abs[i][j];
                    max_row = i;
                    max_col = j;
                }
            }
        }
        max_frq_cols.push_back(max_col);
        FreqCorr_abs[max_row][max_col] = 0.0;
    }

    vector<int> J(NoEst, -1);

    for (int i = 0; i < NoEst; i++)
    {
        J[i] = max_frq_cols[i];
    }

    vector<double> tangle;
    for (double i = 0; i <= 359; i++)
    {
        tangle.push_back(i * 2 * M_PI / 360);
    }

    for (int n = 0; n < NoEst; n++)
    {
        int Wmax_all = J[n];
        cout << "Processing estimate " << n + 1 << " of " << NoEst << " with Wmax_all: " << Wmax_all << endl;
        vector<complex<double>> G(M);
        vector<vector<complex<double>>> PRF(M, vector<complex<double>>(tangle.size()));
        vector<vector<complex<double>>> CICS(M, vector<complex<double>>(tangle.size()));
        for (int pair_count = 0; pair_count < M; pair_count++)
        {
            G[pair_count] = FreqCorr[pair_count][Wmax_all] / abs(FreqCorr[pair_count][Wmax_all]);
            vector<double> td(tangle.size(), 0.0);

            for (int w = 0; w < tangle.size(); w++)
            {
                double td_1 = td_adjacent(l, tangle[w], 1, a, c) * Fs;
                double td_2 = td_adjacent(l, tangle[w], pair_count + 1, a, c) * Fs;
                td[w] = td_1 - td_2;
                PRF[pair_count][w] = exp(complex<double>(0, -2 * M_PI * (Wmax_all + 1) * td[w] / fftsize));
                CICS[pair_count][w] = G[pair_count] * PRF[pair_count][w];
            }
        }

        vector<double> CICS_sum_abs(tangle.size(), 0.0);
        for (int i = 0; i < tangle.size(); i++)
        {
            complex<double> sum(0.0, 0.0);
            for (int j = 0; j < M; j++)
            {
                sum += CICS[j][i];
            }
            CICS_sum_abs[i] = abs(sum);
        }

        auto maxCICS = max_element(CICS_sum_abs.begin(), CICS_sum_abs.end());
        double CICS_angle_index = distance(CICS_sum_abs.begin(), maxCICS);
        double CICS_angle = CICS_angle_index - 1.0;
        double estangle = CICS_angle;
        if (estangle < 0)
        {
            estangle += 360.0;
        }
        estimates.push_back(estangle);
        cout << "Added estimate: " << estangle << endl;
    }
    cout << "Karabasi_estimation_AllbestBinPairs returning " << estimates.size() << " estimates" << endl;
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    return estimates;
}

vector<double> DOAsCurrentFrame(vector<vector<double>> current_data, int fft_size, int ignore_bins, int Nk, double MMC, double zone_step, int K, double c, double MyFs, int NopEst, double radius)
{
    int frame_size = current_data.size();        // Number of time frames (rows)
    int M = current_data[0].size();              // Number of channels (columns)
    double a = 2 * M_PI / M;                     // Angle between microphones
    double l = 2 * radius * sin(M_PI / M);       // Distance between adjacent microphones
    vector<double> window = hamming(frame_size); // Hamming window
    vector<vector<double>> windowed_data(frame_size, vector<double>(M));
    for (int i = 0; i < frame_size; i++)
    {
        for (int j = 0; j < M; j++)
        {
            windowed_data[i][j] = current_data[i][j] * window[i]; // Multiply each channel by the window
        }
    }

    //  For the FFT
    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);
    fftw_plan p = fftw_plan_dft_1d(fft_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Perform FFT on each channel of the windowed data
    vector<vector<complex<double>>> FFT_data(fft_size, vector<complex<double>>(M));
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < fft_size; j++)
        {
            if (j < windowed_data.size()) // Check time index is within bounds
            {
                in[j][0] = windowed_data[j][i]; // Real part
                in[j][1] = 0.0;                 // Imaginary part
            }
            else
            {
                in[j][0] = 0.0; // Zero padding
                in[j][1] = 0.0;
            }
        }
        //  Execute FFT
        fftw_execute(p);

        //  Copy FFt results to FFT_data
        for (int k = 0; k < fft_size / 2 + 1; k++)
        {
            FFT_data[k][i] = complex<double>(out[k][0], out[k][1]);
        }
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    int pair = 0;
    vector<vector<double>> rcoef(Nk - 1, vector<double>(M - 1, 0.0));
    vector<vector<complex<double>>> X_complete;
    vector<vector<complex<double>>> Y_complete;

    for (int i = 0; i < M; i++)
    {
        vector<complex<double>> X(fft_size / 2 + 1);
        vector<complex<double>> Y(fft_size / 2 + 1);

        for (int j = i + 1; j <= M; j++)
        {
            int current_j = (j == M) ? 0 : j;

            for (int l = 0; l < fft_size / 2 + 1; l++)
            {
                if (current_j == 0)
                {
                    X[l] = FFT_data[l][0]; // Column i (microphone i)
                    Y[l] = FFT_data[l][1]; // Column j (microphone j)
                }
                else
                {
                    X[l] = FFT_data[l][i];         // Column i (microphone i)
                    Y[l] = FFT_data[l][current_j]; // Column j (microphone j)
                }
            }
            if (current_j == 0)
            {
                break;
            }
            else
            {
                for (int k = 0; k < Nk - 1; k++)
                {
                    int start = ignore_bins + zone_step * k;
                    int end = start + K;

                    if (end > FFT_data.size())
                        break;

                    double Rxy = 0.0, Rxx = 0.0, Ryy = 0.0;
                    for (int n = start; n < end; n++)
                    {
                        Rxy += abs(X[n]) * abs(Y[n]);
                        Rxx += abs(X[n]) * abs(X[n]);
                        Ryy += abs(Y[n]) * abs(Y[n]);
                    }

                    Rxy = Rxy / static_cast<double>(K);
                    Rxx = Rxx / static_cast<double>(K);
                    Ryy = Ryy / static_cast<double>(K);

                    if (sqrt(Rxx * Ryy) != 0)
                    {
                        rcoef[k][pair] = Rxy / sqrt(Rxx * Ryy);
                    }
                    else
                    {
                        rcoef[k][pair] = 0.0;
                    }
                }
            }
        }
        pair++;
    }
    vector<double> MCCs;
    for (int i = 0; i < rcoef.size(); i++)
    {
        MCCs.push_back(accumulate(rcoef[i].begin(), rcoef[i].end(), 0.0) / rcoef[i].size());
    }

    vector<vector<double>> sorted_zones_descend = locate_zones(MCCs, MCCs.size(), MMC);
    reverse(sorted_zones_descend.begin(), sorted_zones_descend.end());
    int nb = sorted_zones_descend.size();
    cout << "DOAsCurrentFrame: Found " << nb << " zones to process." << endl;
    for (int z = 0; z < nb; z++) {
        cout << "  Zone index: " << z << ", best_zone: " << sorted_zones_descend[z][1] << endl;
    }
    vector<double> estimates;

    if (nb != 0)
    {
        for (int z = 0; z < nb; z++)
        {
            int best_zone = sorted_zones_descend[z][1];
            vector<double> zone_estimates = Karabasi_estimation_AllbestBinPairs(best_zone, NopEst, current_data, frame_size, fft_size, zone_step, K, M, ignore_bins, c, MyFs, radius);
            estimates.insert(estimates.end(), zone_estimates.begin(), zone_estimates.end());
        }
    }
    cout << "DOAsCurrentFrame: Returning " << estimates.size() << " estimates." << endl;
    return estimates;
}