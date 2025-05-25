#include "DOA.hpp"
#include <numeric>
#include <fstream>

int framesize = 2048;
double time_step = framesize * 0.5;
int fftsize = 2048;
int Fs = 44100;
int UpCutOff = 4000;
int freqzonewidth = 344;
int freqzonebins = round(freqzonewidth * ((double)fftsize / Fs));
double zoneoverlap = 0.5;
double zone_step = freqzonebins * zoneoverlap;
int numzones = floor(UpCutOff / (freqzonewidth * zoneoverlap));
double MMCC = 0.8;
int c = 343;
int NumEstPerZone = 2;
int NSources = 2;
double Ra = 0.05;
int ignore_bins = 5;
int BlocksizeSeconds = 1; // in seconds
int BlocksizeFrames = floor((BlocksizeSeconds * Fs) / time_step);
int blockcounter = 0;
int gaussWidth=40;
vector<double> Blockestimates;
vector<int> NoEstFr(BlocksizeFrames + 2, 0);  // Αρχικοποίηση με BlocksizeFrames + 2 στοιχεία

// Function to read a MATLAB matrix from a .mat file and store it as a 2D vector
vector<vector<double>> ReadMatfile(string filen, string varname)
{
    cout << "Starting to read mat file: " << filen << endl;
    vector<vector<double>> data; // Initialization of a 2D vector to store the data

    mat_t *matfp;     // Pointer to the MAT file
    matvar_t *matvar; // Pointer to the MATLAB variable

    // Opening MAT file in read-only mode
    cout << "Opening mat file..." << endl;
    matfp = Mat_Open(filen.c_str(), MAT_ACC_RDONLY);
    if (matfp == NULL)
    {
        cerr << "Error opening file " << filen << endl; // Error message if the file cannot be opened
        return data;                                    // Return the empty data vector
    }

    // Reading the specified  variable from the MAT file
    cout << "Reading variable: " << varname << endl;
    matvar = Mat_VarRead(matfp, varname.c_str());
    if (matvar == NULL)
    {
        cerr << "Error reading variable: " << varname << " from " << filen << endl; // Error message if the variable cannot be read
        Mat_Close(matfp);                                                         // Close the MAT file
        return data;                                                              // Return the empty vector
    }

    // Check if the variable is a 2D matrix
    cout << "Checking matrix dimensions..." << endl;
    if (matvar->rank == 2 && matvar)
    {
        size_t rows = matvar->dims[0];             // Get the number of rows in the matrix
        size_t cols = matvar->dims[1];             // Get the number of columns in the matrix
        cout << "Matrix dimensions: " << rows << "x" << cols << endl;
        double *mat_data = (double *)matvar->data; // Get the data from the matrix and cast it to a double pointer

        // Reshaping the 2D vector to match the matrix dimensions
        data.resize(rows, (vector<double>(cols)));

        // Fill the 2D vector with the matrix data
        cout << "Filling data vector..." << endl;
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < cols; j++)
            {
                data[i][j] = mat_data[i + j * rows]; // Accessing the data in column-major order
            }
        }
    }

    // Free the memory allocated for the MATLAB variable
    cout << "Freeing MATLAB variable..." << endl;
    Mat_VarFree(matvar);

    // Closing MAT file
    cout << "Closing mat file..." << endl;
    Mat_Close(matfp);

    cout << "Finished reading mat file" << endl;
    return data; // Return the populated 2D vector
}

int main(){
    cout << "Starting main..." << endl;
    string filename = "angle_45.mat";
    string varname = "AuData";

    srand(time(NULL));

    cout << "Reading mat file..." << endl;
    vector<vector<double>> AuData = ReadMatfile(filename, varname);
    if (AuData.empty())
    {
        cout << "Error reading file" << endl;
        return 1;
    }

    cout << "Data read successfully. Size: " << AuData.size() << "x" << (AuData.empty() ? 0 : AuData[0].size()) << endl;

    size_t Lcommon = AuData.size();
    size_t Nsteps = floor((Lcommon - framesize) / time_step + 1);
    cout << "Number of steps: " << Nsteps << endl;

    vector<int> hist_bins(359);
    iota(hist_bins.begin(), hist_bins.end(), 1);

    vector<vector<double>> EstimatesFrames(Nsteps, vector<double>(numzones * NumEstPerZone, NAN));
    vector<vector<double>> Estimated_DOAs(Nsteps, vector<double>(NSources, NAN));

    int windowlength = 5;

    cout << "Starting main processing loop..." << endl;
    for (int istep = 0; istep < Nsteps; istep++){
        size_t start = istep * time_step;
        size_t end = min(start + framesize, AuData.size());

        vector<vector<double>> current_data;
        for(int i = start; i < end; i++){
            current_data.push_back(AuData[i]);
        }

        vector<double> current_estimates = DOAsCurrentFrame(current_data, fftsize, ignore_bins, numzones, MMCC, zone_step, freqzonebins, c, Fs, NumEstPerZone, Ra);
        cout << "current_estimates.size(): " << current_estimates.size() << endl;
        cout << "numzones: " << numzones << ", NumEstPerZone: " << NumEstPerZone << ", expected size: " << numzones * NumEstPerZone << endl;
        copy(current_estimates.begin(), current_estimates.end(), EstimatesFrames[istep].begin());

        if(blockcounter < BlocksizeFrames){
            NoEstFr[blockcounter] = current_estimates.size();
            cout << "NoEstFr[" << blockcounter << "]: " << NoEstFr[blockcounter] << endl;
            Blockestimates.insert(Blockestimates.end(), current_estimates.begin(), current_estimates.end());
            if(blockcounter == 0) {
                cout << "First iteration Blockestimates:" << endl;
                for(size_t i = 0; i < Blockestimates.size(); i++) {
                    cout << Blockestimates[i] << endl;
                }
            }
            blockcounter++;
        }else{
            NoEstFr[blockcounter] = current_estimates.size();
            cout << "NoEstFr[" << blockcounter << "]: " << NoEstFr[blockcounter] << endl;
            Blockestimates.insert(Blockestimates.end(), current_estimates.begin(), current_estimates.end());
        }

        vector<double> PlainHist, SmoothedHist;
        DOA_Histogram(Blockestimates, windowlength, PlainHist, SmoothedHist);
        vector<double> block_DOA = processHist(SmoothedHist, NSources, gaussWidth);
        Estimated_DOAs[istep] = block_DOA;
        
        // Print the current estimates
        cout << "\nFrame " << istep << " DOA Estimates:" << endl;
        for(size_t i = 0; i < block_DOA.size(); i++) {
            cout << "Source " << i+1 << ": " << block_DOA[i] << " degrees" << endl;
        }

        if (!NoEstFr.empty()) {
            int removeCount = NoEstFr.front();
            if (removeCount <= Blockestimates.size()) {
                Blockestimates.erase(Blockestimates.begin(), Blockestimates.begin() + removeCount);
            }
            NoEstFr.erase(NoEstFr.begin());
        } else {
            NoEstFr.resize(BlocksizeFrames + 2, 0);
            blockcounter = 0;
        }

        // Print Estimated_DOAs for debugging
        cout << "\nEstimated_DOAs at step " << istep << ":" << endl;
        for(size_t i = 0; i < Estimated_DOAs[istep].size(); i++) {
            cout << "Source " << i << ": " << Estimated_DOAs[istep][i] << endl;
        }
    }

    // Print final Estimated_DOAs
    cout << "\nFinal Estimated_DOAs:" << endl;
    for(size_t i = 0; i < Estimated_DOAs.size(); i++) {
        cout << "Frame " << i << ": ";
        for(size_t j = 0; j < Estimated_DOAs[i].size(); j++) {
            cout << Estimated_DOAs[i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}
