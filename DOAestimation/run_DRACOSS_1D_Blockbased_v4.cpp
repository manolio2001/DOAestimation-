#include "DOA.hpp"
#include <numeric>

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
    vector<vector<double>> data; // Initialization of a 2D vector to store the data

    mat_t *matfp;     // Pointer to the MAT file
    matvar_t *matvar; // Pointer to the MATLAB variable

    // Opening MAT file in read-only mode
    matfp = Mat_Open(filen.c_str(), MAT_ACC_RDONLY);
    if (matfp == NULL)
    {
        cerr << "Error opening file " << filen << endl; // Error message if the file cannot be opened
        return data;                                    // Return the empty data vector
    }

    // Reading the specified  variable from the MAT file
    matvar = Mat_VarRead(matfp, varname.c_str());
    if (matvar == NULL)
    {

        cerr << "Error reading variable: " << varname << "from" << filen << endl; // Error message if the variable cannot be read
        Mat_Close(matfp);                                                         // Close the MAT file
        return data;                                                              // Return the empty vector
    }

    // Check if the variable is a 2D matrix
    if (matvar->rank == 2 && matvar)
    {
        size_t rows = matvar->dims[0];             // Get the number of rows in the matrix
        size_t cols = matvar->dims[1];             // Get the number of columns in the matrix
        double *mat_data = (double *)matvar->data; // Get the data from the matrix and cast it to a double pointer

        // Reshaping the 2D vector to match the matrix dimensions
        data.resize(rows, (vector<double>(cols)));

        // Fill the 2D vector with the matrix data
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < cols; j++)
            {
                data[i][j] = mat_data[i + j * rows]; // Accessing the data in column-major order
            }
        }
    }

    // Free the memory allocated for the MATLAB variable
    Mat_VarFree(matvar);

    // Closing MAT file
    Mat_Close(matfp);

    return data; // Return the populated 2D vector
}

int main(){
    string filename = "angle_45.mat";
    string varname = "AuData";

    srand(time(NULL));

    vector<vector<double>> AuData = ReadMatfile(filename, varname);
    if (AuData.empty())
    {

        cout << "Error reading file" << endl;
    }

    size_t Lcommon = AuData.size();

    size_t Nsteps = floor((Lcommon - framesize) / time_step + 1);

    vector<int> hist_bins(359);
    iota(hist_bins.begin(), hist_bins.end(), 1);

    vector<vector<double>> EstimatesFrames(Nsteps, vector<double>(numzones * NumEstPerZone, NAN));
    vector<vector<double>> Estimated_DOAs(Nsteps, vector<double>(NSources, NAN));

    int windowlength = 5;

    for (int istep = 0; istep < Nsteps; istep++){
        size_t start = istep * time_step;
        size_t end = min(start + framesize, AuData.size());

        vector<vector<double>> current_data;
        for(int i = start; i < end; i++){
            current_data.push_back(AuData[i]);
        }

        vector<double> current_estimates = DOAsCurrentFrame(current_data, fftsize, ignore_bins, numzones, MMCC, zone_step, freqzonebins, c, Fs, NumEstPerZone, Ra);
        
        copy(current_estimates.begin(), current_estimates.end(), EstimatesFrames[istep].begin());

        if(blockcounter < BlocksizeFrames){
            NoEstFr[blockcounter] = current_estimates.size();  // Όπως στο MATLAB
            Blockestimates.insert(Blockestimates.end(), current_estimates.begin(), current_estimates.end());
            blockcounter++;
        }else{
            NoEstFr[blockcounter] = current_estimates.size();  // Όπως στο MATLAB
            Blockestimates.insert(Blockestimates.end(), current_estimates.begin(), current_estimates.end());
        }

        vector<double> PlainHist, SmoothedHist;
        DOA_Histogram(Blockestimates, windowlength, PlainHist, SmoothedHist);
        vector<double> block_DOA = processHist(SmoothedHist, NSources, gaussWidth);
        Estimated_DOAs[istep] = block_DOA;

        int removeCount = NoEstFr.front();
        if(removeCount <= Blockestimates.size()) Blockestimates.erase(Blockestimates.begin(), Blockestimates.begin() + removeCount);
       // NoEstFr.erase(NoEstFr.begin());   
    }

    // // Εκτύπωση του NoEstFr
    // cout << "\nNoEstFr:\n";
    // for(int i = 0; i < NoEstFr.size(); i++) {
    //     cout << NoEstFr[i] << " "<<endl;
    // }

    // return 0;
}
