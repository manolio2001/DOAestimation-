#include "DOA.hpp"
#include <stdlib.h>
#include <time.h>


// Define the variables in one source file
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
double Ra = 0.05;
int ignore_bins = 5;

//Function to read a MATLAB matrix from a .mat file and store it as a 2D vector
vector<vector<double>>ReadMatfile(string filen, string varname){
   vector<vector<double>> data; //Initialization of a 2D vector to store the data

   mat_t *matfp; //Pointer to the MAT file
   matvar_t *matvar; //Pointer to the MATLAB variable

   //Opening MAT file in read-only mode
   matfp = Mat_Open(filen.c_str(),MAT_ACC_RDONLY); 
   if ( matfp == NULL ) {
      cerr<<"Error opening file "<<filen<<endl; //Error message if the file cannot be opened
      return data; //Return the empty data vector
   }  
   //Reading the specified  variable from the MAT file
   matvar = Mat_VarRead(matfp,varname.c_str());
   if ( matvar == NULL ) {
     
      cerr<<"Error reading variable: "<<varname<<"from"<<filen<<endl; //Error message if the variable cannot be read
      Mat_Close(matfp); //Close the MAT file
      return data; //Return the empty vector
   }
   
   //Check if the variable is a 2D matrix
   if(matvar->rank == 2 && matvar){
     size_t rows = matvar->dims[0]; //Get the number of rows in the matrix
     size_t cols = matvar->dims[1]; //Get the number of columns in the matrix
     double *mat_data = (double*)matvar->data; //Get the data from the matrix and cast it to a double pointer
     
     //Reshaping the 2D vector to match the matrix dimensions
     data.resize(rows, (vector<double>(cols)));

     //Fill the 2D vector with the matrix data
     for(size_t i = 0; i < rows; i++){
        for(size_t j = 0; j < cols; j++){
           data[i][j] = mat_data[i + j*rows];//Accessing the data in column-major order
        }
     }
   }
   //Free the memory allocated for the MATLAB variable
   Mat_VarFree(matvar);

   //Closing MAT file
   Mat_Close(matfp);

   return data; //Return the populated 2D vector
}

int main(){

   //Path to the MATLAB file and the variable name
   string filename = "./DOAestimation/angle_45.mat"; //MATLAB file containing the data
   string varname = "AuData"; //Variable name inside the .mat file to extract

   //Seed the random number generator with the current time
   srand(time(NULL));

   //Read the matrix data from the specified file and variable
   vector<vector<double>> AuData = ReadMatfile(filename, varname);
   if(AuData.empty()){
     //If the matrix is empty, print an error message
      cout<<"Error reading file"<<endl;
   }

   //Is the length of the audio data 
   size_t Lcommon = AuData.size();

   //The total number of chuncks of audio contained in AuData
   size_t Nsteps = floor((Lcommon - framesize)/time_step + 1);
   
   //Selects a ramdom frame for processing
   size_t istep = rand() % Nsteps + 1;

   //Extract the current frame from the dataset
   vector<vector<double>> current_data;
   //The starting point from the extraction
   size_t start = (istep - 1) * time_step;
   for(int i = start; i <  framesize + start; i++){

      //Copies the relevant frame data
      current_data.push_back(AuData[i]);
   }

   //Perform DOA estimation on the extracted frame
  vector<double>EstimatesFrames = DOAsCurrentFrame(current_data, fftsize, ignore_bins, numzones, MMCC, zone_step, freqzonebins, c, Fs, NumEstPerZone, Ra);

   return 0;
}