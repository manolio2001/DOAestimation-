#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include "DOA.hpp"

using namespace std;
bool compare(vector<double> &a, vector<double> &b)
{
   return a[0] < b[0];
}

vector<vector<double>> locate_zones(vector<double> &r,int size, double threshold)
{
   vector<vector<double>> candidate_zones;
   for(int i = 1; i < size; i++)
   {
         
         int index = i - 1 ;
         if(r[index] > threshold)
         {
            candidate_zones.push_back({r[index], static_cast<double>(i),static_cast<double> (1)});
         }
      
   }
   sort(candidate_zones.begin(), candidate_zones.end(), compare);
   return candidate_zones;
}