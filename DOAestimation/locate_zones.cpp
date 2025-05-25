#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include "DOA.hpp"

using namespace std;

// r: vector<double>, size = N
// size: int, number of elements in r
// threshold: double
// Returns: vector<vector<double>> where each row is {value, index, 0}, sorted by value (ascending)
vector<vector<double>> locate_zones(vector<double> &r, int size, double threshold) {
    vector<tuple<double, int, int>> candidate_zones;
    for (int i = 0; i < size; ++i) {
        if (r[i] > threshold) {
            candidate_zones.emplace_back(r[i], i, 0); // 0 for column, for compatibility
        }
    }
    // Sort by value (ascending)
    sort(candidate_zones.begin(), candidate_zones.end());
    // Convert to vector<vector<double>>
    vector<vector<double>> candidate_zones_sorted;
    for (const auto& t : candidate_zones) {
        candidate_zones_sorted.push_back({get<0>(t), static_cast<double>(get<1>(t)), static_cast<double>(get<2>(t))});
    }
    return candidate_zones_sorted;
}