#include<bits/stdc++.h>
#include<eigen-unstable/Eigen/Dense>

using namespace Eigen;
using namespace std;

// Reading CSV files into an Eigen MatrixXd variable
template<typename M>
M load_csv (const string & path) {
    ifstream indata;
    indata.open(path);
    string line;
    vector<double> values;
    uint rows = 0;  // Unsigned int variable used. Replace with long long in case of overflow because of larger files
    while (getline(indata, line)) {
        stringstream lineStream(line);
        string cell;
        while (getline(lineStream, cell, ',')) {
            values.push_back(stod(cell));
        }
        ++rows;
    }

    return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
}
