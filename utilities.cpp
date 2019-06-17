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

void vtkwrite(string filename, MatrixXd connectivity, MatrixXd position, MatrixXd u, MatrixXd v, MatrixXd a){
  string path = "/home/hsharsh/fnm/data/";
  path.append(filename);
  ofstream writer(path);

  if(writer.is_open()){
    writer << "# vtk DataFile Version 1.0\n";
    writer << "VTK - Matlab export\n";
    writer << "ASCII\n";

    long precision = 8;

    writer << "\nDATASET POLYDATA\n";
    writer << "POINTS " << position.rows() << " float\n";
    for(int i = 0; i < position.rows(); i++)
        writer << position(i,all) << "\n";

    writer << "\nPOLYGONS " << connectivity.rows() << " " << 5*connectivity.rows() << "\n";
    for(int i = 0; i < connectivity.rows(); i++)
        writer << "4 " << connectivity(i,all).array()-1 << "\n";


    writer << "\nPOINT_DATA " << position.rows() <<  " \n";
    writer << "VECTORS displacement float\n";
    for(int i = 0; i < u.rows();i++)
        writer << u(i,all) << " \n";

    writer << "\nVECTORS velocity float\n";
    for(int i = 0; i < v.rows();i++)
        writer << v(i,all) << "\n";

    writer << "\nVECTORS acceleration float\n";
    for(int i = 0; i < a.rows();i++)
        writer << a(i,all) << "\n";

    // writer << "SCALARS stress float\n";
    // writer << "LOOKUP_TABLE deafult\n";
    //
    // for(int i = 0; i < stress.size();i++)
    //     writer << stress(i,all) << "; \n";

    writer.close();
  }
  else{
    cerr << "Cannot open file";
  }

}
