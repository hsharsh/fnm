#ifndef UTILITIES
#define UTILITIES

  #include<bits/stdc++.h>
  #include<eigen-unstable/Eigen/Dense>

  using namespace Eigen;
  using namespace std;

  vector<double> xgp = {-sqrt(3.0/5.0), 0, sqrt(3.0/5.0)};
  vector<double> wgp = {5.0/9.0, 8.0/9.0, 5.0/9.0};
  int ngp = wgp.size();

  // Reading CSV files into an Eigen MatrixXd variable
  template<typename M, typename type>
  M load_csv (const string & path) {
      ifstream indata;
      indata.open(path);
      string line;
      vector<type> values;
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

  void print_vector(vector<vector<double> > v){
    int ni = v.size(), nj = v[0].size();
    cout << endl;
    for (int i = 0; i < ni; i++){
      for (int j = 0; j < nj; j++){
        cout << v[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }

  void print_vector(vector<vector<int> > v){
    int ni = v.size(), nj = v[0].size();
    cout << endl;
    for (int i = 0; i < ni; i++){
      for (int j = 0; j < nj; j++){
        cout << v[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
  void print_vector(vector<double> v){
    int ni = v.size();
    cout << endl;
    for (int i = 0; i < ni; i++){
        cout << v[i] << " ";
    }
    cout << endl;
  }

  void print_vector(vector<int> v){
    int ni = v.size();
    cout << endl;
    for (int i = 0; i < ni; i++){
        cout << v[i] << " ";
    }
    cout << endl;
  }
  void vtkwrite(string filename, vector<vector<int> > fl_connectivity, long s, MatrixXd position, MatrixXd u, MatrixXd v, MatrixXd a){
    string path = "/home/hsharsh/fnm/data/";
    path.append(filename);
    ofstream writer(path);

    if(writer.is_open()){
      writer << "# vtk DataFile Version 1.0\n";
      writer << "VTK - Matlab export\n";
      writer << "ASCII\n";

      int precision = 8;

      writer << "\nDATASET POLYDATA\n";
      writer << "POINTS " << position.rows() << " float\n";
      for(int i = 0; i < position.rows(); i++)
        writer << position(i,all) << "\n";

        writer << "\nPOLYGONS " << fl_connectivity.size() << " " << s+fl_connectivity.size() << "\n"; //<< connectivity.rows() << 5*connectivity.rows()+4*fl_connectivity.size() << "\n";
        // for(int i = 0; i < connectivity.rows(); i++)
        //   writer << "4 " << connectivity(i,all) << "\n";

        for(int i = 0; i < fl_connectivity.size(); i++){
          writer << fl_connectivity[i].size() << " ";
          for (int j = 0; j < fl_connectivity[i].size(); j++){
            writer << fl_connectivity[i][j] << " ";
          }
          writer << "\n";
        }

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
#endif
