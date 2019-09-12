#ifndef UTILITIES
#define UTILITIES

  #include<bits/stdc++.h>
  #include<algorithm>
  #include<sys/stat.h>
  #include<eigen-unstable/Eigen/Dense>
  #include<omp.h>

  using namespace Eigen;
  using namespace std;

  #define pi 			3.141592653593
  #define eps 		0.0000001

  vector<double> xgp = {-sqrt(3.0/5.0), 0, sqrt(3.0/5.0)};
  vector<double> wgp = {5.0/9.0, 8.0/9.0, 5.0/9.0};
  int ngp = wgp.size();
  double ar_tol,sy;
  int cracked = 0;
  string path = "/home/hsharsh/fnm/data";

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

  void print(vector<vector<double> > v){
    int ni = v.size();
    cout << endl;
    for (int i = 0; i < ni; ++i){
      for (int j = 0; j < v[i].size(); ++j){
        cout << v[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }

  void print(vector<vector<int> > v){
    int ni = v.size();
    cout << endl;
    for (int i = 0; i < ni; ++i){
      for (int j = 0; j < v[i].size(); ++j){
        cout << v[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
  void print(vector<double> v){
    int ni = v.size();
    cout << endl;
    for (int i = 0; i < ni; ++i){
        cout << v[i] << " ";
    }
    cout << endl;
  }

  void print(vector<int> v){
    int ni = v.size();
    cout << endl;
    for (int i = 0; i < ni; ++i){
        cout << v[i] << " ";
    }
    cout << endl;
  }

  bool load_config(double &dt, double &tmax, double &E, double &nu, double &rho, double &alpha, double &sy, double &ar_tol, int &sampling_rate, int &tc, int &rf, int &nlayers){
    ifstream cFile("parameters.cfg");
    if (cFile.is_open()){
        cout << "Running with parameters:" << endl;
        string line;
        while(getline(cFile, line)){
        line.erase(remove_if(line.begin(), line.end(), ::isspace),line.end());
        if(line[0] == '#' || line.empty())
          continue;
        auto delimiterPos = line.find("=");
        string name = line.substr(0, delimiterPos);
        string value = line.substr(delimiterPos + 1);
        cout << name << ":\t" << value << endl;
        if(name == "dt")
          dt = stof(value);
        if(name == "tmax")
          tmax = stof(value);
        if(name == "E")
          E = stof(value);
        if(name == "nu")
          nu = stof(value);
        if(name == "rho")
          rho = stof(value);
        if(name == "alpha")
          alpha = stof(value);
        if(name == "sy")
          sy = stof(value);
        if(name == "ar_tol")
          ar_tol = stof(value);
        if(name == "sampling_rate")
          sampling_rate = stoi(value);
        if(name == "tc")
          tc = stof(value);
        if(name == "rf")
          rf = stoi(value);
        if(name == "nlayers")
          nlayers = stoi(value);
      }
      cout << endl;
      return 1;
    }
    return 0;
  }

  void vtkwrite(string filename, vector<vector<int> > fl_connectivity, long s, MatrixXd position, MatrixXd u, MatrixXd v, MatrixXd a, VectorXd stress, MatrixXd fi){
    string p = path;
    p.append(filename);
    ofstream writer(p);

    if(writer.is_open()){
      writer << "# vtk DataFile Version 1.0\n";
      writer << "VTK - Matlab export\n";
      writer << "ASCII\n";

      int precision = 8;

      writer << "\nDATASET POLYDATA\n";
      writer << "POINTS " << position.rows() << " float\n";
      for(int i = 0; i < position.rows(); ++i)
        writer << position(i,all) << "\n";

        writer << "\nPOLYGONS " << fl_connectivity.size() << " " << s+fl_connectivity.size() << "\n"; //<< connectivity.rows() << 5*connectivity.rows()+4*fl_connectivity.size() << "\n";
        // for(int i = 0; i < connectivity.rows(); ++i)
        //   writer << "4 " << connectivity(i,all) << "\n";

        for(int i = 0; i < fl_connectivity.size(); ++i){
          writer << fl_connectivity[i].size() << " ";
          for (int j = 0; j < fl_connectivity[i].size(); ++j){
            writer << fl_connectivity[i][j] << " ";
          }
          writer << "\n";
        }

      writer << "\nPOINT_DATA " << position.rows() <<  " \n";
      writer << "VECTORS displacement float\n";
      for(int i = 0; i < u.rows();++i)
        writer << u(i,all) << " \n";

      writer << "\nVECTORS velocity float\n";
      for(int i = 0; i < v.rows();++i)
        writer << v(i,all) << "\n";

      writer << "\nVECTORS acceleration float\n";
      for(int i = 0; i < a.rows();++i)
        writer << a(i,all) << "\n";

      writer << "\nVECTORS reaction float\n";
      for(int i = 0; i < fi.rows();++i)
        writer << fi(i,all) << "\n";

      writer << "SCALARS stress float\n";
      writer << "LOOKUP_TABLE deafult\n";

      for(int i = 0; i < stress.size();++i)
          writer << stress(i) << "\n";

      writer.close();
    }
    else{
      cerr << "Cannot open file" << endl;
      exit(0);
    }
  }

  void write_j(string filename, double t, double j_integral){
    string p = path;
    p.append(filename);
    ofstream writer(p,ofstream::out | ofstream::app);

    if(writer.is_open()){
      writer << t << " " << j_integral << "\n";
      writer.close();
    }
    else{
      cerr << "Cannot open file" << endl;
      exit(0);
    }
  }

  // if(system("exec rm -r /home/hsharsh/fnm/data/*"))
  //   cerr << "Error clearing old data" << endl;
#endif
