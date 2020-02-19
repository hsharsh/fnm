#include "utilities.cpp"

inline MatrixXd constitutive(double E, double nu){
  MatrixXd D(3,3);
  D << 1-nu, nu, 0,
        nu, 1-nu, 0,
        0, 0, (1-2*nu);
  D = E/((1+nu)*(1-2*nu))*D.array();
  return D;
}

void setproperties(vector<pair<double,double> > &matprop, double E, double nu){
  for(int i = 0; i < matprop.size(); ++i){
    matprop[i].first = E;
    matprop[i].second = nu;
  }
}
