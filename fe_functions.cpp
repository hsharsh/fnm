#include "utilities.cpp"

MatrixXd quad_kl(MatrixXd &xv, double &area, double E, double nu){
  MatrixXd kl = MatrixXd::Zero(4*2,4*2);

  for(int i = 0; i < ngp; i++){
    for(int j = 0; j < ngp; j++){
        double r = xgp[i], s = xgp[j];
        MatrixXd B0(2,4), jac(2,2), B1(3,4), B3 = MatrixXd::Zero(4,8), B2 = MatrixXd::Zero(4,4);
        B0 << -(1-s)/4,  (1-s)/4, (1+s)/4, -(1+s)/4,
              -(1-r)/4,  -(1+r)/4,  (1+r)/4,  (1-r)/4;
        jac = (B0*xv);

        B1 << 1, 0, 0, 0,
              0, 0, 0, 1,
              0, 1, 1, 0;

        B2(seq(0,1),seq(0,1)) = jac.inverse();
        B2(seq(2,3),seq(2,3)) = jac.inverse();

        // Define B3
        B3(seq(0,1),seq(0,last,2)) = B0;
        B3(seq(2,3),seq(1,last,2)) = B0;

        MatrixXd B = B1*B2*B3;
        MatrixXd D(3,3);
        D << 1-nu, 0, 0,
              0, 1-nu, 0,
              0, 0, (1-2*nu)/2;
        D = E/((1+nu)*(1-2*nu))*D.array();

        kl = kl + B.transpose()*D*B*jac.determinant() * wgp[i] * wgp[j];
        area = area + jac.determinant()* wgp[i] * wgp[j];
    }
  }
  return kl;
}

MatrixXd quad_ml(MatrixXd &xv, double &area, double rho){
  MatrixXd ml = MatrixXd::Zero(4*2,4*2);

  double area_elem = 0;
  for(int i = 0; i < ngp; i++){
    for(int j = 0; j < ngp; j++){
        double r = xgp[i], s = xgp[j];
        MatrixXd B0(2,4), jac(2,2);
        B0 << -(1-s)/4,  (1-s)/4, (1+s)/4, -(1+s)/4,
              -(1-r)/4,  -(1+r)/4,  (1+r)/4,  (1-r)/4;
        jac = (B0*xv);

        area_elem = area_elem + jac.determinant()* wgp[i] * wgp[j];
    }
  }
  area = area + area_elem;
  ml = rho*(area_elem/4)*(MatrixXd::Identity(4*2,4*2)).array();

  return ml;
}

void assemble_mg(VectorXd &mg, MatrixXd &x, MatrixXi &conn, double rho){
  double area = 0;
  long nelm = conn.rows();
  for(int i = 0; i < nelm; i++){
    if(true){ // Change this to reflect True for elements which don't have floating nodes activated
      VectorXi nodes = conn(i,all);
      MatrixXd xv = x(nodes,all);
      MatrixXd ml = quad_ml(xv,area,rho);

      vector<int> dof = {   nodes(0)*2, nodes(0)*2+1,
                            nodes(1)*2, nodes(1)*2+1,
                            nodes(2)*2, nodes(2)*2+1,
                            nodes(3)*2, nodes(3)*2+1};

      mg(dof) = mg(dof) + ml.diagonal();
    }
    else{
      // Calls for floating nodes
    }
  }
  // cout << area << endl;
}

void assemble_fi(VectorXd &fi, VectorXd &un, MatrixXd &x, MatrixXi &conn, double E, double nu){
  double area = 0;
  long nelm = conn.rows();
  for(int i = 0; i < nelm; i++){
    if(true){ // Change this to reflect True for elements which don't have floating nodes activated
      VectorXi nodes = conn(i,all);
      MatrixXd xv = x(nodes,all);
      MatrixXd kl = quad_kl(xv,area,E,nu);
      vector<int> dof = {   nodes(0)*2, nodes(0)*2+1,
                            nodes(1)*2, nodes(1)*2+1,
                            nodes(2)*2, nodes(2)*2+1,
                            nodes(3)*2, nodes(3)*2+1};
      VectorXd u = un(dof);

      fi(dof) = fi(dof) + (kl*u);
    }
    else{
      // Calls for floating nodes
    }
  }
}

void assemble_lcg(VectorXd &lcg, VectorXd &vn, MatrixXd &x, MatrixXi &conn, double eta){
  double area = 0;
  long nelm = conn.rows();
  for(int i = 0; i < nelm; i++){
    if(true){ // Change this to reflect True for elements which don't have floating nodes activated
      VectorXi nodes = conn(i,all);
      MatrixXd xv = x(nodes,all);
      MatrixXd cl = (MatrixXd::Ones(4*2,4*2)).array()*eta;
      vector<int> dof = {   nodes(0)*2, nodes(0)*2+1,
                            nodes(1)*2, nodes(1)*2+1,
                            nodes(2)*2, nodes(2)*2+1,
                            nodes(3)*2, nodes(3)*2+1};
      VectorXd v = vn(dof);

      lcg(dof) = lcg(dof) + (cl*v);
    }
    else{
      // Calls for floating nodes
    }
  }
}
