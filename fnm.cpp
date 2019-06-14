#include "utilities.cpp"

vector<double> xgp = {-sqrt(3.0/5.0), 0, sqrt(3.0/5.0)};
vector<double> wgp = {5.0/9.0, 8.0/9.0, 5.0/9.0};
long ngp = wgp.size();

void boundary_conditions(VectorXd &vn, VectorXd &vn1){
  double bc;
  VectorXd left = VectorXd::LinSpaced(21,1,421);
  bc = -0.01;
  for(int i=0; i < left.size();i++){
    long node_bc = left(i);
    long idof = node_bc*2-1;
    vn(idof) = bc;
    vn1(idof) = bc;
  }

  VectorXd right = VectorXd::LinSpaced(21,21,441);
  bc = 0.01;
  for(int i=0; i < right.size();i++){
    long node_bc = right(i);
    long idof = node_bc*2-1;
    vn(idof) = bc;
    vn1(idof) = bc;
  }
}

MatrixXd quad_kl(MatrixXd &xv, double &area, double E, double nu){
  MatrixXd kl = MatrixXd::Zero(4*2,4*2);

  for(int i = 0; i < ngp; i++){
    for(int j = 0; j < ngp; j++){
        double r = xgp[i], s = xgp[j];
        MatrixXd B1(2,4), jac(2,2), B0(3,4), B3 = MatrixXd::Zero(4,8), Bjac = MatrixXd::Zero(4,4);
        B1 << -(1-s)/4,  (1-s)/4, (1+s)/4, -(1+s)/4,
              -(1-r)/4,  -(1+r)/4,  (1+r)/4,  (1-r)/4;
        jac = (B1*xv).transpose();  // Check this for problems

        B0 << 1, 0, 0, 0,
              0, 0, 0, 1,
              0, 1, 1, 0;
        // Define Bjac
        // Bjac(1:2,1:2) = inv(jac);
        // Bjac(3:4,3:4) = inv(jac);

        // Define B3
        // B3(1:2, 1:2:end) = B1;
        // B3(3:4, 2:2:end) = B1;
        MatrixXd B = B0*Bjac*B3;
        MatrixXd D(3,3);
        D << 1-nu, 0, 0,
              0, 1-nu, 0,
              0, 0, 1-2*nu;
        D = E/((1+nu)*(1-2*nu))*D.array();

        kl = kl + B.transpose()*D*B*jac.determinant() * wgp[i] * wgp[j];
        area = area + jac.determinant()* wgp[i] * wgp[j];
    }
  }
  return kl;
}


void assemble_fi(VectorXd &fi, VectorXd &un, MatrixXd &x, MatrixXd &conn, double E, double nu){
  double area = 0;
  long nelm = conn.rows();
  for(int i = 0; i < nelm; i++){
    if(true){ // Change this to reflect True for elements which don't have floating nodes activated
      VectorXd nodes = conn(0,all);
      MatrixXd xv = x(nodes,all);
      MatrixXd kl = quad_kl(xv,area,E,nu);
      vector<double> dof = {nodes(0)*2-1, nodes(0)*2,
                            nodes(1)*2-1, nodes(1)*2,
                            nodes(2)*2-1, nodes(2)*2,
                            nodes(3)*2-1, nodes(3)*2};
      VectorXd u = un(dof);

      fi(dof) = fi(dof) + (kl*u);
    }
    else{
      // Calls for floating nodes
    }
  }
  cout << area << endl;
}

int main(){
  MatrixXd conn = load_csv<MatrixXd>("/home/hsharsh/fnm/elements.inp");
  MatrixXd x = load_csv<MatrixXd>("/home/hsharsh/fnm/nodes.inp");
  conn = conn(all,seq(1,last));
  x = x(all,seq(1,last));

  long nnod = x.rows();
  long nelm = conn.rows();
  long nodi = nnod*2 + 1;
  double E = 1, nu = 0, rho = 1, eta = 0;

  VectorXd un = VectorXd::Zero(nnod*2), un1 = VectorXd::Zero(nnod*2);
  VectorXd vn = VectorXd::Zero(nnod*2), vn1 = VectorXd::Zero(nnod*2), an1 = VectorXd::Zero(nnod*2);

  boundary_conditions(vn, vn1);


  double t = 0, tmax = 0.04;
  double dt = 0.05;

  while(t <= tmax){
    VectorXd  mg = VectorXd::Ones((nodi-1));
    VectorXd  fi = VectorXd::Zero((nodi-1));
    VectorXd  fg = VectorXd::Zero((nodi-1));
    VectorXd lcg = VectorXd::Zero((nodi-1));

    // Define crack
      // Crack definition

    // Add floating nodes to the global matrices
      // Floating nodes definition

    // Linearized Global Stiffness matrix assembly
    assemble_fi(fi, un, x, conn, E, nu);

    // Linearized Global damping matrix
      // Damping matrix assembly

    // Linearized Global Mass matrix assembly
    // assemble_mg();

    // Enforce boundary conditions
    boundary_conditions(vn,vn1);

    // Solver
    an1 = mg.array().inverse()*(fg-fi-lcg).array();
    vn1 = vn + an1*dt;
    un1 = un + vn1*dt;

    t = t+dt;
  }
}
