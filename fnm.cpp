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


        Bjac(seq(0,1),seq(0,1)) = jac.inverse();
        Bjac(seq(2,3),seq(2,3)) = jac.inverse();


        // Define B3
        B3(seq(0,1),seq(0,last,2)) = B1;
        B3(seq(2,3),seq(1,last,2)) = B1;

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
  // cout << kl << endl << endl;
  return kl;
}

MatrixXd quad_ml(MatrixXd &xv, double &area, double rho){
  MatrixXd ml = MatrixXd::Zero(4*2,4*2);

  double area_elem = 0;
  for(int i = 0; i < ngp; i++){
    for(int j = 0; j < ngp; j++){
        double r = xgp[i], s = xgp[j];
        MatrixXd B1(2,4), jac(2,2);
        B1 << -(1-s)/4,  (1-s)/4, (1+s)/4, -(1+s)/4,
              -(1-r)/4,  -(1+r)/4,  (1+r)/4,  (1-r)/4;
        jac = (B1*xv).transpose();  // Check this for problems

        area_elem = area_elem + jac.determinant()* wgp[i] * wgp[j];
    }
  }
  area = area + area_elem;
  ml = rho*(area_elem/4)*(MatrixXd::Identity(4*2,4*2)).array();

  // cout << ml << endl;
  return ml;
}

void assemble_mg(VectorXd &mg, MatrixXd &x, MatrixXd &conn, double rho){
  double area = 0;
  long nelm = conn.rows();
  for(int i = 0; i < nelm; i++){
    if(true){ // Change this to reflect True for elements which don't have floating nodes activated
      VectorXd nodes = conn(0,all);
      MatrixXd xv = x(nodes,all);
      MatrixXd ml = quad_ml(xv,area,rho);
      vector<double> dof = {nodes(0)*2-1, nodes(0)*2,
                            nodes(1)*2-1, nodes(1)*2,
                            nodes(2)*2-1, nodes(2)*2,
                            nodes(3)*2-1, nodes(3)*2};

      mg(dof) = mg(dof) + ml.diagonal();
    }
    else{
      // Calls for floating nodes
    }
  }
  // cout << area << endl;
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
      cout << u << endl;

      fi(dof) = fi(dof) + (kl*u);
    }
    else{
      // Calls for floating nodes
    }
  }
  // cout << fi.head(15) << endl;
  // cout << area << endl;
}

int main(){
  MatrixXd elements = load_csv<MatrixXd>("/home/hsharsh/fnm/elements.inp");
  MatrixXd nodes = load_csv<MatrixXd>("/home/hsharsh/fnm/nodes.inp");

  MatrixXd conn = elements(all,seq(1,last));
  MatrixXd x = nodes(all,seq(1,last));

  long nnod = x.rows();
  long nelm = conn.rows();
  long nodi = nnod*2 + 1;
  double E = 1, nu = 0, rho = 1, eta = 0;

  VectorXd un = VectorXd::Zero(nnod*2), un1 = VectorXd::Zero(nnod*2);
  VectorXd vn = VectorXd::Zero(nnod*2), vn1 = VectorXd::Zero(nnod*2), an1 = VectorXd::Zero(nnod*2);

  boundary_conditions(vn, vn1);


  double t = 0, tmax = 0.05;
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
    assemble_mg(mg, x, conn, rho);

    // Enforce boundary conditions
    boundary_conditions(vn,vn1);

    // Solver

    cout << fi.head(15) << endl;
    an1 = mg.array().inverse()*(fg-fi-lcg).array();
    vn1 = vn + an1*dt;

    cout << an1.head(15) << endl;

    boundary_conditions(vn,vn1);
    // cout << vn1.head(15);
    un1 = un + vn1*dt;
    // cout << un1.head(15);

    MatrixXd xdef = MatrixXd::Zero(x.rows(),x.cols()+1);

    MatrixXd u = MatrixXd::Zero((int)un1.rows()/2,3), v = MatrixXd::Zero((int)vn1.rows()/2,3), a = MatrixXd::Zero((int)an1.rows()/2,3);
    xdef(all,0) = x(all,0) + un(seq(0,last,2));
    xdef(all,1) = x(all,1) + un(seq(1,last,2));

    u(all,0) = un1(seq(0,last,2));    u(all,1) = un1(seq(1,last,2));
    v(all,0) = vn1(seq(0,last,2));    v(all,1) = vn1(seq(1,last,2));
    a(all,0) = an1(seq(0,last,2));    a(all,1) = an1(seq(1,last,2));


    string filename = "x0";   filename.append(to_string((long)(t*1e5)));   filename.append(".vtk");
    vtkwrite(filename,conn,xdef,u,v,a);

    cout << "Time: " << t << endl;
    t = t+dt;
    un = un1;
    vn = vn1;
  }
}
