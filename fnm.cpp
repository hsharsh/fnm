#include "utilities.cpp"
#include "fe_functions.cpp"
#include "bound_cond.cpp"
#include "fn_functions.cpp"

int main(int argc, char* argv[]){
  double tmax = 0.05, dt = 0.05;
  if(argc == 2){
    dt = atof(argv[1]);
  }
  else if(argc == 3){
    dt = atof(argv[1]);
    tmax = atof(argv[2]);
  }
  MatrixXi elements = load_csv<MatrixXi,int>("/home/hsharsh/fnm/elements.inp");
  MatrixXd nodes = load_csv<MatrixXd,double>("/home/hsharsh/fnm/nodes.inp");

  // All coefficients are decreased by one for consistency with 0-indexing
  MatrixXi conn = elements(all,seq(1,last)).array()-1;
  MatrixXd x = nodes(all,seq(1,last));

  int nnod = x.rows();
  int nelm = conn.rows();
  // int nodi = nnod*2 + 1;
  double E = 1, nu = 0, rho = 1, eta = 0;

  // nnod*4 to include the maximum number of floating nodes considering only type 1 and type 2 kind of division
  VectorXd un = VectorXd::Zero(nnod*4), un1 = VectorXd::Zero(nnod*4);
  VectorXd vn = VectorXd::Zero(nnod*4), vn1 = VectorXd::Zero(nnod*4), an1 = VectorXd::Zero(nnod*4);

  vector<int> discont(nnod,0);
  map <int,element> fn_elements;

  boundary_conditions(vn, vn1);

  double t = 0, n = 1, nf = 1;
  while(t <= tmax){
    VectorXd  mg = VectorXd::Zero(nnod*2);
    VectorXd  fi = VectorXd::Zero(nnod*2);
    VectorXd  fg = VectorXd::Zero(nnod*2);
    VectorXd lcg = VectorXd::Zero(nnod*2);

    // Linearized Global Stiffness matrix assembly
    assemble_fi(fi, un, x, conn, E, nu);

    // Linearized Global damping matrix
    assemble_lcg(lcg, vn, x, conn, eta);

    // Linearized Global Mass matrix assembly
    assemble_mg(mg, x, conn, rho);

    boundary_conditions(vn,vn1);

    // cout << "fi" << endl;
    // cout << fi.head(15) << endl << endl;
    //
    // cout << "mg" << endl;
    // cout << mg.head(15) << endl << endl;
    // Solver

    an1(seq(0,nnod*2-1)) = mg.array().inverse()*(fg-fi-lcg).array();
    vn1 = vn + an1*dt;

    boundary_conditions(vn,vn1);
    // if(t < 4.0){
    //   temporary_bc(vn,vn1);
    // }

    un1 = un + vn1*dt;

    // File writing operations
    if(n >= nf){
      MatrixXd xdef = MatrixXd::Zero(nnod,4);

      MatrixXd u = MatrixXd::Zero(nnod,4), v = MatrixXd::Zero(nnod,4), a = MatrixXd::Zero(nnod,4);
      xdef(all,0) = x(seq(0,nnod-1),0) + un(seq(0,nnod*2-1,2));
      xdef(all,1) = x(seq(0,nnod-1),1) + un(seq(1,nnod*2-1,2));

      u(all,0) = un1(seq(0,nnod*2-1,2));    u(all,1) = un1(seq(1,nnod*2-1,2));
      v(all,0) = vn1(seq(0,nnod*2-1,2));    v(all,1) = vn1(seq(1,nnod*2-1,2));
      a(all,0) = an1(seq(0,nnod*2-1,2));    a(all,1) = an1(seq(1,nnod*2-1,2));


      string filename = "x0";   filename.append(to_string((int)(t*1e5)));   filename.append(".vtk");
      vtkwrite(filename,conn,xdef,u,v,a);
      n = 1;
    }
    else
      n++;

    // Define crack
      // Crack definition

    // Add floating nodes to the global matrices
    floating_nodes(discont, fn_elements, conn, x, nnod);
    cout << "Time: " << t << endl;
    t = t+dt;
    un = un1;
    vn = vn1;
  }
}
