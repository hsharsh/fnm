#include "utilities.cpp"
#include "fn_functions.cpp"
#include "fe_functions.cpp"
#include "bound_cond.cpp"
#include "crack_def.cpp"
int main(int argc, char* argv[]){
  double dt = 0.05, tmax = 1.0;
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
  MatrixXd x = MatrixXd::Zero(nodes.rows()*2+elements.rows()*4,2);
  x(seq(0,nodes.rows()-1),all) = nodes(all,seq(1,last));
  MatrixXi conn = elements(all,seq(1,last)).array()-1;

  int nnod = nodes.rows();
  int ndof = 2*nnod;
  int nelm = elements.rows();
  double E = 1, nu = 0, rho = 1, eta = 0;

  // nnod*4 to include the maximum number of floating nodes considering only type 1 and type 2 kind of division
  VectorXd un = VectorXd::Zero(nnod*2+nelm*8), un1 = VectorXd::Zero(nnod*2+nelm*8);
  VectorXd vn = VectorXd::Zero(nnod*2+nelm*8), vn1 = VectorXd::Zero(nnod*2+nelm*8), an1 = VectorXd::Zero(nnod*2+nelm*8);

  vector<int> discont(nnod,0);
  map <int,element> fn_elements;

  boundary_conditions(vn, vn1);

  double t = 0, n = 1, nf = 1;
  while(t <= tmax){
    cout << "Time: " << t << endl;

    // VectorXd mg = VectorXd::Zero(ndof);
    VectorXd mg = VectorXd::Ones(ndof);
    VectorXd fi = VectorXd::Zero(ndof);
    VectorXd fg = VectorXd::Zero(ndof);
    VectorXd lcg = VectorXd::Zero(ndof);

    // // Linearized Global Stiffness matrix assembly
    assemble_fi(fi, un, x, conn, discont, fn_elements, E, nu);
    //
    // // Linearized Global damping matrix
    assemble_lcg(lcg, vn, x, conn, discont, fn_elements, eta);
    //
    // // Linearized Global Mass matrix assembly
    assemble_mg(mg, x, conn, discont, fn_elements, rho);

    boundary_conditions(vn,vn1);

    // cout << "fi" << endl;
    // cout << fi.head(15) << endl << endl;
    //
    // cout << "mg" << endl;
    // cout << mg.head(15) << endl << endl;
    // Solver

    an1(seq(0,ndof-1)) = mg.array().inverse()*(fg-fi-lcg).array();
    vn1 = vn + an1*dt;

    boundary_conditions(vn,vn1);
    // if(t < 4.0){
    //   temporary_bc(vn,vn1);
    // }

    un1 = un + vn1*dt;

    // File writing operations
    if(n >= nf){
      MatrixXd xdef = MatrixXd::Zero(ndof/2,3);

      MatrixXd u = MatrixXd::Zero(ndof/2,3), v = MatrixXd::Zero(ndof/2,3), a = MatrixXd::Zero(ndof/2,3);
      xdef(all,0) = x(seq(0,(ndof/2)-1),0) + un(seq(0,ndof-1,2));
      xdef(all,1) = x(seq(0,(ndof/2)-1),1) + un(seq(1,ndof-1,2));

      u(all,0) = un1(seq(0,ndof-1,2));    u(all,1) = un1(seq(1,ndof-1,2));
      v(all,0) = vn1(seq(0,ndof-1,2));    v(all,1) = vn1(seq(1,ndof-1,2));
      a(all,0) = an1(seq(0,ndof-1,2));    a(all,1) = an1(seq(1,ndof-1,2));


      string filename = "x0";   filename.append(to_string((int)(t*1e5)));   filename.append(".vtk");
      vtkwrite(filename,conn,xdef,u,v,a);
      n = 1;
    }
    else
      n++;

    // Define crack
    if(t == 0.1){
      cout << "Crack added" << endl;
      crack_def(discont,fn_elements);
    }

    // Add floating nodes to the global matrices
    floating_nodes(discont, fn_elements, conn, x, ndof);

    t = t+dt;
    un = un1;
    vn = vn1;
  }
}
