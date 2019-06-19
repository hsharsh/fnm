#include "utilities.cpp"
#include "fe_functions.cpp"

// Note that with 0-indexed arrays, first dof is node*(#dofs), next is node*(#dofs)+1, and so on
void boundary_conditions(VectorXd &vn, VectorXd &vn1){
  double bc;
  // Syntax: LinSpaced(increment,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXd bottom(38);
  bottom <<     3,   4,   9,  16,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76, 214,
215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 351, 352, 353, 354, 355,
356, 357, 358, 359, 360, 361;
  bottom = bottom.array()-1;
  bc = -0.01;
  for(int i=0; i < bottom.size();i++){
    long node_bc = bottom(i);
    long idof = node_bc*2+1;
    vn(idof) = bc;
    vn1(idof) = bc;
  }

  VectorXd top(38);
  top <<   14,  15,  19,  21, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 464,
 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 533, 534, 535, 536, 537,
 538, 539, 540, 541, 542, 543;
  top = top.array()-1;
  bc = 0.01;
  for(int i=0; i < top.size();i++){
    long node_bc = top(i);
    long idof = node_bc*2+1;
    vn(idof) = bc;
    vn1(idof) = bc;
  }
}



void temporary_bc(VectorXd &vn, VectorXd &vn1){
  double bc;
  VectorXd left = VectorXd::LinSpaced(21,0,420);
  bc = 0.01;
  for(int i=0; i < left.size();i++){
    long node_bc = left(i);
    long idof = node_bc*2;
    vn(idof) = bc;
    vn1(idof) = bc;
  }
}

int main(int argc, char* argv[]){
  double tmax = 10, dt = 0.05;
  if(argc == 2){
    dt = atof(argv[1]);
  }
  else if(argc == 3){
    dt = atof(argv[1]);
    tmax = atof(argv[2]);
  }

  MatrixXd elements = load_csv<MatrixXd>("/home/hsharsh/fnm/elements.inp");
  MatrixXd nodes = load_csv<MatrixXd>("/home/hsharsh/fnm/nodes.inp");

  // All coefficients are decreased by one for consistency with 0-indexing
  MatrixXd conn = elements(all,seq(1,last)).array()-1;
  MatrixXd x = nodes(all,seq(1,last));

  long nnod = x.rows();
  long nelm = conn.rows();
  long nodi = nnod*2 + 1;
  double E = 1, nu = 0, rho = 1, eta = 0;

  VectorXd un = VectorXd::Zero(nnod*2), un1 = VectorXd::Zero(nnod*2);
  VectorXd vn = VectorXd::Zero(nnod*2), vn1 = VectorXd::Zero(nnod*2), an1 = VectorXd::Zero(nnod*2);

  // boundary_conditions(vn, vn1);

  double t = 0, n = 1, nf = 1;
  while(t <= tmax){
    VectorXd  mg = VectorXd::Zero((nodi-1));
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
    an1 = mg.array().inverse()*(fg-fi-lcg).array();
    vn1 = vn + an1*dt;

    boundary_conditions(vn,vn1);
    // if(t < 4.0){
    //   temporary_bc(vn,vn1);
    // }

    un1 = un + vn1*dt;

    if(n >= nf){
      MatrixXd xdef = MatrixXd::Zero(x.rows(),x.cols()+1);

      MatrixXd u = MatrixXd::Zero((int)un1.rows()/2,3), v = MatrixXd::Zero((int)vn1.rows()/2,3), a = MatrixXd::Zero((int)an1.rows()/2,3);
      xdef(all,0) = x(all,0) + un(seq(0,last,2));
      xdef(all,1) = x(all,1) + un(seq(1,last,2));

      u(all,0) = un1(seq(0,last,2));    u(all,1) = un1(seq(1,last,2));
      v(all,0) = vn1(seq(0,last,2));    v(all,1) = vn1(seq(1,last,2));
      a(all,0) = an1(seq(0,last,2));    a(all,1) = an1(seq(1,last,2));


      string filename = "x0";   filename.append(to_string((long)(t*1e5)));   filename.append(".vtk");
      vtkwrite(filename,conn,xdef,u,v,a);
      n = 1;
    }
    else
      n++;

    // File writing operations

    cout << "Time: " << t << endl;
    t = t+dt;
    un = un1;
    vn = vn1;
  }
}
