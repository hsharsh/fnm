#include "utilities.cpp"
#include "fn_functions.cpp"
#include "fe_functions.cpp"
#include "bound_cond.cpp"
#include "crack_def.cpp"

int main(int argc, char* argv[]){
  double dt, tmax, E, nu, rho, alpha, nf, tc, rf;
  if(system("exec rm -r /home/hsharsh/fnm/data/*"))
    cerr << "Error clearing old data" << endl;
  MatrixXi elements = load_csv<MatrixXi,int>("/home/hsharsh/fnm/elements.inp");
  MatrixXd nodes = load_csv<MatrixXd,double>("/home/hsharsh/fnm/nodes.inp");

  // All coefficients are decreased by one for consistency with 0-indexing
  MatrixXd x = MatrixXd::Zero(nodes.rows()*2+elements.rows()*4,2);
  x(seq(0,nodes.rows()-1),all) = nodes(all,seq(1,last));
  MatrixXi conn = elements(all,seq(1,last)).array()-1;

  int nnod = nodes.rows();
  int ndof = 2*nnod;
  int nelm = elements.rows();

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
      if(name == "nf")
        nf = stoi(value);
      if(name == "tc")
        tc = stof(value);
      if(name == "rf")
        rf = stoi(value);
    }
    cout << endl;
  }
  else{
    cerr << "Couldn't open config file for reading." << endl;
    return 0;
  }

  // nnod*4 to include the maximum number of floating nodes considering only type 1 and type 2 kind of division
  VectorXd un = VectorXd::Zero(nnod*2+nelm*8), un1 = VectorXd::Zero(nnod*2+nelm*8);
  VectorXd vn = VectorXd::Zero(nnod*2+nelm*8), vn1 = VectorXd::Zero(nnod*2+nelm*8);
  VectorXd an1 = VectorXd::Zero(nnod*2+nelm*8), stress = VectorXd::Zero(nnod*2+nelm*8);

  vector<int> discont(nnod,0);
  map <int,element> fn_elements;

  // boundary_conditions(vn,vn1);
  // boundary_conditions(un,un1,vn,vn1);

  double t = 0, n = nf, ti = 0;
  bool crack_active = 1;
  while(t <= tmax){
    cout << "Time: " << t << endl;

    VectorXd mg = VectorXd::Zero(ndof);
    VectorXd fi = VectorXd::Zero(ndof);
    VectorXd fg = VectorXd::Zero(ndof);
    VectorXd lcg = VectorXd::Zero(ndof);


    // Linearized Global Stiffness matrix assembly
    assemble_fi(fi, un, x, conn, discont, fn_elements, E, nu);

    // Linearized Global Mass matrix assembly
    assemble_mg(mg, x, conn, discont, fn_elements, rho, ndof);

    // Linearized Global damping matrix
    assemble_lcg(lcg, vn, x, conn, discont, fn_elements, rho, alpha);

    // BC for fg
    boundary_conditions(un,un1,vn,vn1,fg);
    // Solver
    an1(seq(0,ndof-1)) = mg.array().inverse()*(fg-fi-lcg).array();
    vn1 = vn + an1*dt;

    // BC for velocity
    boundary_conditions(un,un1,vn,vn1,fg);

    // if(t < 4.0){
    //   temporary_bc(vn,vn1);
    // }

    un1 = un + vn1*dt;

    // BC for displacement
    boundary_conditions(un,un1,vn,vn1,fg);

    // Define crack
    if(abs(t-0) < 1e-5){
      cout << "Crack intialized" << endl;
      crack_def(discont,fn_elements);
    }

    if(crack_active){
      int c = stress_based_crack(discont, fn_elements, conn, x, un1, ndof, E, nu);
      if (c == 1){
        crack_active = 0;
      }
    }
    else{
      if(ti < tc)
        ti+=dt;
      else{
        ti = 0;
        crack_active = 1;
      }
    }


    // Add floating nodes to the global matrices
    floating_nodes(discont, fn_elements, conn, x, un1, ndof);

    // Remove elements which have area smaller than a certain tolerance (defined in utilities.cpp)
    remove_singular_elements(fn_elements,x);

    // Compute nodal properties
    compute_properties(stress, discont, fn_elements, conn, x, un1, ndof, E, nu);

    // File writing operations
    vector<vector<int> > fl_conn;
    MatrixXd conn_w;
    long s = 0;
    if(n >= nf){
      MatrixXd xdef = MatrixXd::Zero(ndof/2,3);

      MatrixXd u = MatrixXd::Zero(ndof/2,3), v = MatrixXd::Zero(ndof/2,3), a = MatrixXd::Zero(ndof/2,3), f = MatrixXd::Zero(ndof/2,3);

      VectorXd str = stress(seq(0,(ndof/2)-1));

      xdef(all,0) = x(seq(0,(ndof/2)-1),0) + un1(seq(0,ndof-1,2));
      xdef(all,1) = x(seq(0,(ndof/2)-1),1) + un1(seq(1,ndof-1,2));

      u(all,0) = un1(seq(0,ndof-1,2));    u(all,1) = un1(seq(1,ndof-1,2));
      v(all,0) = vn1(seq(0,ndof-1,2));    v(all,1) = vn1(seq(1,ndof-1,2));
      a(all,0) = an1(seq(0,ndof-1,2));    a(all,1) = an1(seq(1,ndof-1,2));
      f(all,0) = fi(seq(0,ndof-1,2));     f(all,1) = fi(seq(1,ndof-1,2));

      for (int i = 0; i < nelm; ++i){
        if(discont[i]){
          for(int j = 0; j < fn_elements[i].conn.size(); ++j){
            if(fn_elements[i].active[j]){
              fl_conn.push_back(fn_elements[i].conn[j]);
              s+=3;
            }

          }

        }
        else{
          vector<int> elem_conn(conn.cols());
          for(int j = 0; j < conn.cols(); ++j){
            elem_conn[j] = conn(i,j);
          }
          fl_conn.push_back(elem_conn);
          s+=4;
        }
      }
      string filename = "x0";   filename.append(to_string((int)(t*1e5)));   filename.append(".vtk");
      vtkwrite(filename,fl_conn,s,xdef,u,v,a,str,f);
      n = 1;
    }
    else
      ++n;

    if(cracked == 1){
      cracked = 2;
      dt/=rf;
    }

    t = t+dt;
    un = un1;
    vn = vn1;
  }
}
