#include "utilities.cpp"
#include "materials.cpp"
#include "fn_functions.cpp"
#include "cze.cpp"
#include "fe_functions.cpp"
#include "bound_cond.cpp"
#include "j_int_formulation.cpp"
#include "crack_def.cpp"

int main(int argc, char* argv[]){
  double dt, tmax, E, nu, rho, alpha, tc;
  int srate, rf, nlyrs;
  bool init_c = 1;

  // Creating a folder with time-stamp for saving data
  path.append(to_string((int)time(0)).append("/"));
  // experimental::filesystem::create_directories(path); // Make this work
  mkdir(path.c_str(),0777);

  // Copy the boundary_conditions
  ifstream  src("bound_cond.cpp");
  ofstream  dst(path+(string)"boundary_conditions.cpp");
  dst << src.rdbuf();
  src.close();
  dst.close();

  MatrixXi elements = load_csv<MatrixXi,int>(project_directory+"elements.inp");
  MatrixXd nodes = load_csv<MatrixXd,double>(project_directory+"nodes.inp");

  // All coefficients are decreased by one for consistency with 0-indexing
  MatrixXd x = MatrixXd::Zero(nodes.rows()*2+elements.rows()*4,2);
  x(seq(0,nodes.rows()-1),all) = nodes(all,seq(1,last));
  MatrixXi conn = elements(all,seq(1,last)).array()-1;
  // vector<pair<pair<int,int>,pair<int,int> > > crack_tip;

  int nnod = nodes.rows();
  int ndof = 2*nnod;
  int nelm = elements.rows();

  if(!load_config(tmax, dt, E, nu, rho, alpha, sy, ar_tol, tc, j_tol, delta, smax, ashear, srate, rf, nlyrs, init_c)){
    cerr << "Couldn't open config file for reading." << endl;
    return 0;
  }

  // Compute neighbours for J-integral computation
  vector<vector<int> > neighbours(nnod*2+nelm*8);
  compute_neighbours(neighbours, conn);

  // nnod*4 to include the maximum number of floating nodes considering only type 1 and type 2 kind of division
  VectorXd un = VectorXd::Zero(nnod*2+nelm*8), un1 = VectorXd::Zero(nnod*2+nelm*8);
  VectorXd vn = VectorXd::Zero(nnod*2+nelm*8), vn1 = VectorXd::Zero(nnod*2+nelm*8);
  VectorXd an1 = VectorXd::Zero(nnod*2+nelm*8), stress = VectorXd::Zero(nnod*2+nelm*8);

  vector<int> discont(nnod,0);
  /* Legend for discont values
  0 -> Standard Quad
  1 -> Cracked/Unprocessed
  2 -> Cracked/Processed
  3 -> Crack tip/Transition element
  4 -> Element adajacent to crack tip changed for compatibility/Transition element
  5 -> 3/Processed
  6 -> 4/Processed
  */
  map <int,element> fn_elements;
  map <pair<int,int>,double> cparam;
  map <pair<int,int>,pair<int,int> > allocated_nodes;
  vector<pair<double,double> > matprop(nelm);
  vector<pair<int,int> > cze(nelm);

  setproperties(matprop, E, nu);

  setcze(cze);
  // boundary_conditions(vn,vn1);
  // boundary_conditions(un,un1,vn,vn1);

  // bool first_iteration = 0;
  double t = 0, n = srate, ti = 0;
  bool crack_active = 1;

  // Define crack
  if(init_c){
    cout << "Crack intialized" << endl;
    crack_def(discont,fn_elements, conn, cparam);
  }

  // Add floating nodes to the global matrices
  floating_nodes(discont, fn_elements, allocated_nodes, conn, x, un1, ndof);

  // Remove elements which have area smaller than a certain tolerance (defined in utilities.cpp)
  remove_singular_elements(fn_elements,x);



  while(t <= tmax){
    // cout << "Time: " << t << endl;

    VectorXd mg = VectorXd::Zero(ndof);
    VectorXd fi = VectorXd::Zero(ndof);
    VectorXd fg = VectorXd::Zero(ndof);
    VectorXd lcg = VectorXd::Zero(ndof);


    // Linearized Global Stiffness matrix assembly
    assemble_fi(fi, un, x, conn, discont, fn_elements, matprop, cze);

    // Linearized Global Mass matrix assembly
    assemble_mg(mg, x, conn, discont, fn_elements, cze, rho, ndof);

    // Linearized Global damping matrix
    assemble_lcg(lcg, vn, x, conn, discont, fn_elements, cze, rho, alpha);

    // BC for fg
    // if(t < 1){
    //   boundary_conditions(un,un1,vn,vn1,fg);
    // }
    boundary_conditions(un,un1,vn,vn1,fg);

    // Solver
    an1(seq(0,ndof-1)) = mg.array().inverse()*(fg-fi-lcg).array();
    vn1 = vn + an1*dt;

    // BC for velocity
    // if(t < 1){
    //   boundary_conditions(un,un1,vn,vn1,fg);
    // }
    boundary_conditions(un,un1,vn,vn1,fg);

    // if(t < 4.0){
    //   temporary_bc(vn,vn1);
    // }

    un1 = un + vn1*dt;

    // BC for displacement
    // if(t < 1){
    //   boundary_conditions(un,un1,vn,vn1,fg);
    // }
    boundary_conditions(un,un1,vn,vn1,fg);

    // Define crack
    // if(abs(t-4.0) < 1e-8 && init_c){
    //   cout << "Crack intialized" << endl;
    //   crack_def_new(discont,fn_elements, conn, cparam);
    // }

    // if(abs(t-0) > 1e-5 && first_iteration == 0){
    //   crack_tip.push_back(make_pair(make_pair(910,911),make_pair(428,449)));
    //   first_iteration = 1;
    // }

    // Crack propagation

    if(crack_active){
      double j_integral = compute_j(neighbours, conn, x, un1, discont, fn_elements, cparam, nnod, E, nu, nlyrs);
      if(j_integral == NAN){
        cerr << "No J-integral computed." << endl;
      }
      j_integral = abs(j_integral); // Converting j-integral calculated opposeite to crack direction to positive.
      int c = j_based_crack(discont, neighbours, fn_elements, cparam, conn, x, un1, ndof, E, nu, j_integral, nlyrs);
      if (c == 1){
        crack_active = 0;
      }
    }

    if(!crack_active){
      if(ti < tc)
        ti+=dt;
      else{
        ti = 0;
        crack_active = 1;
      }
    }

    // Add floating nodes to the global matrices
    floating_nodes(discont, fn_elements, allocated_nodes, conn, x, un1, ndof);

    // Remove elements which have area smaller than a certain tolerance (defined in utilities.cpp)
    remove_singular_elements(fn_elements,x);

    // Compute nodal properties
    compute_properties(stress, discont, fn_elements, conn, x, un1, ndof, E, nu, cze);

    // File writing operations
    vector<vector<int> > fl_conn;
    MatrixXd conn_w;
    long s = 0;
    // cout << "x: " << endl;
    // cout << x(seq(1,ndof/2),all) << endl << endl;
    if(n >= srate){
      write_j("j_int.m",t,compute_j(neighbours, conn, x, un1, discont, fn_elements, cparam, nnod, E, nu, nlyrs));
      cout << "Time: " << t << endl;

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
        if(cze[i].first == 2){
          continue;
        }
        if(discont[i]){
          for(int j = 0; j < fn_elements[i].conn.size(); ++j){
            if(fn_elements[i].active[j]){
              fl_conn.push_back(fn_elements[i].conn[j]);
              s+=fn_elements[i].conn[j].size();
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
      string filename = "x0";   filename.append(to_string((long long)(t*1e5)));   filename.append(".vtk");
      vtkwrite(filename,fl_conn,s,xdef,u,v,a,str,f);
      n = 1;
    }
    else
      ++n;

    if(cracked == 1){
      cracked = 2;
      dt/=rf;
      srate = 5;
    }
    // cout << dt << endl;
    t = t+dt;
    un = un1;
    vn = vn1;
  }
}
