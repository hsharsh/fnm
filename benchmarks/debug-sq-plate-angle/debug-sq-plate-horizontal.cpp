// Boundary conditions for simple case used for debugging
// Activate only the boundary_conditions function. No temporary bc is needed

void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  VectorXi bottom = VectorXi::LinSpaced(4,1,4);
  bottom = bottom.array()-1;
  bcx = 0.5e-4;
  bcy = 0.86602545e-4;
  for(int i=0; i < bottom.size();++i){
    int node_bc = bottom(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof) = bcx;
    vn1(idof+1) = bcy;
  }
  // VectorXi top = VectorXi::LinSpaced(4,13,16);
  // top = top.array()-1;
  // bcx = 0;
  // bcy = 1e-4;
  // for(int i=0; i < top.size();++i){
  //   int node_bc = top(i);
  //   int idof = node_bc*2;
  //   vn(idof) = bcx;
  //   vn(idof+1) = bcy;
  //   vn1(idof) = bcx;
  //   vn1(idof+1) = bcy;
  // }
  VectorXi right = VectorXi::LinSpaced(4,4,16);
  right = right.array()-1;
  bcx = -0.5e-4;
  bcy = -0.86602545e-4;
  for(int i=0; i < right.size();++i){
    int node_bc = right(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof) = bcx;
    vn1(idof+1) = bcy;
  }
}

void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
  vector <int> cracked{3};
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 3;
    fn_elements[cracked[i]].edge = {NAN, -0.5, NAN, 0.5};
    VectorXi nodes = conn(cracked[i],all);
    for(int j = 0; j < 4; ++j){
      if(!isnan(fn_elements[cracked[i]].edge[j])){
        if(cparam.find(make_pair(nodes(j),nodes((j+1)%4))) == cparam.end() && cparam.find(make_pair(nodes((j+1)%4),nodes(j))) == cparam.end()){
          cparam[make_pair(nodes(j),nodes((j+1)%4))] = abs(fn_elements[cracked[i]].edge[j]);
        }
      }
    }
  }

  // Transition element
  discont[4] = 4;
}

tmax = 0.03
dt = 0.02
E = 1
nu = 0
rho = 1
alpha = 0
sy = 6e8
ar_tol = 1e-8
rf = 2
tc = 0.1
srate = 1
nlyrs = 1
init_c = 1
j_tol = 100

// Angle plate path for abaqus
 4, 1, 17, 9, 11, 10, 20, 18, 13, 12, 19, 15, 6, 5, 16, 8, 2, 3, 21
 
// Old shiz below this line
