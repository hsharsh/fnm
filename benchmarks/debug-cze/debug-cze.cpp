// Boundary conditions for simple case used for debugging
// Activate only the boundary_conditions function. No temporary bc is needed

void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  VectorXi left(4);
  left << 1, 3, 5, 7;
  left = left.array()-1;
  bcx = 0;
  for(int i=0; i < left.size();++i){
    int node_bc = left(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn1(idof) = bcx;
  }

  VectorXi bottom(2);
  bottom << 1, 2;
  bottom = bottom.array()-1;
  bcx = 0;
  bcy = 0;
  for(int i=0; i < bottom.size();++i){
    int node_bc = bottom(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof) = bcx;
    vn1(idof+1) = bcy;
  }

  VectorXi top(2);
  top << 7, 8;
  top = top.array()-1;
  bcx = 0;
  bcy = 0.1;
  for(int i=0; i < top.size();++i){
    int node_bc = top(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof) = bcx;
    vn1(idof+1) = bcy;
  }
}

void setcze(vector<pair<int,int> > &cze){
  for(int i = 0; i < cze.size(); ++i){
    cze[i].first = 0;
    cze[i].second = -1;
  }

  VectorXi cohesive(1);
  cohesive << 2;
  cohesive = cohesive.array()-1;

  for(int i=0; i < cohesive.size();++i){
    cze[cohesive[i]].first = 1;
    cze[cohesive[i]].second = 1;
  }
}

tmax = 10
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
nlyrs = 2
init_c = 0
j_tol = 100
delta = 0.2
smax = 0.001
ashear = 1
