// Boundary conditions for simple case used for debugging
// Activate only the boundary_conditions function. No temporary bc is needed

void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  VectorXi right(8);
  right << 1,   4,   6,   7,  83,  84, 119, 120;
  right = right.array()-1;
  bcx = 0;
  bcy = 0;
  for(int i=0; i < right.size();++i){
    int node_bc = right(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof) = bcx;
    vn1(idof+1) = bcy;
  }

  VectorXi leftbottom(1);
  leftbottom << 10;
  leftbottom = leftbottom.array()-1;
  bcx = 0;
  bcy = -0.001;
  for(int i=0; i < leftbottom.size();++i){
    int node_bc = leftbottom(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof) = bcx;
    vn1(idof+1) = bcy;
  }

  VectorXi lefttop(1);
  lefttop << 12;
  lefttop = lefttop.array()-1;
  bcx = 0;
  bcy = 0.001;
  for(int i=0; i < lefttop.size();++i){
    int node_bc = lefttop(i);
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

  VectorXi cohesive = VectorXi::LinSpaced(35,211,245);
  cohesive = cohesive.array()-1;

  for(int i=0; i < cohesive.size();++i){
    cze[cohesive[i]].first = 1;
    cze[cohesive[i]].second = 1;
  }

  VectorXi remove = VectorXi::LinSpaced(15,336,350);
  remove = remove.array()-1;

  for(int i=0; i < remove.size();++i){
    cze[remove[i]].first = 2;
    cze[remove[i]].second = -1;
  }
}

tmax = 1
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
delta = 0.001
smax = 1e-7
ashear = 0.1
