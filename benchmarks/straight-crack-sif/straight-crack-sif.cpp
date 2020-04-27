// Boundary conditions for plate with a hole and streching on top and bottom
// Activate only the boundary_conditions function. No temporary bc is needed
void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  // Syntax: LinSpaced(#elements,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXi bottom = VectorXi::LinSpaced(10,1,10);
  bottom = bottom.array()-1;
  bcx = 0;
  bcy = -1e6/9;
  for(int i=0; i < bottom.size();++i){
    long node_bc = bottom(i);
    long idof = node_bc*2;
    // vn(idof) = bcx;
    // vn1(idof) = bcx;
    // vn(idof+1) = bcy;
    // vn1(idof+1) = bcy;
    fg(idof+1) = bcy;
  }

  VectorXi top = VectorXi::LinSpaced(10,371,380);
  top = top.array()-1;
  bcx = 0;
  bcy = 1e6/9;
  for(int i=0; i < top.size();++i){
    long node_bc = top(i);
    long idof = node_bc*2;
    // vn(idof) = bcx;
    // vn1(idof) = bcx;
    // vn(idof+1) = bcy;
    // vn1(idof+1) = bcy;
    fg(idof+1) = bcy;
  }

  VectorXi left = VectorXi::LinSpaced(38,1,371);
  left = left.array()-1;
  bcx = 0;
  for(int i=0; i < left.size();++i){
    long node_bc = left(i);
    long idof = node_bc*2;
    un(idof) = bcx;
    un1(idof) = bcx;
  }

  VectorXi right = VectorXi::LinSpaced(38,10,380);
  right = right.array()-1;
  bcx = 0;
  for(int i=0; i < right.size();++i){
    long node_bc = right(i);
    long idof = node_bc*2;
    un(idof) = bcx;
    un1(idof) = bcx;
  }
}

void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
  vector <int> cracked;
  for (int i = 163; i <= 164; ++i){
    cracked.push_back(i-1);
  }
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {NAN, 0.5, NAN, 0.5};
    VectorXi nodes = conn(cracked[i],all);
    for(int j = 0; j < 4; ++j){
      if(!isnan(fn_elements[cracked[i]].edge[j])){
        if(cparam.find(make_pair(nodes(j),nodes((j+1)%4))) == cparam.end() && cparam.find(make_pair(nodes((j+1)%4),nodes(j))) == cparam.end()){
          cparam[make_pair(nodes(j),nodes((j+1)%4))] = abs(fn_elements[cracked[i]].edge[j]);
        }
      }
    }
  }
  cracked.clear();
  cracked.push_back(164);
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
  discont[165] = 4;
}

tmax = 0.01
dt = 1e-7
E = 200e9
nu = 0.3
rho = 8000
alpha = 0
sy = 6e8
ar_tol = 1e-8
rf = 1
tc = 0.1
srate = 10
nlyrs = 2
init_c = 1
j_tol = 100
