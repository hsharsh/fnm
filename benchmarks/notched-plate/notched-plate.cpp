// Boundary conditions for plate with a hole and streching on top and bottom
// Activate only the boundary_conditions function. No temporary bc is needed
void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  // Syntax: LinSpaced(#elements,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXi bottom(38);
  bottom <<  7,   8,  15,  16, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
120, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 316, 317,
318, 319, 320, 321, 322, 323;
  bottom = bottom.array()-1;
  bcx = 0;
  bcy = -1e-4;
  for(int i=0; i < bottom.size();++i){
    long node_bc = bottom(i);
    long idof = node_bc*2;
    vn(idof) = bcx;
    vn1(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof+1) = bcy;
  }

  VectorXi top(39);
  top <<   1,   4,  11,  12,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,
 69, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 274, 275,
276, 277, 278, 279, 280, 281, 282;
  top = top.array()-1;
  bcx = 0;
  bcy = 1e-4;
  for(int i=0; i < top.size();++i){
    long node_bc = top(i);
    long idof = node_bc*2;
    vn(idof) = bcx;
    vn1(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof+1) = bcy;
  }
}

void crack_def(vector<int> &discont, map<int,element> &fn_elements){
  VectorXi cracked(5);
  cracked << 784, 783, 782, 781, 780;
  cracked = cracked.array()-1;
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {NAN, 0.5, NAN, 0.5};
  }
}

void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
  VectorXi cracked(5);
  cracked << 784, 783, 782, 781, 780;
  cracked = cracked.array()-1;
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {NAN, 0.5, NAN, 0.5};
    VectorXi nodes = conn(cracked[i],all);
    for(int j = 0; j < 4; ++j){
      if(!isnan(fn_elements[cracked[i]].edge[j])){
        if(cparam.find(make_pair(nodes(j),nodes((j+1)%4))) == cparam.end() && cparam.find(make_pair(nodes((j+1)%4),nodes(j))) == cparam.end()){
          cparam[make_pair(nodes(j),nodes((j+1)%4))] = fn_elements[cracked[i]].edge[j];
        }
      }
    }
  }
}

tmax = 0.03
dt = 2e-2
E = 1
nu = 0
rho = 1
alpha = 0.05
sy = 1
ar_tol = 2.5e-4
rf = 2
tc = 0.1
sampling_rate = 100
