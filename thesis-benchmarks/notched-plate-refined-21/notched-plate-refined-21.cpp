// Boundary conditions for plate with a hole and streching on top and bottom
// Activate only the boundary_conditions function. No temporary bc is needed
void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  // Syntax: LinSpaced(#elements,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXi bottom(47);
  bottom <<     21,  22,  25,  28, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317,
 318, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 481, 482,
 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497;
  bottom = bottom.array()-1;
  bcx = 0;
  bcy = 0;
  for(int i=0; i < bottom.size();++i){
    long node_bc = bottom(i);
    long idof = node_bc*2;
    vn(idof) = bcx;
    vn1(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof+1) = bcy;
  }

  VectorXi top(47);
  top <<      1,   5,  10,  12,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,
 87,  88,  89,  90,  91, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153,
154, 155, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181;
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

void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
  vector <int> cracked;
  for(int i = 2549; i <= 2552; ++i){
    cracked.push_back(i-1);
  }
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {NAN, 0.3228, NAN, 0.6772};
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
  cracked.push_back(2552);
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 3;
    fn_elements[cracked[i]].edge = {NAN, -0.3228, NAN, 0.6772};
    VectorXi nodes = conn(cracked[i],all);
    for(int j = 0; j < 4; ++j){
      if(!isnan(fn_elements[cracked[i]].edge[j])){
        if(cparam.find(make_pair(nodes(j),nodes((j+1)%4))) == cparam.end() && cparam.find(make_pair(nodes((j+1)%4),nodes(j))) == cparam.end()){
          cparam[make_pair(nodes(j),nodes((j+1)%4))] = abs(fn_elements[cracked[i]].edge[j]);
        }
      }
    }
  }
  discont[2553] = 4;
}

tmax = 15
dt = 5e-4
E = 1
nu = 0
rho = 1
alpha = 1
sy = 6e8
ar_tol = 1e-6
rf = 1
tc = 0.1
srate = 10
nlyrs = 2
init_c = 1
j_tol = 1e-7
