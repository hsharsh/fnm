// Boundary conditions for plate with a hole and streching on top and bottom
// Activate only the boundary_conditions function. No temporary bc is needed
void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  // Syntax: LinSpaced(#elements,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXi bottom(91);
  bottom <<     21,  22,  25,  28, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605,
 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 860, 861,
 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877,
 878, 879, 880, 881, 882, 883, 884, 885, 951, 952, 953, 954, 955, 956, 957, 958,
 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 971, 972, 973, 974,
 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985;
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

  VectorXi top(91);
  top <<         1,   5,  10,  12, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132,
   133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148,
   149, 150, 151, 152, 153, 154, 155, 264, 265, 266, 267, 268, 269, 270, 271, 272,
   273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288,
   289, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330,
   331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341;
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
  for(int i = 9878; i <= 9886; ++i){
    cracked.push_back(i-1);
  }
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {NAN, 0.0976, NAN, 0.9024};
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
  cracked.push_back(9886);
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 3;
    fn_elements[cracked[i]].edge = {NAN, -0.0976, NAN, 0.9024};
    VectorXi nodes = conn(cracked[i],all);
    for(int j = 0; j < 4; ++j){
      if(!isnan(fn_elements[cracked[i]].edge[j])){
        if(cparam.find(make_pair(nodes(j),nodes((j+1)%4))) == cparam.end() && cparam.find(make_pair(nodes((j+1)%4),nodes(j))) == cparam.end()){
          cparam[make_pair(nodes(j),nodes((j+1)%4))] = abs(fn_elements[cracked[i]].edge[j]);
        }
      }
    }
  }
  discont[9887] = 4;
}

tmax = 15
dt = 5e-4
E = 1
nu = 0
rho = 1
alpha = 0
sy = 6e8
ar_tol = 2.5e-4
rf = 1
tc = 0.1
srate = 10
nlyrs = 2
init_c = 1
j_tol = 1e-7

// Old shiz below this line.

void crack_def(vector<int> &discont, map<int,element> &fn_elements){
  vector<int> cracked;
  for(int i = 9878; i <= 9887; ++i)
    cracked.push_back(i-1);

  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {NAN, 0.0976, NAN, 0.9024};
  }
}
