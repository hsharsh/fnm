// Boundary conditions for simple case used for debugging
// Activate only the boundary_conditions function. No temporary bc is needed

void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  VectorXi bottom = VectorXi::LinSpaced(4,1,4);
  bottom = bottom.array()-1;
  bcx = 0;
  bcy = -1e-2;
  for(int i=0; i < bottom.size();++i){
    int node_bc = bottom(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof) = bcx;
    vn1(idof+1) = bcy;
  }
  VectorXi top = VectorXi::LinSpaced(4,13,16);
  top = top.array()-1;
  bcx = 0;
  bcy = 1e-2;
  for(int i=0; i < top.size();++i){
    int node_bc = top(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof) = bcx;
    vn1(idof+1) = bcy;
  }
  VectorXi right = VectorXi::LinSpaced(4,4,16);
  right = right.array()-1;
  bcx = 0;
  for(int i=0; i < right.size();++i){
    int node_bc = right(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn1(idof) = bcx;
  }
}
// Growing Crack
void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
  vector <int> cracked{3};

  // Standard
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

  //Crack tip
  cracked.clear();
  cracked.push_back(4);
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

  //Transition
  discont[5] = 4;
}


// Complete Crack
void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
  vector <int> cracked{3,4,5};

  // Standard
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
}

// Mixed growing crack

void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
  vector <int> cracked{0};

  // Standard
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {NAN, NAN, 0.5, 0.5};
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
  cracked.push_back(3);
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = { 0.5, 0.5, NAN, NAN};
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
  cracked.push_back(4);
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = { NAN, 0.5, NAN, 0.5};
    VectorXi nodes = conn(cracked[i],all);
    for(int j = 0; j < 4; ++j){
      if(!isnan(fn_elements[cracked[i]].edge[j])){
        if(cparam.find(make_pair(nodes(j),nodes((j+1)%4))) == cparam.end() && cparam.find(make_pair(nodes((j+1)%4),nodes(j))) == cparam.end()){
          cparam[make_pair(nodes(j),nodes((j+1)%4))] = abs(fn_elements[cracked[i]].edge[j]);
        }
      }
    }
  }

  //Crack tip
  cracked.clear();
  cracked.push_back(5);
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 3;
    fn_elements[cracked[i]].edge = {NAN, NAN, -0.5, 0.5};
    VectorXi nodes = conn(cracked[i],all);
    for(int j = 0; j < 4; ++j){
      if(!isnan(fn_elements[cracked[i]].edge[j])){
        if(cparam.find(make_pair(nodes(j),nodes((j+1)%4))) == cparam.end() && cparam.find(make_pair(nodes((j+1)%4),nodes(j))) == cparam.end()){
          cparam[make_pair(nodes(j),nodes((j+1)%4))] = abs(fn_elements[cracked[i]].edge[j]);
        }
      }
    }
  }

  //Transition
  discont[8] = 4;
}

// Mixed full crack

void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
  vector <int> cracked{0,5};

  // Standard
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {NAN, NAN, 0.5, 0.5};
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
  cracked.push_back(3);
  cracked.push_back(8);
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = { 0.5, 0.5, NAN, NAN};
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
  cracked.push_back(4);
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = { NAN, 0.5, NAN, 0.5};
    VectorXi nodes = conn(cracked[i],all);
    for(int j = 0; j < 4; ++j){
      if(!isnan(fn_elements[cracked[i]].edge[j])){
        if(cparam.find(make_pair(nodes(j),nodes((j+1)%4))) == cparam.end() && cparam.find(make_pair(nodes((j+1)%4),nodes(j))) == cparam.end()){
          cparam[make_pair(nodes(j),nodes((j+1)%4))] = abs(fn_elements[cracked[i]].edge[j]);
        }
      }
    }
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
nlyrs = 1
init_c = 1
j_tol = 100

// Square plate path for abaqus
 4, 1, 17, 9, 11, 10, 20, 18, 13, 12, 19, 15, 6, 5, 16, 8, 2, 3, 21
