// Boundary conditions for square plate with velocity on both sides
// Activate the temporary_bc for 4 seconds. Initiate crack just before 15 seconds.

void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  // VectorXi left(9);
  // left << 203,  204,  205, 1384, 1385, 1386, 1408, 1409, 1410;
  VectorXi left(4);
  left << 239,  240,  241, 2833;
  left = left.array()-1;
  bcx = 0;
  bcy = 0;
  for(int i=0; i < left.size();++i){
    int node_bc = left(i);
    int idof = node_bc*2;
    un(idof) = bcx;
    un(idof+1) = bcy;
    un1(idof) = bcx;
    un1(idof+1) = bcy;
  }
  // VectorXi right(9);
  // right <<  90,   91,   92, 1180, 1181, 1188, 1189, 1196, 1197;
  VectorXi right(4);
  right <<  440,  441,  442, 5938;
  right = right.array()-1;
  bcy = 0;
  for(int i=0; i < right.size();++i){
    int node_bc = right(i);
    int idof = node_bc*2;
    un(idof+1) = bcy;
    un1(idof+1) = bcy;
  }

  VectorXi top(4);
  top << 1099, 10495, 10498, 10501;
  top = top.array()-1;
  bcy = -1e-4;
  for(int i=0; i < top.size();++i){
    int node_bc = top(i);
    int idof = node_bc*2;
    vn(idof+1) = bcy;
    vn1(idof+1) = bcy;
  }
}

void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
  vector <int> cracked;
  for(int i = 5475; i <= 5493; ++i)
    cracked.push_back(i-1);
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {NAN, 0.1905, NAN, 0.8095};
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
  cracked.push_back(5474);
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 3;
    fn_elements[cracked[i]].edge = {NAN, -0.1905, NAN, 0.8095};
    VectorXi nodes = conn(cracked[i],all);
    for(int j = 0; j < 4; ++j){
      if(!isnan(fn_elements[cracked[i]].edge[j])){
        if(cparam.find(make_pair(nodes(j),nodes((j+1)%4))) == cparam.end() && cparam.find(make_pair(nodes((j+1)%4),nodes(j))) == cparam.end()){
          cparam[make_pair(nodes(j),nodes((j+1)%4))] = abs(fn_elements[cracked[i]].edge[j]);
        }
      }
    }
  }
  discont[5473] = 4;
}

tmax = 300
dt = 1e-2
E = 1
nu = 0.3
rho = 1
alpha = 1
sy = 100
ar_tol = 1e-4
rf = 2
tc = 0.1
srate = 10
nlyrs = 2
init_c = 1
j_tol = 1e-8


// Old shiz below this line

 void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
   vector<int> cracked;
   for(int i = 5474; i <= 5493; ++i)
     cracked.push_back(i-1);
   for (int i = 0; i < cracked.size(); ++i){
     discont[cracked[i]] = 1;
     fn_elements[cracked[i]].edge = {NAN, 0.1905, NAN, 0.8095};
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

 tmax = 200
 dt = 2e-2
 E = 1
 nu = 0
 rho = 1
 alpha = 0
 sy = 1
 ar_tol = 2.5e-4
 rf = 2
 tc = 0.1
 srate = 10
 nlyrs = 2
 init_c = 1
 j_tol = 1e-7
