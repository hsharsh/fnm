#include "utilities.cpp"

void compute_neighbours(vector<vector<int> > &neighbours, MatrixXi &conn){
    int nelm = conn.rows();
    int nodes_elem = conn.cols();

    for(int i = 0; i < nelm; ++i){
      for(int j = 0; j < nodes_elem; ++j)
        neighbours[conn(i,j)].push_back(i);
    }
}

double quad_j(MatrixXd &xv, VectorXd &u, VectorXd &lq, double E, double nu, pair<double,double> direction){
  double lj = 0;
  for (int j = 0; j < ngp; ++j){
    for (int k = 0; k < ngp; ++k){
      double r = xgp[j], s =  xgp[k];
      MatrixXd B0(2,4), jac(2,2), jacr(2,2), B1(3,4), B2 = MatrixXd::Zero(4,4), Br = MatrixXd::Zero(4,4), B3 = MatrixXd::Zero(4,8);

      B0 << -(1-s)/4,  (1-s)/4, (1+s)/4, -(1+s)/4,
            -(1-r)/4,  -(1+r)/4,  (1+r)/4,  (1-r)/4;

      jac = (B0*xv);

      double ct = direction.first/sqrt(pow(direction.first,2)+pow(direction.second,2));
      double st = direction.second/sqrt(pow(direction.first,2)+pow(direction.second,2));

      jacr << ct, st,
              -st, ct;
      Br(seq(0,1),seq(0,1)) = jacr;
      Br(seq(2,3),seq(2,3)) = jacr;

      B1 << 1, 0, 0, 0,
            0, 0, 0, 1,
            0, 0.5, 0.5, 0;

      B2(seq(0,1),seq(0,1)) = jac.inverse();
      B2(seq(2,3),seq(2,3)) = jac.inverse();

      // Define B3
      B3(seq(0,1),seq(0,last,2)) = B0;
      B3(seq(2,3),seq(1,last,2)) = B0;

      MatrixXd Bu = Br*B2*B3;
      MatrixXd Bj = jacr*jac.inverse()*B0;
      MatrixXd B = B1*Bu;

      MatrixXd strain = B*u;

      MatrixXd D = constitutive(E, nu);

      MatrixXd stress = D*strain;
      MatrixXd du = Bu*u;
      MatrixXd dq = Bj*lq;
      double w = 0.5*(stress(0)*strain(0) + stress(1)*strain(1) + 2*stress(2)*strain(2));
      lj += ( (stress(0)*du(0)*dq(0) + stress(2)*du(2)*dq(0) + stress(2)*du(0)*dq(1) + stress(1)*du(2)*dq(1) ) - w*dq(0) )*wgp[j]*wgp[k]*jac.determinant();
    }
  }
  return lj;
}

double tri_j(MatrixXd &xv, VectorXd &u, VectorXd &lq, double E, double nu, pair<double,double> direction){
  double lj = 0;

  double area_elem = 0.5*(xv(0,0)*xv(1,1)-xv(1,0)*xv(0,1) + xv(1,0)*xv(2,1)-xv(2,0)*xv(1,1) + xv(2,0)*xv(0,1)-xv(0,0)*xv(2,1));
  MatrixXd B0(2,3), jacr(2,2), B2 = MatrixXd::Zero(4,6), B1(3,4), Br = MatrixXd::Zero(4,4);
  B0 << xv(1,1)-xv(2,1), xv(2,1)-xv(0,1), xv(0,1)-xv(1,1),
          xv(2,0)-xv(1,0), xv(0,0)-xv(2,0), xv(1,0)-xv(0,0);
  B0 = B0.array()/(2*area_elem);

  double ct = direction.first/sqrt(pow(direction.first,2)+pow(direction.second,2));
  double st = direction.second/sqrt(pow(direction.first,2)+pow(direction.second,2));

  jacr << ct, st,
          -st, ct;

  Br(seq(0,1),seq(0,1)) = jacr;
  Br(seq(2,3),seq(2,3)) = jacr;

  B1 << 1, 0, 0, 0,
        0, 0, 0, 1,
        0, 0.5, 0.5, 0;

  B2(seq(0,1),seq(0,last,2)) = B0;
  B2(seq(2,3),seq(1,last,2)) = B0;

  MatrixXd B = B1*B2;

  MatrixXd strain = B*u;

  MatrixXd D = constitutive(E, nu);

  MatrixXd stress = D*strain;
  MatrixXd du = Br*B2*u;
  MatrixXd dq = jacr*B0*lq;

  double w = 0.5*(stress(0)*strain(0) + stress(1)*strain(1) + 2*stress(2)*strain(2));
  lj += ( (stress(0)*du(0)*dq(0) + stress(2)*du(2)*dq(0) + stress(2)*du(0)*dq(1) + stress(1)*du(2)*dq(1) ) - w*dq(0) )*area_elem;

  return lj;
}

double compute_j(vector<vector<int> > &neighbours, MatrixXi &conn, MatrixXd &x, VectorXd &un1, vector<int> &discont,map <int,element> &fn_elements, map <pair<int,int>,double> &cparam, int nnod, double E, double nu, int nlayers){

  double j_int = 0;
  int nelm = conn.rows();

  int tip_element = -1;
  for(int i = 0; i < nelm; ++i){
    if(discont[i] == 5){
      tip_element = i;
    }
  }

  // pair<double,double> direction = make_pair(0.866025,-0.5);
  pair<double,double> direction = make_pair(1,0);

  // cout << "NN element -> " << tip_element << endl;
  if(tip_element == -1){
    // cerr << "Crack tip is indeterminate. J-intergral not computed" << endl;
    return NAN;
  }

  set<int> domain_elem, inner_nodes, outer_nodes;

  // Seed domain with the tip_element and outer with ndoes of tip_element
  for(int j = 0; j < conn.cols(); ++j)
    outer_nodes.insert(conn(tip_element,j));
  domain_elem.insert(tip_element);

  // Loop for adding layers
  for (int ni = 0; ni < nlayers; ++ni){
    // Insert all the neighbours of outer_nodes to domain_elem
    for (set<int>::iterator it = outer_nodes.begin(); it != outer_nodes.end(); ++it){
      for (int k = 0; k < neighbours[*it].size(); ++k){
        domain_elem.insert(neighbours[*it][k]);
      }
    }

    // Add outer_nodes to the existing inner_nodes
    for (set<int>::iterator it = outer_nodes.begin(); it != outer_nodes.end(); ++it){
      inner_nodes.insert(*it);
    }

    // Add nodes which are not present in inner nodes to outer_nodes
    outer_nodes = set<int>();
    for (set<int>::iterator it = domain_elem.begin(); it != domain_elem.end(); ++it){
      for (int k = 0; k < conn.cols(); ++k){
        if(inner_nodes.count(conn(*it,k)) == 0){
          outer_nodes.insert(conn(*it,k));
        }
      }
    }

  }

  // for (set<int>::iterator it = domain_elem.begin(); it != domain_elem.end(); ++it){
  //   cout << *it << " ";
  // }
  // cout << endl;

  VectorXd q = -1*VectorXd::Ones(nelm*4+nnod);
  // Assign q values to nodes

  // cout << "\nInner Nodes: ";
  for (set<int>::iterator it = inner_nodes.begin(); it != inner_nodes.end(); ++it){
    q[*it] = 1;
    // cout << (*it)+1 << " ";

  }

  // cout << "\nOuter Nodes: ";
  for (set<int>::iterator it = outer_nodes.begin(); it != outer_nodes.end(); ++it){
    q[*it] = 0;
    // cout << (*it)+1 << " ";
  }

  // cout << "\nDomain Elem: ";
  for (set<int>::iterator it = domain_elem.begin(); it != domain_elem.end(); ++it){
    // cout << (*it)+1 << " ";
    if(discont[*it] && discont[*it] != 6){
      vector<vector<int> > lconn = fn_elements[*it].conn;
      bool out = 0;

      for (int j = 0; j < lconn.size(); ++j){
        vector<int> nodes = lconn[j];
        for (int k = 0; k < nodes.size(); ++k){
          if(outer_nodes.count(nodes[k]) != 0)
            out = 1;
        }
      }
      for (int j = 0; j < lconn.size(); ++j){
        vector<int> nodes = lconn[j];
        for (int k = 0; k < nodes.size(); ++k){
          if(nodes[k] >= nnod && !out){
            q(nodes[k]) = 1;

          }
          else if(nodes[k] >= nnod && out){
            q(nodes[k]) = 0;
          }
        }
      }
    }
  }

  vector<double> xpos,ypos;
  for (set<int>::iterator it = domain_elem.begin(); it != domain_elem.end(); ++it){
    if(discont[*it]){
      vector<vector<int> > lconn = fn_elements[*it].conn;
      bool out = 0;

      for (int j = 0; j < lconn.size(); ++j){
        vector<int> nodes = lconn[j];
        for (int k = 0; k < nodes.size(); ++k){
          if(nodes[k] >= nnod){
            q(nodes[k]) = 1;
            xpos.push_back(x(nodes[k],0));
            ypos.push_back(x(nodes[k],1));
          }
        }
      }
    }
  }
  // print(xpos);
  // print(ypos);

  double sy = 0, sx = 0, sxx = 0, sxy = 0;
  int n = xpos.size();
  for (int i = 0; i < n; i++){
    sx += xpos[i];
    sy += ypos[i];
    sxx += xpos[i]*xpos[i];
    sxy += xpos[i]*ypos[i];
  }

  if(abs(n*sxx-sx*sx) < eps){
    direction = make_pair(0,1);
  }
  else{
    double m = (n*sxy-sx*sy)/(n*sxx-sx*sx);
    direction = make_pair(1/sqrt(1+m*m),m/sqrt(1+m*m));
  }
  // cout << n*sxx-sx*sx << endl;
  // cout << sx << " "<< sy << " "<< sxy << " "<< sxx << " "<< sxy << endl;
  // cout << direction.first << " " << direction.second << endl;

  // cout << "\nQ-vector: " << endl;
  // cout << q << endl;
  for (set<int>::iterator it = domain_elem.begin(); it != domain_elem.end(); ++it){
    int i = *it;
    if(!discont[i] || discont[i] == 6){
      VectorXi nodes = conn(i,all);
      MatrixXd xv = x(nodes,all);
      // cout << "Element " << i << endl;
      vector<int> dof = {   nodes(0)*2, nodes(0)*2+1,
                            nodes(1)*2, nodes(1)*2+1,
                            nodes(2)*2, nodes(2)*2+1,
                            nodes(3)*2, nodes(3)*2+1};
      VectorXd u = un1(dof);
      VectorXd lq = q(nodes);
      // cout << "Element " << i << ":\n"<< lq << endl;
      // cout << "Element " << i << ": ";
      // double temp = j_int;
      j_int += quad_j(xv, u, lq, E, nu, direction);
      // cout << "Element " << i << ": "<< j_int-temp << endl;

    }
    else{
      vector<vector<int> > lconn = fn_elements[i].conn;

      for (int j = 0; j < lconn.size(); ++j){
        if(fn_elements[i].active[j]){
          vector<int> nodes = lconn[j];
          MatrixXd xv = x(nodes,all);

          vector<int> dof = { nodes[0]*2, nodes[0]*2+1,
                              nodes[1]*2, nodes[1]*2+1,
                              nodes[2]*2, nodes[2]*2+1};
          VectorXd u = un1(dof);
          VectorXd lq = q(nodes);
          // cout << "Element " << i << "-subelement "<< j << ":\n" << lq << endl;
          // double temp = j_int;
          j_int += tri_j(xv, u, lq, E, nu, direction);
          // cout << "Element " << i << "-subelement "<< j << ": " << j_int-temp << endl;

        }
      }

    }
  }
  return j_int;
}
