#include "utilities.cpp"

void compute_properties(VectorXd &stress, vector<int> &discont, map <int,element> &fn_elements, MatrixXi &conn, MatrixXd &x, VectorXd &un1, int &ndof, double E, double nu){
  int nelm = conn.rows();
  int nnod = ndof/2;

  vector<vector<double> > str_v(nnod);

  vector<int> rr = {1, -1, -1, 1};
  vector<int> ss = {1, 1, -1, -1};

  for(int i = 0; i < nelm; ++i){
    if(!discont[i]){
      VectorXi nodes = conn(i,all);
      MatrixXd xv = x(nodes,all);

      vector<int> dof = {   nodes(0)*2, nodes(0)*2+1,
                            nodes(1)*2, nodes(1)*2+1,
                            nodes(2)*2, nodes(2)*2+1,
                            nodes(3)*2, nodes(3)*2+1};
      VectorXd u = un1(dof);

      for (int j = 0; j < rr.size(); ++j){
        double r = rr[j], s =  ss[j];
        MatrixXd B0(2,4), jac(2,2), B1(3,4), B2 = MatrixXd::Zero(4,4), B3 = MatrixXd::Zero(4,8);

        B0 << -(1-s)/4,  (1-s)/4, (1+s)/4, -(1+s)/4,
              -(1-r)/4,  -(1+r)/4,  (1+r)/4,  (1-r)/4;

        jac = (B0*xv);

        B1 << 1, 0, 0, 0,
              0, 0, 0, 1,
              0, 1, 1, 0;

        B2(seq(0,1),seq(0,1)) = jac.inverse();
        B2(seq(2,3),seq(2,3)) = jac.inverse();

        // Define B3
        B3(seq(0,1),seq(0,last,2)) = B0;
        B3(seq(2,3),seq(1,last,2)) = B0;

        MatrixXd B = B1*B2*B3;
        MatrixXd strain = B*u;

        MatrixXd D(3,3);
        D << 1-nu, 0, 0,
              0, 1-nu, 0,
              0, 0, (1-2*nu)/2;
        D = E/((1+nu)*(1-2*nu))*D.array();

        VectorXd str = D*strain;
        str_v[nodes(j)].push_back(sqrt((pow(str(0)-str(1),2) + pow(str(0),2) + pow(str(1),2) + 6*(pow(str(2),2)))/2));
      }
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

          double area_elem = 0.5*(xv(0,0)*xv(1,1)-xv(1,0)*xv(0,1) + xv(1,0)*xv(2,1)-xv(2,0)*xv(1,1) + xv(2,0)*xv(0,1)-xv(0,0)*xv(2,1));
          MatrixXd B0(2,3), B2 = MatrixXd::Zero(4,6), B1(3,4);
          B0 << xv(1,1)-xv(2,1), xv(2,1)-xv(0,1), xv(0,1)-xv(1,1),
                  xv(2,0)-xv(1,0), xv(0,0)-xv(2,0), xv(1,0)-xv(0,0);
          B0 = B0.array()/(2*area_elem);

          B1 << 1, 0, 0, 0,
                0, 0, 0, 1,
                0, 1, 1, 0;

          B2(seq(0,1),seq(0,last,2)) = B0;
          B2(seq(2,3),seq(1,last,2)) = B0;

          MatrixXd B = B1*B2;
          MatrixXd strain = B*u;

          MatrixXd D(3,3);
          D << 1-nu, 0, 0,
                0, 1-nu, 0,
                0, 0, (1-2*nu)/2;
          D = E/((1+nu)*(1-2*nu))*D.array();
          MatrixXd str = D*strain;
          double str_mises = sqrt((pow(str(0)-str(1),2) + pow(str(0),2) + pow(str(1),2) + 6*(pow(str(2),2)))/2);
          for (int k = 0; k < 3; ++k)
            str_v[nodes[k]].push_back(str_mises);
        }
      }
    }
  }
  for (int i = 0; i < nnod; ++i){
    if(str_v[i].size() == 0)
      stress(i) = 0.0;
    else
    stress(i) = accumulate(str_v[i].begin(), str_v[i].end(), 0.0) / str_v[i].size();
  }
}

MatrixXd tri_kl(MatrixXd &xv, double &area, double E, double nu){
  MatrixXd kl = MatrixXd::Zero(3*2,3*2);
  double area_elem = 0.5*(xv(0,0)*xv(1,1)-xv(1,0)*xv(0,1) + xv(1,0)*xv(2,1)-xv(2,0)*xv(1,1) + xv(2,0)*xv(0,1)-xv(0,0)*xv(2,1));
  MatrixXd B0(2,3), B2 = MatrixXd::Zero(4,6), B1(3,4);
  B0 << xv(1,1)-xv(2,1), xv(2,1)-xv(0,1), xv(0,1)-xv(1,1),
          xv(2,0)-xv(1,0), xv(0,0)-xv(2,0), xv(1,0)-xv(0,0);
  B0 = B0.array()/(2*area_elem);

  B1 << 1, 0, 0, 0,
        0, 0, 0, 1,
        0, 1, 1, 0;

  B2(seq(0,1),seq(0,last,2)) = B0;
  B2(seq(2,3),seq(1,last,2)) = B0;

  MatrixXd B = B1*B2;

  MatrixXd D(3,3);
  D << 1-nu, 0, 0,
        0, 1-nu, 0,
        0, 0, (1-2*nu)/2;
  D = E/((1+nu)*(1-2*nu))*D.array();

  kl = kl + B.transpose()*D*B*area_elem;
  area += area_elem;

  return kl;
}

MatrixXd quad_kl(MatrixXd &xv, double &area, double E, double nu){
  MatrixXd kl = MatrixXd::Zero(4*2,4*2);

  for(int i = 0; i < ngp; ++i){
    for(int j = 0; j < ngp; ++j){
        double r = xgp[i], s = xgp[j];
        MatrixXd B0(2,4), jac(2,2), B1(3,4), B2 = MatrixXd::Zero(4,4), B3 = MatrixXd::Zero(4,8);

        B0 << -(1-s)/4,  (1-s)/4, (1+s)/4, -(1+s)/4,
              -(1-r)/4,  -(1+r)/4,  (1+r)/4,  (1-r)/4;

        jac = (B0*xv);

        B1 << 1, 0, 0, 0,
              0, 0, 0, 1,
              0, 1, 1, 0;

        B2(seq(0,1),seq(0,1)) = jac.inverse();
        B2(seq(2,3),seq(2,3)) = jac.inverse();

        // Define B3
        B3(seq(0,1),seq(0,last,2)) = B0;
        B3(seq(2,3),seq(1,last,2)) = B0;

        MatrixXd B = B1*B2*B3;

        MatrixXd D(3,3);
        D << 1-nu, 0, 0,
              0, 1-nu, 0,
              0, 0, (1-2*nu)/2;
        D = E/((1+nu)*(1-2*nu))*D.array();

        kl = kl + B.transpose()*D*B*jac.determinant() * wgp[i] * wgp[j];
        area += jac.determinant()* wgp[i] * wgp[j];
    }
  }
  return kl;
}

MatrixXd tri_ml(MatrixXd &xv, double &area, double rho){
  MatrixXd ml = MatrixXd::Zero(3*2,3*2);

  double area_elem = 0.5*(xv(0,0)*xv(1,1)-xv(1,0)*xv(0,1) + xv(1,0)*xv(2,1)-xv(2,0)*xv(1,1) + xv(2,0)*xv(0,1)-xv(0,0)*xv(2,1));

  ml = rho*(area_elem/3)*(MatrixXd::Identity(3*2,3*2)).array();
  area += area_elem;

  return ml;
}

MatrixXd quad_ml(MatrixXd &xv, double &area, double rho){
  MatrixXd ml = MatrixXd::Zero(4*2,4*2);

  double area_elem = 0;
  for(int i = 0; i < ngp; ++i){
    for(int j = 0; j < ngp; ++j){
        double r = xgp[i], s = xgp[j];
        MatrixXd B0(2,4), jac(2,2);
        B0 << -(1-s)/4,  (1-s)/4, (1+s)/4, -(1+s)/4,
              -(1-r)/4,  -(1+r)/4,  (1+r)/4,  (1-r)/4;
        jac = (B0*xv);

        area_elem = area_elem + jac.determinant()* wgp[i] * wgp[j];
    }
  }

  ml = rho*(area_elem/4)*(MatrixXd::Identity(4*2,4*2)).array();
  area += area_elem;

  return ml;
}

void assemble_mg(VectorXd &mg, MatrixXd &x, MatrixXi &conn, vector<int> &discont, map <int,element> fn_elements, double rho, int ndof){
  double area = 0;
  int nelm = conn.rows();
  for(int i = 0; i < nelm; ++i){
    if(!discont[i]){
      VectorXi nodes = conn(i,all);
      MatrixXd xv = x(nodes,all);
      MatrixXd ml = quad_ml(xv,area,rho);

      vector<int> dof = {   nodes(0)*2, nodes(0)*2+1,
                            nodes(1)*2, nodes(1)*2+1,
                            nodes(2)*2, nodes(2)*2+1,
                            nodes(3)*2, nodes(3)*2+1};

      mg(dof) = mg(dof) + ml.diagonal();
    }
    else{
      vector<vector<int> > lconn = fn_elements[i].conn;
      for (int j = 0; j < lconn.size(); ++j){
        if(fn_elements[i].active[j]){
          vector<int> nodes = lconn[j];
          MatrixXd xv = x(nodes,all);
          MatrixXd ml = tri_ml(xv,area,rho);

          vector<int> dof = { nodes[0]*2, nodes[0]*2+1,
                              nodes[1]*2, nodes[1]*2+1,
                              nodes[2]*2, nodes[2]*2+1};

          mg(dof) = mg(dof) + ml.diagonal();
        }
      }
    }
  }
  for (int i = 0; i < ndof; ++i){
    if(mg(i) == 0)
      mg(i) = 1.0;
  }
  // cout << area << endl;
}

void assemble_fi(VectorXd &fi, VectorXd &un, MatrixXd &x, MatrixXi &conn, vector<int> &discont, map <int,element> fn_elements, double E, double nu){
  double area = 0;
  int nelm = conn.rows();
  for(int i = 0; i < nelm; ++i){
    if(!discont[i]){
      VectorXi nodes = conn(i,all);
      MatrixXd xv = x(nodes,all);
      MatrixXd kl = quad_kl(xv,area,E,nu);

      vector<int> dof = {   nodes(0)*2, nodes(0)*2+1,
                            nodes(1)*2, nodes(1)*2+1,
                            nodes(2)*2, nodes(2)*2+1,
                            nodes(3)*2, nodes(3)*2+1};
      VectorXd u = un(dof);

      fi(dof) = fi(dof) + (kl*u);
    }
    else{

      vector<vector<int> > lconn = fn_elements[i].conn;
      // cout << "lconn" << endl;
      // print_vector(lconn);
      // cout << endl;
      for (int j = 0; j < lconn.size(); ++j){
        if(fn_elements[i].active[j]){
          vector<int> nodes = lconn[j];
          MatrixXd xv = x(nodes,all);
          MatrixXd kl = tri_kl(xv,area,E,nu);

          vector<int> dof = { nodes[0]*2, nodes[0]*2+1,
                              nodes[1]*2, nodes[1]*2+1,
                              nodes[2]*2, nodes[2]*2+1};
          VectorXd u = un(dof);
          // cout << endl;
          // cout << "xv" << endl;
          // cout << xv << endl;
          // cout << "kl" << endl;
          // cout << kl << endl;
          fi(dof) = fi(dof) + (kl*u);

        }
      }
    }
  }
}

void assemble_lcg(VectorXd &lcg, VectorXd &vn, MatrixXd &x, MatrixXi &conn, vector<int> &discont, map <int,element> fn_elements, double eta){
  double area = 0;
  int nelm = conn.rows();
  for(int i = 0; i < nelm; ++i){
    if(!discont[i]){
      VectorXi nodes = conn(i,all);
      MatrixXd xv = x(nodes,all);
      MatrixXd cl = (MatrixXd::Ones(4*2,4*2)).array()*eta;
      vector<int> dof = {   nodes(0)*2, nodes(0)*2+1,
                            nodes(1)*2, nodes(1)*2+1,
                            nodes(2)*2, nodes(2)*2+1,
                            nodes(3)*2, nodes(3)*2+1};
      VectorXd v = vn(dof);

      lcg(dof) = lcg(dof) + (cl*v);
    }
    else{
      vector<vector<int> > lconn = fn_elements[i].conn;
      for (int j = 0; j < lconn.size(); ++j){
        if(fn_elements[i].active[j]){
          VectorXi nodes = conn(i,all);
          MatrixXd xv = x(nodes,all);
          MatrixXd cl = (MatrixXd::Ones(3*2,3*2)).array()*eta;

          vector<int> dof = { nodes[0]*2, nodes[0]*2+1,
                              nodes[1]*2, nodes[1]*2+1,
                              nodes[2]*2, nodes[2]*2+1};
          VectorXd v = vn(dof);

          lcg(dof) = lcg(dof) + (cl*v);

        }
      }
    }
  }
}
