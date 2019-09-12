#include "utilities.cpp"

// Don't forget to make dicsont corresponding to the element as "1" to activate the floating nodes.
void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
  vector <int> cracked;
  for (int i = 401; i <= 408; ++i){
    cracked.push_back(i-1);
  }
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

int stress_based_crack(vector<int> &discont, map <int,element> &fn_elements, map <pair<int,int>,double> &cparam, MatrixXi &conn, MatrixXd &x, VectorXd &un1, int &ndof, double E, double nu){
  int ci = 0;
  int nelm = conn.rows();
  vector<double> xgp = {0};
  vector<double> wgp = {2};
  int ngp = wgp.size();

  for(int i = 0; i < nelm; ++i){
    if(!discont[i]){

      VectorXi nodes = conn(i,all);
      MatrixXd xv = x(nodes,all);
      bool connected = 0;
      for(int j = 0; j < 4; j++){
        if(cparam.find(make_pair(nodes(j),nodes((j+1)%4))) != cparam.end()){
          connected = 1;
          break;
        }
        if(cparam.find(make_pair(nodes((j+1)%4),nodes(j))) != cparam.end()){
          connected = 1;
          break;
        }
      }
      if(!connected)
        continue;

      vector<int> dof = {   nodes(0)*2, nodes(0)*2+1,
                            nodes(1)*2, nodes(1)*2+1,
                            nodes(2)*2, nodes(2)*2+1,
                            nodes(3)*2, nodes(3)*2+1};
      VectorXd u = un1(dof);

      double area_elem = 0.0;
      VectorXd str = VectorXd::Zero(3);

      for (int j = 0; j < ngp; ++j){
        for (int k = 0; k < ngp; ++k){
          double r = xgp[j], s =  xgp[k];
          MatrixXd B0(2,4), jac(2,2), B1(3,4), B2 = MatrixXd::Zero(4,4), B3 = MatrixXd::Zero(4,8);

          B0 << -(1-s)/4,  (1-s)/4, (1+s)/4, -(1+s)/4,
                -(1-r)/4,  -(1+r)/4,  (1+r)/4,  (1-r)/4;

          jac = (B0*xv);

          B1 << 1, 0, 0, 0,
                0, 0, 0, 1,
                0, 0.5, 0.5, 0;

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
                0, 0, (1-2*nu);
          D = E/((1+nu)*(1-2*nu))*D.array();

          str = str.array() + (D*strain).array()*wgp[j]*wgp[k];
          area_elem += jac.determinant()*wgp[j]*wgp[k];
        }
      }

      // double smises = sqrt((pow(str(0)-str(1),2) + pow(str(0),2) + pow(str(1),2) + 6*(pow(str(2),2)))/2)/area_elem;
      double se = (pow(str(0),2)+pow(str(1),2)+pow(str(2),2) -2*nu*(str(0)*str(1)+str(1)*str(2)+str(2)*str(0)))/(pow(area_elem,2));

      // if(max(eig_values(0),eig_values(1)) > sy){
      if(se > sy){
        element fn;

        double cx = (xv(0,0)+xv(1,0)+xv(2,0)+xv(3,0))/4, cy = (xv(0,1)+xv(1,1)+xv(2,1)+xv(3,1))/4;

        double r = 0, s = 0;

        for(int j = 0; j < 4; j++){
          if(cparam.find(make_pair(nodes(j),nodes((j+1)%4))) != cparam.end()){
            double e = cparam[make_pair(nodes(j),nodes((j+1)%4))];
            cx = xv(j,0)*(1-e)+xv((j+1)%4,0)*e;
            cy = xv(j,1)*(1-e)+xv((j+1)%4,1)*e;
            if(j == 0){
              r = e; s = 1;
            }
            else if(j == 1){
              r = -1; s = e;
            }
            else if(j == 2){
              r = e; s = -1;
            }
            else{
              r = 1; s = e;
            }
            break;
          }
          if(cparam.find(make_pair(nodes((j+1)%4),nodes(j))) != cparam.end()){
            double e = 1-cparam[make_pair(nodes((j+1)%4),nodes(j))];
            cx = xv(j,0)*(1-e)+xv((j+1)%4,0)*e;
            cy = xv(j,1)*(1-e)+xv((j+1)%4,1)*e;
            break;
            if(j == 0){
              r = 1-e; s = 1;
            }
            else if(j == 1){
              r = -1; s = 1-e;
            }
            else if(j == 2){
              r = 1-e; s = -1;
            }
            else{
              r = 1; s = 1-e;
            }
            break;
          }
        }

        MatrixXd B0(2,4), jac(2,2), B1(3,4), B2 = MatrixXd::Zero(4,4), B3 = MatrixXd::Zero(4,8);

        B0 << -(1-s)/4,  (1-s)/4, (1+s)/4, -(1+s)/4,
              -(1-r)/4,  -(1+r)/4,  (1+r)/4,  (1-r)/4;

        jac = (B0*xv);

        B1 << 1, 0, 0, 0,
              0, 0, 0, 1,
              0, 0.5, 0.5, 0;

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
              0, 0, (1-2*nu);
        D = E/((1+nu)*(1-2*nu))*D.array();

        str = str.array() + (D*strain).array();

        MatrixXd sig(2,2);
        sig << str(0), str(2),
                str(2), str(1);

        SelfAdjointEigenSolver<MatrixXd> es(sig);
        VectorXd eig_values = es.eigenvalues();

        VectorXd direction = VectorXd::Zero(2);

        if(eig_values(0) > eig_values(1))
          direction = es.eigenvectors().col(0);
        else
          direction = es.eigenvectors().col(1);

        MatrixXd T(2,2);
        T <<  cos(pi/2), -sin(pi/2),
              sin(pi/2), cos(pi/2);

        direction = T*direction;

        double dx = direction(0), dy = direction(1);


        for (int j = 0; j < 4; j++){
          fn.edge[j] = (dx*(cy-xv(j,1))-dy*(cx-xv(j,0)))/(dx*(xv((j+1)%4,1)-xv(j,1))-dy*(xv((j+1)%4,0)-xv(j,0)));
          if(fn.edge[j] > 1.0 || fn.edge[j] < 0.0){
            fn.edge[j] = NAN;
          }
          if(!isnan(fn.edge[j])){
            if(cparam.find(make_pair(nodes(j),nodes((j+1)%4))) == cparam.end() && cparam.find(make_pair(nodes((j+1)%4),nodes(j))) == cparam.end()){
              cparam[make_pair(nodes(j),nodes((j+1)%4))] = fn.edge[j];
            }
          }
        }
        cout << "Crack added at element " << i << endl;
        if(cracked == 0){
          cracked = 1;
        }
        ci = 1;
        discont[i] = 1;
        fn_elements[i] = fn;
      }
    }
  }
  return ci;
}
