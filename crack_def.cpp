#include "utilities.cpp"

void crack_def(vector<int> &discont, map<int,element> &fn_elements){
  vector <int> cracked;
  for (int i = 2; i < 9; i+=3){
    cracked.push_back(i-1);
  }
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {0.5, NAN, 0.5, NAN};
  }
}

void stress_based_crack(vector<int> &discont, map <int,element> &fn_elements, MatrixXi &conn, MatrixXd &x, VectorXd &un1, int &ndof, double E, double nu){
  int nelm = conn.rows();
  vector <double> xgp = {0};
  vector <double> wgp = {2};
  int ngp = wgp.size();

  for(int i = 0; i < nelm; ++i){
    if(!discont[i]){
      VectorXi nodes = conn(i,all);
      MatrixXd xv = x(nodes,all);

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

          str = str.array() + (D*strain).array()*wgp[i]*wgp[j];
          area_elem += jac.determinant()*wgp[i]*wgp[j];
        }
      }

      double stress = sqrt((pow(str(0)-str(1),2) + pow(str(0),2) + pow(str(1),2) + 6*(pow(str(2),2)))/2)/area_elem;

      if(max(str(0),str(1)) > sy){
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

        element fn;

        double cx = (xv(0,0)+xv(1,0)+xv(2,0)+xv(3,0))/4, cy = (xv(0,1)+xv(1,1)+xv(2,1)+xv(3,1))/4;
        for (int j = 0; j < 4; j++){
          fn.edge[j] = (dx*(cy-xv(j,1))-dy*(cx-xv(j,0)))/(dx*(xv((j+1)%4,1)-xv(j,1))-dy*(xv((j+1)%4,0)-xv(j,0)));
          cout << fn.edge[j] << endl;
          if(fn.edge[j] > 1.0 || fn.edge[j] < 0.0){
            fn.edge[j] = NAN;
          }
        }

        discont[i] = 1;
        fn_elements[i] = fn;
      }
    }
  }
}
