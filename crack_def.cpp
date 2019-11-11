#include "utilities.cpp"

// Don't forget to make dicsont corresponding to the element as "1" to activate the floating nodes.
void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
  vector<int> cracked;
  for(int i = 11129; i <= 11148; ++i)
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

int j_based_crack(vector<int> &discont, vector<vector<int> > &neighbours, map <int,element> &fn_elements, map <pair<int,int>,double> &cparam, MatrixXi &conn, MatrixXd &x, VectorXd &un1, int &ndof, double E, double nu, double j_integral, int nlayers){
  int ci = 0;
  // vector<double> xgp = {0};
  // vector<double> wgp = {2};
  int nelm = conn.rows();
  int ngp = wgp.size();

  for(int i = 0; i < nelm; ++i){
    if(!discont[i]){

      VectorXi nodes = conn(i,all);
      MatrixXd xv = x(nodes,all);

      // Move to next element if the current one is not connected to a pre-existing criterion
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

      if(j_integral > j_tol){
        element fn;

        // Find the position of the crack-tip
        double cx = (xv(0,0)+xv(1,0)+xv(2,0)+xv(3,0))/4, cy = (xv(0,1)+xv(1,1)+xv(2,1)+xv(3,1))/4;


        for(int j = 0; j < 4; j++){
          if(cparam.find(make_pair(nodes(j),nodes((j+1)%4))) != cparam.end()){
            double e = cparam[make_pair(nodes(j),nodes((j+1)%4))];
            cx = xv(j,0)*(1-e)+xv((j+1)%4,0)*e;
            cy = xv(j,1)*(1-e)+xv((j+1)%4,1)*e;
            break;
          }
          if(cparam.find(make_pair(nodes((j+1)%4),nodes(j))) != cparam.end()){
            double e = 1-cparam[make_pair(nodes((j+1)%4),nodes(j))];
            cx = xv(j,0)*(1-e)+xv((j+1)%4,0)*e;
            cy = xv(j,1)*(1-e)+xv((j+1)%4,1)*e;
            break;
          }
        }

        double dx, dy;


        int tip_element = i;

        // Find the domain elements
        set<int> domain_elem, outer_nodes;

        // Seed domain with the tip_element and outer with ndoes of tip_element
        for(int j = 0; j < conn.cols(); ++j)
          outer_nodes.insert(conn(tip_element,j));
        domain_elem.insert(tip_element);

        // Loop for adding layers
        for (int ni = 0; ni < nlayers; ++ni){
          // Insert all the neighbours of outer_nodes to to domain_elem
          for (set<int>::iterator it = outer_nodes.begin(); it != outer_nodes.end(); ++it){
            for (int k = 0; k < neighbours[*it].size(); ++k){
              domain_elem.insert(neighbours[*it][k]);
            }
          }
        }

        double maxsig = -1;
        for (set<int>::iterator it = domain_elem.begin(); it != domain_elem.end(); ++it){
          int l = *it;
          if(!discont[l]){
            VectorXi nodes = conn(l,all);
            MatrixXd xv = x(nodes,all);
            vector<int> dof = {   nodes(0)*2, nodes(0)*2+1,
                                  nodes(1)*2, nodes(1)*2+1,
                                  nodes(2)*2, nodes(2)*2+1,
                                  nodes(3)*2, nodes(3)*2+1};
            VectorXd u = un1(dof);

            // cout << "Element " << l << endl;
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

                MatrixXd D = constitutive(E, nu);
                MatrixXd stress = D*strain;
                MatrixXd str(2,2);
                str << stress(0), stress(2),
                  stress(2), stress(1);

                double x = 0.25*((1-r)*(1-s)*xv(0,0)+(1+r)*(1-s)*xv(1,0)+(1+r)*(1+s)*xv(2,0)+(1-r)*(1+s)*xv(3,0));
                double y = 0.25*((1-r)*(1-s)*xv(0,1)+(1+r)*(1-s)*xv(1,1)+(1+r)*(1+s)*xv(2,1)+(1-r)*(1+s)*xv(3,1));
                double theta = atan2(y-cy,x-cx);


                // cout << "Theta: " << theta << endl;

                MatrixXd T(2,2);
                T << cos(theta), -sin(theta),
                    sin(theta), cos(theta);

                // cout << "T" << endl;
                // cout << T << endl;
                MatrixXd strt = T.transpose()*str*T;

                // cout << "r: " << r << " s: "<< s<< endl;
                // cout << "Str" << endl << str << endl;
                // cout << "Str-rotated" << endl << strt << endl;
                if(strt(1,1) > maxsig){
                  maxsig = strt(1,1);
                  dx = (x-cx)/sqrt(pow(x-cx,2)+pow(y-cy,2));
                  dy = (y-cy)/sqrt(pow(x-cx,2)+pow(y-cy,2));
                  // cout << xv << endl;
                  // cout << "x-: " << x-cx << " y-: "<< y-cy << endl;
                  // cout << "x: " << x << " y: "<< y<< endl;

                }
              }
            }
          }
        }
        // cout << endl;
        cout << "dx: " << dx << " dy: "<< dy << endl;
        // cout << "cx: " << cx << " cy: "<< cy << endl;
        cout << endl;


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
    if(ci == 1){
      return ci;
    }

  }
  return ci;
}
