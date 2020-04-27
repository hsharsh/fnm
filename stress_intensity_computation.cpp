#include "utilities.cpp"

int compute_du(vector<vector<int> > &neighbours, MatrixXi &conn, MatrixXd &x,  vector<int> &discont, map <int,element> &fn_elements, int nnod, double E, double nu, double dk, int nlayers, VectorXd &du, int mode){

  int nelm = conn.rows();

  int tip_element = -1, crack_tip = -1;
  for(int i = 0; i < nelm; ++i){
    if(discont[i] == 6){
      vector<vector<int> > lconn = fn_elements[i].conn;
      for (int j = 0; j < lconn.size(); ++j){
        vector<int> nodes = lconn[j];
        for (int k = 0; k < nodes.size(); ++k){
          if(nodes[k] >= nnod){
            crack_tip = nodes[k];
          }
        }
      }
    }
    if(discont[i] == 6)
      tip_element = i;
  }

  // cout << "Crack tip element -> " << tip_element << ", " << "Crack tip -> " << crack_tip << endl;
  if(tip_element == -1 || crack_tip == -1){
    cerr << "Crack tip is indeterminate. Stress-intensities will not computed" << endl;
    return 1;
  }

  set<int> domain_elem, outer_nodes, inner_nodes;

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
  double cx = x(crack_tip,0), cy = x(crack_tip,1);
  // cout << cx << " " << cy << endl;

  double ctheta = 0;
  int npoints = 0;
  for (set<int>::iterator it = domain_elem.begin(); it != domain_elem.end(); ++it){
    if(discont[*it]){
      vector<vector<int> > lconn = fn_elements[*it].conn;
      for (int j = 0; j < lconn.size(); ++j){
        vector<int> nodes = lconn[j];
        for (int k = 0; k < nodes.size(); ++k){
          if(nodes[k] >= nnod){
            // cout << x(nodes[k],0) << " " << x(nodes[k],1) << endl;
            if(abs(x(nodes[k],0)-cx) > eps || abs(x(nodes[k],1)-cy) > eps){
              ctheta += atan2(cy-x(nodes[k],1),cx-x(nodes[k],0));
              npoints++;
              // cout << atan2(cy-x(nodes[k],1),cx-x(nodes[k],0)) << " ";
            }
          }
        }
      }
    }
  }
  // cout << endl << npoints << endl;
  ctheta = ctheta/npoints;
  // cout << endl <<  "Angle change: " << ctheta*180/pi << endl;
  // cout << endl;
  for (set<int>::iterator it = domain_elem.begin(); it != domain_elem.end(); ++it){
    if(!discont[*it] || discont[*it] == 6){
      VectorXi nodes = conn(*it,all);
      for (int k = 0; k < nodes.size(); ++k){
        // cout << nodes[k] << " ";
        double dx, dy, r, theta, G, mu;
        dx = x(nodes[k],0)-cx;
        dy = x(nodes[k],1)-cy;
        r = sqrt(dx*dx+dy*dy);
        theta = atan2(dy,dx)-ctheta;

        if(theta > pi)
          theta -= 2*pi;
        else if(theta < -pi)
          theta += 2*pi;

        // cout << dx+cx << " " << dy + cy << " : ";
        // cout << theta*180/pi << " " << atan2(dy,dx)*180/pi << endl;
        G = E/(1+nu); // G doesn't have a 1/2 factor as the strains being used for the computation are "true strains" and "not" engineering strains.
        mu = 3-4*nu; // Plane strain
        double u, v;
        if(mode == 1){
          // du(nodes[k]*2) = (dk/(4*G))*(sqrt(r/(2*pi)))*(2*(mu-1)*cos(theta/2)+2*sin(theta)*sin(theta/2));
          // du(nodes[k]*2+1) = (dk/(4*G))*(sqrt(r/(2*pi)))*(2*(mu+1)*sin(theta/2)-2*sin(theta)*cos(theta/2));
          u = (dk/(G))*(sqrt(r/(2*pi)))*cos(theta/2)*(mu-1+2*sin(theta/2)*sin(theta/2));
          v = (dk/(G))*(sqrt(r/(2*pi)))*sin(theta/2)*(mu+1-2*cos(theta/2)*cos(theta/2));
        }
        else{
          // du(nodes[k]*2) = (dk/(4*G))*(sqrt(r/(2*pi)))*(2*(mu+1)*sin(theta/2)+2*sin(theta)*cos(theta/2));
          // du(nodes[k]*2+1) = (dk/(4*G))*(sqrt(r/(2*pi)))*(-2*(mu-1)*cos(theta/2)+2*sin(theta)*sin(theta/2));
          u = (dk/(G))*(sqrt(r/(2*pi)))*sin(theta/2)*(mu+1+2*cos(theta/2)*cos(theta/2));
          v = -(dk/(G))*(sqrt(r/(2*pi)))*cos(theta/2)*(mu-1-2*sin(theta/2)*sin(theta/2));
        }
        du(nodes[k]*2)    = u*cos(ctheta)-v*sin(ctheta);
        du(nodes[k]*2+1)  = v*cos(ctheta)+u*sin(ctheta);
        // cout << nodes[k] << endl;
        // cout << x(nodes[k],0) << " " << x(nodes[k],1) << " " << du(nodes[k]*2) << " " << du(nodes[k]*2+1) << endl;
      }
    }
    else{
      vector<vector<int> > lconn = fn_elements[*it].conn;
      for (int j = 0; j < lconn.size(); ++j){
        vector<int> nodes = lconn[j];
        double dx = 0, dy = 0, lctheta = 0;

        // cout << "Element " << *it << " Subelement " << j << endl;

        for (int k = 0; k < nodes.size(); ++k){
          dx += x(nodes[k],0)-cx;
          dy += x(nodes[k],1)-cy;
          // cout << "x: " << dx << " " << dy << endl;
        }

        lctheta += (atan2(dy,dx)-ctheta);
        if(lctheta > pi)
          lctheta -= 2*pi;
        else if(lctheta < -pi)
          lctheta += 2*pi;

        // cout << "lcTheta: " << lctheta << endl;
        for (int k = 0; k < nodes.size(); ++k){
          // cout << nodes[k] << " ";
          double dx, dy, r, theta, G, mu;
          dx = x(nodes[k],0)-cx;
          dy = x(nodes[k],1)-cy;
          r = sqrt(dx*dx+dy*dy);
          theta = atan2(dy,dx)-ctheta;

          // cout << "Theta: " << theta << endl;

          if(theta > pi)
            theta -= 2*pi;
          else if(theta < -pi)
            theta += 2*pi;

          // cout << "Theta: " << theta << endl;
          if(lctheta < 0){
            theta = -abs(theta);
          }
          else{
            theta = abs(theta);
          }

          // cout << "Modified Theta: " << theta << endl;

          G = E/(1+nu); // G doesn't have a 1/2 factor as the strains being used for the computation are "true strains" and "not" engineering strains.
          mu = 3-4*nu; // Plane strain
          double u, v;
          if(mode == 1){
            // du(nodes[k]*2) = (dk/(4*G))*(sqrt(r/(2*pi)))*(2*(mu-1)*cos(theta/2)+2*sin(theta)*sin(theta/2));
            // du(nodes[k]*2+1) = (dk/(4*G))*(sqrt(r/(2*pi)))*(2*(mu+1)*sin(theta/2)-2*sin(theta)*cos(theta/2));
            u = (dk/(G))*(sqrt(r/(2*pi)))*cos(theta/2)*(mu-1+2*sin(theta/2)*sin(theta/2));
            v = (dk/(G))*(sqrt(r/(2*pi)))*sin(theta/2)*(mu+1-2*cos(theta/2)*cos(theta/2));
          }
          else{
            // du(nodes[k]*2) = (dk/(4*G))*(sqrt(r/(2*pi)))*(2*(mu+1)*sin(theta/2)+2*sin(theta)*cos(theta/2));
            // du(nodes[k]*2+1) = (dk/(4*G))*(sqrt(r/(2*pi)))*(-2*(mu-1)*cos(theta/2)+2*sin(theta)*sin(theta/2));
            u = (dk/(G))*(sqrt(r/(2*pi)))*sin(theta/2)*(mu+1+2*cos(theta/2)*cos(theta/2));
            v = -(dk/(G))*(sqrt(r/(2*pi)))*cos(theta/2)*(mu-1-2*sin(theta/2)*sin(theta/2));
          }
          du(nodes[k]*2)    = u*cos(ctheta)-v*sin(ctheta);
          du(nodes[k]*2+1)  = v*cos(ctheta)+u*sin(ctheta);
          // cout << nodes[k] << endl;
          // cout << x(nodes[k],0) << " " << x(nodes[k],1) << " " << du(nodes[k]*2) << " " << du(nodes[k]*2+1) << endl;
        }
      }
    }
  }
  return 0;
}

pair<double,double> compute_K(vector<vector<int> > &neighbours, MatrixXi &conn, MatrixXd &x, VectorXd &un1, vector<int> &discont,map <int,element> &fn_elements, map <pair<int,int>,double> &cparam, int nnod, double E, double nu, int nlayers, double dk){

  int nelm = conn.rows();

  double J, J1, J2, dJ1, dJ2, K1, K2, H;
  H = E/(1-nu*nu);

  VectorXd du = VectorXd::Zero(nnod*2+nelm*8);

  // Vertical crack
  // du << -0.001000000000, 0.000162021000, -0.001207420000, 0.000140505000, 0.001207420000, 0.000140505000, 0.001000000000, 0.000162021000, -0.001000000000, 0.000297811000, -0.000586148000, -0.000097815500, 0.000586148000, -0.000097815400, 0.001000000000, 0.000297811000, -0.001000000000, 0.000034872800, -0.000195687000, -0.000012365100, 0.000195688000, -0.000012365100, 0.001000000000, 0.000034872800, -0.001000000000, 0.000000000000, -0.000406759000, -0.000000000000, 0.000406760000, 0.000000000000, 0.001000000000, -0.000000000000, 0.000000000008, -0.000302176000, 0.000000000000, 0.000000000000, -0.001766770000, -0.000745558000, 0.001766770000, -0.000745558000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000;

  // Horizontal crack
  // du << 0.000162021000, -0.001000000000, 0.000297811000, -0.001000000000, 0.000034872800, -0.001000000000, 0.000000000000, -0.001000000000, 0.000140505000, -0.001207420000, -0.000097815500, -0.000586148000, -0.000012365100, -0.000195687000, 0.000000000000, -0.000406759000, 0.000140505000, 0.001207420000, -0.000097815400, 0.000586148000, -0.000012365100, 0.000195688000, -0.000000000000, 0.000406760000, 0.000162021000, 0.001000000000, 0.000297811000, 0.001000000000, 0.000034872800, 0.001000000000, -0.000000000000, 0.001000000000, -0.000302176000, 0.000000000008, 0.000000000000, 0.000000000000, -0.000745558000, 0.001766770000, -0.000745558000, -0.001766770000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000;

  // un1 = du;
  // cout << "J: " << compute_j(neighbours, conn, x, du, discont, fn_elements, cparam, nnod, E, nu, nlayers) << ", Computed J: " << dk*dk/H << endl;

// cout << "J-integral value for final iteration: " << compute_j(neighbours, conn, x, du, discont, fn_elements, cparam, nnod, E, nu, nlayers) << endl << endl;
  J = compute_j(neighbours, conn, x, un1, discont, fn_elements, cparam, nnod, E, nu, nlayers);


  // Computation for K1
  compute_du(neighbours, conn, x, discont, fn_elements, nnod, E, nu, dk, nlayers, du, 1);
  // un1 = du;


  VectorXd u1 = un1+du;
  J1 = compute_j(neighbours, conn, x, u1, discont, fn_elements, cparam, nnod, E, nu, nlayers);
  dJ1 = J1 - J;

  // cout << "du: " << du << endl;
  K1 = (H/2)*(dJ1/dk) - dk/2;
  // cout << "J(K1,0): " << compute_j(neighbours, conn, x, du, discont, fn_elements, cparam, nnod, E, nu, nlayers) << ", Theoretical J(K1,0): " << dk*dk/H << endl << endl;

  // Computation for K2
  compute_du(neighbours, conn, x, discont, fn_elements, nnod, E, nu, dk, nlayers, du, 2);

  VectorXd u2 = un1+du;
  J2 = compute_j(neighbours, conn, x, u2, discont, fn_elements, cparam, nnod, E, nu, nlayers);
  dJ2 = J2 - J;

  K2 = (H/2)*(dJ2/dk) - dk/2;
  // cout << "J(0,K2):" << compute_j(neighbours, conn, x, du, discont, fn_elements, cparam, nnod, E, nu, nlayers) << ", Theoretical J(0,K2): " << dk*dk/H << endl << endl;

  // cout << x << endl;
  // cout << un1.size() << endl;
  // cout << "J: " << J << ", K1: " << K1 << ", K2: " << K2 << ", J(K1,K2):" << (K1*K1+K2*K2)/H << endl;

  return make_pair(K1,K2);
}
