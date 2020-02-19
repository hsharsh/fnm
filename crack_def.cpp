#include "utilities.cpp"

// Don't forget to make dicsont corresponding to the element as "1" to activate the floating nodes.
void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
  vector <int> cracked;
  for (int i = 401; i <= 419; ++i){
    cracked.push_back(i-1);
  }
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
  cracked.clear();
  cracked.push_back(419);
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
  discont[420] = 4;
}

int j_based_crack(vector<int> &discont, vector<vector<int> > &neighbours, map <int,element> &fn_elements, map <pair<int,int>,double> &cparam, MatrixXi &conn, MatrixXd &x, VectorXd &un1, int &ndof, double E, double nu, double j_integral, int nlayers){
  int ci = 0; // Crack-Initiated: Return paramater indicating whether a crack was added or not

  // vector<double> xgp = {0};
  // vector<double> wgp = {2};
  int nelm = conn.rows();
  int ngp = wgp.size();

  for(int i = 0; i < nelm; ++i){
    if(discont[i] == 6){

      VectorXi nodes = conn(i,all);
      MatrixXd xv = x(nodes,all);

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
          if(!discont[l] || discont[l] == 6){
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

                // Restricting crack to happen at a maximum of max_deviation angle (in radians).
                // double max_deviation = pi/3;
                // if( (abs(theta) < max_deviation || (abs(theta) < pi && abs(theta) > 2*max_deviation)) && cracked )
                //   continue;
                // cout << "Theta: " << theta*180/pi << endl;

                MatrixXd T(2,2);
                T << cos(theta), -sin(theta),
                    sin(theta), cos(theta);

                // cout << "T" << endl;
                // cout << T << endl;
                MatrixXd strt = T.transpose()*str*T;

                // cout << "r: " << r << " s: "<< s<< endl;
                // cout << "Str" << endl << str << endl;
                // cout << "Str-rotated" << endl << strt(1,1) << endl << endl;
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
        // cout << endl;

        for (int j = 0; j < 4; ++j){
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

        int previous_tip = -1, next_transition = -1, tip_edge = -1;
        // Finish the code
        for (int j = 0; j < 4; ++j){
          vector<int> nn = neighbours[nodes[j]];
          vector<int> nn1 = neighbours[nodes[(j+1)%4]];
          bool correct_tip_edge = 0;
          map <int, int> common;

          if(!isnan(fn.edge[j])){
            correct_tip_edge = 1;
          }

          bool adjacent_to_first_node = 0;
          for (int k = 0; k < nn.size(); ++k){
            if(discont[nn[k]] == 5){
              previous_tip = nn[k];
              adjacent_to_first_node = 1;
            }
          }

          for (int k = 0; k < nn1.size(); ++k){
            if(discont[nn1[k]] == 5){
              if(adjacent_to_first_node)
                correct_tip_edge = 0;
            }
          }

          if(!isnan(fn.edge[j])){
            for (int k = 0; k < nn.size(); ++k)
                common[nn[k]]++;
            for (int k = 0; k < nn1.size(); ++k)
                common[nn1[k]]++;
          }

          for (auto it:common){
            if(it.second == 2 && it.first != i && correct_tip_edge && discont[it.first] == 0){
              next_transition = it.first;
            }
          }

          if(correct_tip_edge){
            tip_edge = j;
          }
        }
        // print(fn.edge);
        // cout << previous_tip << endl;
        // cout << next_transition << endl;
        // cout << tip_edge << endl;


        // Remove the connectivity of the previous tip element and set it for reprocessing.
        discont[previous_tip] = 1;
        fn_elements[previous_tip].conn.clear();
        fn_elements[previous_tip].processed = 0;
        fn_elements[previous_tip].active.clear();
        /*
        If there is a transition element
          Set the next transition element for processing.
          Set the current crack tip for processing and reassign the floating node element properties
          Change the sign of the edge dataset to reflect the crack tip
        Else
          Put a complete crack with no transition and tip elements
        */
        if(next_transition != -1){
          fn.edge[tip_edge] = -fn.edge[tip_edge];
          discont[next_transition] = 4;
          discont[i] = 3;
          fn_elements[i] = fn;
        }
        else{
          discont[i] = 1;
          fn_elements[i] = fn;
        }

        cout << "Crack added at element " << i << endl << endl;
        ci = 1;
        // Global variable - Represents if the given geometry has cracked further from initial state
        if(!cracked)
          cracked = 1;
      }
    }
  }
  return ci;
}
