#include "utilities.cpp"

// Don't forget to make dicsont corresponding to the element as "1" to activate the floating nodes.
void crack_def(vector<int> &discont, map<int,element> &fn_elements, MatrixXi &conn, map <pair<int,int>,double> &cparam){
  vector <int> cracked;
  for (int i = 401; i <= 407; ++i){
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
  cracked.push_back(407);
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
  discont[408] = 4;
}

// Computing crack propagation by the formulation for maximum tangential stress formulation given in ABAQUS documentation.

int max_tangential_crack(vector<int> &discont, vector<vector<int> > &neighbours, map <int,element> &fn_elements, map <pair<int,int>,double> &cparam, MatrixXi &conn, MatrixXd &x, pair<double,double> K, double j_integral){
  int ci = 0; // Crack-Initiated: Return paramater indicating whether a crack was added or not

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

        // Compute the crack direction

        // Compute the new crack propagation direction using K1, K2
          double K1 = K.first, K2 = K.second;
          cout << "K1: " << K1 << ", K2: " << K2 << endl;
          double cpd = acos((3*pow(K2,2)+sqrt(pow(K1,4)+8*pow(K1,2)*pow(K2,2)))/(pow(K1,2)+9*pow(K2,2)));
          if(K2 > 0){
            cpd = -cpd;
          }
          cout << "CPD: " << cpd*180/pi << endl;
        // Define dx, dy in terms of the theta+orignal_direction_theta

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
