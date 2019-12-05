#include "utilities.cpp"

void setcze(vector<pair<int,int> > &cze){
  for(int i = 0; i < cze.size(); ++i){
    cze[i].first = 0;
    cze[i].second = -1;
  }

  VectorXi cohesive = VectorXi::LinSpaced(20,181,200);
  cohesive = cohesive.array()-1;

  for(int i=0; i < cohesive.size();++i){
    cze[cohesive[i]].first = 1;
    cze[cohesive[i]].second = 1;
  }

  VectorXi remove = VectorXi::LinSpaced(20,81,100);
  remove = remove.array()-1;

  for(int i=0; i < remove.size();++i){
    cze[remove[i]].first = 2;
    cze[remove[i]].second = -1;
  }
}

VectorXd fi_cze(MatrixXd &xv, VectorXd &u, int orientation){
  VectorXd fi = VectorXd::Zero(4*2);
  bool fail = 0;
  if (orientation == 1){
    double l = (sqrt(pow(u(2)+xv(1,0)-u(0)-xv(0,0),2)+pow(u(3)+xv(1,1)-u(1)-xv(0,1),2))+sqrt(pow(u(4)+xv(2,0)-u(6)-xv(3,0),2)+pow(u(5)+xv(2,1)-u(7)-xv(3,1),2)))/2;
    double theta = (atan2(u(3)+xv(1,1)-u(1)-xv(0,1),u(2)+xv(1,0)-u(0)-xv(0,0))+atan2(u(5)+xv(2,1)-u(7)-xv(3,1),u(4)+xv(2,0)-u(6)-xv(3,0)))/2;
    // cout << "Theta: " << theta << endl;
    MatrixXd Q(2,2), A(2,4);
    Q << cos(theta), -sin(theta),
        sin(theta), cos(theta);

    A << -1, 0, 1, 0,
          0, -1, 0, 1;

    for(int i = 0; i < ngp; ++i){
      MatrixXd N = MatrixXd::Zero(4,8), T(2,1), U(2,1);
      double r = xgp[i];
      N << (1-r)/2,  0,  (1+r)/2,  0,   0,  0,  0,  0,
            0,  (1-r)/2,  0,  (1+r)/2,  0,  0,  0,  0,
            0,  0,  0,  0,  (1-r)/2,  0,  (1+r)/2,  0,
            0,  0,  0,  0,  0,  (1-r)/2,  0,  (1+r)/2;

      U = Q*A*N*u;

      double us = U(0,0), un = U(1,0);

      cout << "Un: " << un << " Us: " << us << endl;
      if(un < 0){
        T << 0,
              0;
      }
      else if(un < delta){
        T << (27/4)*smax*ashear*(us/delta)*(1-2*(un/delta)+pow(un/delta,2)),
        (27/4)*smax*((un/delta)*(1-2*(un/delta)+pow(un/delta,2))+ashear*pow(us/delta,2)*((un/delta)-1));
        // T << (27/4)*smax*ashear*(us/delta)*(1-2*abs(un/delta)+pow(un/delta,2)),
        //     (27/4)*smax*(abs(un/delta)*(1-2*abs(un/delta)+pow(un/delta,2))+sgn(un)*ashear*pow(us/delta,2)*(abs(un/delta)-1));
      }
      else{
        fail = 1;
      }

      fi = fi + (N.transpose()*A.transpose()*Q.transpose()*T)*(1/l)*wgp[i];
    }
  }
  else{
    double l = (sqrt(pow(u(2)+xv(1,0)-u(4)-xv(2,0),2)+pow(u(3)+xv(1,1)-u(5)-xv(2,1),2))+sqrt(pow(u(0)+xv(0,0)-u(6)-xv(3,0),2)+pow(u(1)+xv(0,1)-u(7)-xv(3,1),2)))/2;
    double theta = (atan2(u(3)+xv(1,1)-u(5)-xv(2,1),u(2)+xv(1,0)-u(4)-xv(2,0))+atan2(u(1)+xv(0,1)-u(7)-xv(3,1),u(0)+xv(0,0)-u(6)-xv(3,0)))/2;
    // cout << "Theta: " << theta << endl;
    MatrixXd Q(2,2), A(2,4);
    Q << cos(theta), -sin(theta),
        sin(theta), cos(theta);

    A << -1, 0, 1, 0,
          0, -1, 0, 1;

    for(int i = 0; i < ngp; ++i){
      MatrixXd N = MatrixXd::Zero(4,8), T(2,1), U(2,1);
      double r = xgp[i];
      N <<  0,  0, (1-r)/2,  0, (1+r)/2,  0,  0,  0,
            0,  0,  0,  (1-r)/2,  0,  (1+r)/2,  0,  0,
            (1-r)/2,  0,  0,  0,  0,  0,  (1+r)/2,  0,
            0,  (1-r)/2,  0,  0,  0,  0,  0,  (1+r)/2;

      U = Q*A*N*u;

      double un = U(0,0), us = U(1,0);

      if(un < 0){
        T << 0,
              0;
      }
      else if(un < delta){
        T << (27/4)*smax*ashear*(us/delta)*(1-2*(un/delta)+pow(un/delta,2)),
            (27/4)*smax*((un/delta)*(1-2*(un/delta)+pow(un/delta,2))+ashear*pow(us/delta,2)*((un/delta)-1));
      }
      else{
        fail = 1;
      }

      fi = fi + (N.transpose()*A.transpose()*Q.transpose()*T)*(1/l)*wgp[i];
    }
  }
  cout << fi << endl;
  if(fail){
    fi << NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN;
    return fi;
  }
  return fi;
}
