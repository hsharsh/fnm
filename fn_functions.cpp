#include "utilities.cpp"

struct element{
  vector<double> edge = {NAN,NAN,NAN,NAN};
  vector<vector<int> > conn;
  bool processed = 0;
  vector<bool> active;
};

inline int count_nan(vector<double> input){
  int count = 0;
  for(int i = 0; i < input.size(); ++i){
    if (isnan(input[i]))
      ++count;
  }
  return count;
}

inline int det_type(vector<double> edge){
  int type = 1;
  if(count_nan(edge) == 2){
    if(isnan(edge[0]) && isnan(edge[2]))
      type = 2;
    else if(isnan(edge[1]) && isnan(edge[3]))
      type = 2;
  }
  else if(count_nan(edge) == 0 & (edge[0] > 1.0 || edge[1] > 1.0 || edge[2] > 1.0 || edge[3] > 1.0))
    type = 3;
  else if(count_nan(edge) == 0)
    type = 4;
  else{
    // Code for T-crack and intersecting cracks
    cout << "Error determining type. Check for inconsistent assignment of edge variables" << endl;
  }
  return type;
}

inline MatrixXd point_interpolation(VectorXd x1, VectorXd x2, double e){
  MatrixXd xv(1,2);
  xv << x1(0)*(1-e) + x2(0)*e, x1(1)*(1-e) + x2(1)*e;
  return xv;
}

// inline MatrixXd point_interpolation(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
//   MatrixXd xv(1,2);
//   cout << x1 << " " << y1 << endl;
//   cout << x2 << " " << y2 << endl;
//   cout << x3 << " " << y3 << endl;
//   cout << x4 << " " << y4 << endl;
//
//   xv << ((x1*x2*y4+x2*x3*y1+x3*x4*y2+x4*x1*y3)-(x1*y2*x4+x2*y3*x1+x3*y4*x2+x4*y1*x3))/((x1*y4+x2*y1+x3*y2+x4*y3)-(x1*y2+x2*y3+x3*y4+x4*y1)),
//               ((y1*x2*y4+y2*x3*y1+y3*x4*y2+y4*x1*y3)-(y1*y2*x4+y2*y3*x1+y3*y4*x2+y4*y1*x3))/((x1*y4+x2*y1+x3*y2+x4*y3)-(x1*y2+x2*y3+x3*y4+x4*y1));
//   xv << (x1*x3*(y2-y4)+x1*x4*(y3-y2)+x2*x3*(y4-y1)+x2*x4*(y1-y3))/((x1-x2)*(y3-y4)+x3*(y2-y1)+x4*(y1-y2)),
//         (x1*y2*y3-x1*y2*y4+x2*y1*(y4-y3)-x3*y1*y4+x3*y2*y4+x4*y3*(y1-y2))/((x1-x2)*(y3-y4)+x3*(y2-y1)+x4*(y1-y2));
//
//   cout << xv << endl;
//   return xv;
// }

inline double disp_interpolation(double u1, double u2, double e){
  return u1*(1-e) + u2*e;
}

// inline double disp_interpolation(double u1, double u2, double u3, double u4, double v1, double v2, double v3, double v4, int direction){
//   if (direction == 1)
//     return (u1*u3*(v2-v4)+u1*u4*(v3-v2)+u2*u3*(v4-v1)+u2*u4*(v1-v3))/((u1-u2)*(v3-v4)+u3*(v2-v1)+u4*(v1-v2));
//   else
//     return (u1*v2*v3-u1*v2*v4+u2*v1*(v4-v3)-u3*v1*v4+u3*v2*v4+u4*v3*(v1-v2))/((u1-u2)*(v3-v4)+u3*(v2-v1)+u4*(v1-v2));
// }

inline double distance(VectorXd x1, VectorXd x2){
  return sqrt(pow(x1(0)-x2(0),2)+pow(x1(1)-x2(1),2));
}
void partition(element &elem, int lmn, MatrixXi &conn, map <pair<int,int>,pair<int,int> > &allocated_nodes, MatrixXd &x, VectorXd &un1, int &ndof){
  VectorXi nodes = conn(lmn,all);
  MatrixXd xv = x(nodes,all);
  MatrixXd ux = un1(nodes.array()*2);
  MatrixXd uy = un1(nodes.array()*2+1);

  int nnod = ndof/2;
  int type = det_type(elem.edge);
  if (type == 1){
    int start = -1;
    for(start = 0; start < 4; ++start){
      if(!isnan(elem.edge[start]) && !isnan(elem.edge[(start+1)%4]))
        break;
    }
    start = (start+1)%4;

    vector<int> old_nodes{start,(start+1)%4,(start+2)%4,(start+3)%4};
    vector<int> fl_nodes{nnod,nnod+1,nnod+2,nnod+3};

    // Code to change the nodes assigned
    if(allocated_nodes.find(make_pair(nodes(old_nodes[0]),nodes(old_nodes[1]))) != allocated_nodes.end()){
      fl_nodes[0] = allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[1]))].first;
      fl_nodes[1] = allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[1]))].second;
    }
    else{
      allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[1]))] = make_pair(fl_nodes[0],fl_nodes[1]);
      allocated_nodes[make_pair(nodes(old_nodes[1]),nodes(old_nodes[0]))] = make_pair(fl_nodes[1],fl_nodes[0]);
    }

    if(allocated_nodes.find(make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))) != allocated_nodes.end()){
      fl_nodes[2] = allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))].first;
      fl_nodes[3] = allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))].second;
    }
    else{
      allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))] = make_pair(fl_nodes[2],fl_nodes[3]);
      allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[3]))] = make_pair(fl_nodes[3],fl_nodes[2]);
    }


    x(fl_nodes[0],all) = point_interpolation(xv(old_nodes[0],all),xv(old_nodes[1],all),elem.edge[old_nodes[0]]-eps);
    x(fl_nodes[1],all) = point_interpolation(xv(old_nodes[0],all),xv(old_nodes[1],all),elem.edge[old_nodes[0]]+eps);//x(fl_nodes[0],all);

    un1(fl_nodes[0]*2) = disp_interpolation(ux(old_nodes[0]),ux(old_nodes[1]),elem.edge[start]);
    un1(fl_nodes[0]*2+1) = disp_interpolation(uy(old_nodes[0]),uy(old_nodes[1]),elem.edge[start]);
    un1(fl_nodes[1]*2) = un1(fl_nodes[0]*2);
    un1(fl_nodes[1]*2+1) = un1(fl_nodes[0]*2+1);

    x(fl_nodes[2],all) = point_interpolation(xv(old_nodes[3],all),xv(old_nodes[0],all),elem.edge[old_nodes[3]]-eps);
    x(fl_nodes[3],all) = point_interpolation(xv(old_nodes[3],all),xv(old_nodes[0],all),elem.edge[old_nodes[3]]+eps);//x(fl_nodes[2],all);

    un1(fl_nodes[2]*2) = disp_interpolation(ux(old_nodes[3]),ux(old_nodes[0]),elem.edge[old_nodes[3]]);
    un1(fl_nodes[2]*2+1) = disp_interpolation(uy(old_nodes[3]),uy(old_nodes[0]),elem.edge[old_nodes[3]]);
    un1(fl_nodes[3]*2) = un1(fl_nodes[2]*2);
    un1(fl_nodes[3]*2+1) = un1(fl_nodes[2]*2+1);

    elem.conn.push_back({nodes(old_nodes[0]), fl_nodes[0] , fl_nodes[3]});
    elem.conn.push_back({fl_nodes[1], nodes(old_nodes[1]), nodes(old_nodes[2])});
    elem.conn.push_back({fl_nodes[1], nodes(old_nodes[2]), fl_nodes[2]});
    elem.conn.push_back({fl_nodes[2], nodes(old_nodes[2]), nodes(old_nodes[3])});

    elem.active = {1,1,1,1};

    ndof+=8;
  }
  else if(type == 2){
    int start = 0;
    if (isnan(elem.edge[1]) && isnan(elem.edge[3])){
      start = 1;
    }

    vector<int> old_nodes{start,(start+1)%4,(start+2)%4,(start+3)%4};
    vector<int> fl_nodes{nnod,nnod+1,nnod+2,nnod+3};

    // Code to change the nodes assigned
    if(allocated_nodes.find(make_pair(nodes(old_nodes[1]),nodes(old_nodes[2]))) != allocated_nodes.end()){
      fl_nodes[0] = allocated_nodes[make_pair(nodes(old_nodes[1]),nodes(old_nodes[2]))].first;
      fl_nodes[1] = allocated_nodes[make_pair(nodes(old_nodes[1]),nodes(old_nodes[2]))].second;
    }
    else{
      allocated_nodes[make_pair(nodes(old_nodes[1]),nodes(old_nodes[2]))] = make_pair(fl_nodes[0],fl_nodes[1]);
      allocated_nodes[make_pair(nodes(old_nodes[2]),nodes(old_nodes[1]))] = make_pair(fl_nodes[1],fl_nodes[0]);
    }

    if(allocated_nodes.find(make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))) != allocated_nodes.end()){
      fl_nodes[2] = allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))].first;
      fl_nodes[3] = allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))].second;
    }
    else{
      allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))] = make_pair(fl_nodes[2],fl_nodes[3]);
      allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[3]))] = make_pair(fl_nodes[3],fl_nodes[2]);
    }


    x(fl_nodes[0],all) = point_interpolation(xv(old_nodes[1],all),xv(old_nodes[2],all),elem.edge[old_nodes[1]]-eps);
    x(fl_nodes[1],all) = point_interpolation(xv(old_nodes[1],all),xv(old_nodes[2],all),elem.edge[old_nodes[1]]+eps);//x(fl_nodes[0],all);
    un1(fl_nodes[0]*2) = disp_interpolation(ux(old_nodes[1]),ux(old_nodes[2]),elem.edge[old_nodes[1]]);
    un1(fl_nodes[0]*2+1) = disp_interpolation(uy(old_nodes[1]),uy(old_nodes[2]),elem.edge[old_nodes[1]]);
    un1(fl_nodes[1]*2) = un1(fl_nodes[0]*2);
    un1(fl_nodes[1]*2+1) = un1(fl_nodes[0]*2+1);

    x(fl_nodes[2],all) = point_interpolation(xv(old_nodes[3],all),xv(old_nodes[0],all),elem.edge[old_nodes[3]]-eps);
    x(fl_nodes[3],all) = point_interpolation(xv(old_nodes[3],all),xv(old_nodes[0],all),elem.edge[old_nodes[3]]+eps);//x(fl_nodes[2],all);
    un1(fl_nodes[2]*2) = disp_interpolation(ux(old_nodes[3]),ux(old_nodes[0]),elem.edge[old_nodes[3]]);
    un1(fl_nodes[2]*2+1) = disp_interpolation(uy((old_nodes[3])%4),uy(old_nodes[0]),elem.edge[old_nodes[3]]);
    un1(fl_nodes[3]*2) = un1(fl_nodes[2]*2);
    un1(fl_nodes[3]*2+1) = un1(fl_nodes[2]*2+1);

    if(distance(xv(old_nodes[0],all),x(fl_nodes[0],all)) < distance(x(fl_nodes[3],all),xv(old_nodes[1],all))){
      elem.conn.push_back({nodes(old_nodes[0]), nodes(old_nodes[1]), fl_nodes[0]});
      elem.conn.push_back({fl_nodes[0], fl_nodes[3], nodes(old_nodes[0])});
    }
    else{
      elem.conn.push_back({nodes(old_nodes[0]), nodes(old_nodes[1]), fl_nodes[3]});
      elem.conn.push_back({nodes(old_nodes[1]), fl_nodes[0], fl_nodes[3]});
    }

    if(distance(x(fl_nodes[2],all),xv(old_nodes[2],all)) < distance(xv(old_nodes[3],all),x(fl_nodes[1],all))){
      elem.conn.push_back({fl_nodes[2], fl_nodes[1], nodes(old_nodes[2])});
      elem.conn.push_back({nodes(old_nodes[2]), nodes(old_nodes[3]), fl_nodes[2]});
    }
    else{
      elem.conn.push_back({fl_nodes[2], fl_nodes[1], nodes(old_nodes[3])});
      elem.conn.push_back({fl_nodes[1], nodes(old_nodes[2]), nodes(old_nodes[3])});
    }
    elem.active = {1,1,1,1};

    ndof+=8;
  }
}

void partition_crack_tip(element &elem, int lmn, MatrixXi &conn, map <pair<int,int>,pair<int,int> > &allocated_nodes, MatrixXd &x, VectorXd &un1, int &ndof){
  VectorXi nodes = conn(lmn,all);
  MatrixXd xv = x(nodes,all);
  MatrixXd ux = un1(nodes.array()*2);
  MatrixXd uy = un1(nodes.array()*2+1);
  int crack_tip_edge = -1;

  for(int i = 0; i < 4; i++){
    if(elem.edge[i] < 0){
      crack_tip_edge = i;
      elem.edge[i] = abs(elem.edge[i]);
      break;
    }
  }
  if (crack_tip_edge == -1){
    cerr << "Crack tip not initialized/updated in the edge-dataset" << endl;
  }
  int nnod = ndof/2;
  int type = det_type(elem.edge);
  if (type == 1){
    int start = -1;
    for(start = 0; start < 4; ++start){
      if(!isnan(elem.edge[start]) && !isnan(elem.edge[(start+1)%4]))
        break;
    }
    start = (start+1)%4;

    vector<int> old_nodes{start,(start+1)%4,(start+2)%4,(start+3)%4};
    vector<int> fl_nodes{nnod,nnod+1,nnod+2,nnod+3};

    // Code to change the nodes assigned
    if(allocated_nodes.find(make_pair(nodes(old_nodes[0]),nodes(old_nodes[1]))) != allocated_nodes.end()){
      fl_nodes[0] = allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[1]))].first;
      fl_nodes[1] = allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[1]))].second;
    }
    else{
      allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[1]))] = make_pair(fl_nodes[0],fl_nodes[1]);
      allocated_nodes[make_pair(nodes(old_nodes[1]),nodes(old_nodes[0]))] = make_pair(fl_nodes[1],fl_nodes[0]);
    }

    if(allocated_nodes.find(make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))) != allocated_nodes.end()){
      fl_nodes[2] = allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))].first;
      fl_nodes[3] = allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))].second;
    }
    else{
      allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))] = make_pair(fl_nodes[2],fl_nodes[3]);
      allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[3]))] = make_pair(fl_nodes[3],fl_nodes[2]);
    }


    x(fl_nodes[0],all) = point_interpolation(xv(old_nodes[0],all),xv(old_nodes[1],all),elem.edge[old_nodes[0]]-eps);
    x(fl_nodes[1],all) = point_interpolation(xv(old_nodes[0],all),xv(old_nodes[1],all),elem.edge[old_nodes[0]]+eps);//x(fl_nodes[0],all);

    un1(fl_nodes[0]*2) = disp_interpolation(ux(old_nodes[0]),ux(old_nodes[1]),elem.edge[start]);
    un1(fl_nodes[0]*2+1) = disp_interpolation(uy(old_nodes[0]),uy(old_nodes[1]),elem.edge[start]);
    un1(fl_nodes[1]*2) = un1(fl_nodes[0]*2);
    un1(fl_nodes[1]*2+1) = un1(fl_nodes[0]*2+1);

    x(fl_nodes[2],all) = point_interpolation(xv(old_nodes[3],all),xv(old_nodes[0],all),elem.edge[old_nodes[3]]-eps);
    x(fl_nodes[3],all) = point_interpolation(xv(old_nodes[3],all),xv(old_nodes[0],all),elem.edge[old_nodes[3]]+eps);//x(fl_nodes[2],all);

    un1(fl_nodes[2]*2) = disp_interpolation(ux(old_nodes[3]),ux(old_nodes[0]),elem.edge[old_nodes[3]]);
    un1(fl_nodes[2]*2+1) = disp_interpolation(uy(old_nodes[3]),uy(old_nodes[0]),elem.edge[old_nodes[3]]);
    un1(fl_nodes[3]*2) = un1(fl_nodes[2]*2);
    un1(fl_nodes[3]*2+1) = un1(fl_nodes[2]*2+1);

    if(start == crack_tip_edge){
      elem.conn.push_back({nodes(old_nodes[0]), fl_nodes[0] , fl_nodes[3]});
      elem.conn.push_back({fl_nodes[0], nodes(old_nodes[1]), nodes(old_nodes[2])});
      elem.conn.push_back({fl_nodes[0], nodes(old_nodes[2]), fl_nodes[2]});
      elem.conn.push_back({fl_nodes[2], nodes(old_nodes[2]), nodes(old_nodes[3])});
      allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[1]))].second = -allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[1]))].second;
      allocated_nodes[make_pair(nodes(old_nodes[1]),nodes(old_nodes[0]))].first = -allocated_nodes[make_pair(nodes(old_nodes[1]),nodes(old_nodes[0]))].first;
    }
    else{
      elem.conn.push_back({nodes(old_nodes[0]), fl_nodes[0] , fl_nodes[2]});
      elem.conn.push_back({fl_nodes[1], nodes(old_nodes[1]), nodes(old_nodes[2])});
      elem.conn.push_back({fl_nodes[1], nodes(old_nodes[2]), fl_nodes[2]});
      elem.conn.push_back({fl_nodes[2], nodes(old_nodes[2]), nodes(old_nodes[3])});

      allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))].second = -allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))].second;
      allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[3]))].first = -allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[3]))].first;
    }

    elem.active = {1,1,1,1};

    ndof+=8;
  }
  else if(type == 2){
    int start = 0;
    if (isnan(elem.edge[1]) && isnan(elem.edge[3])){
      start = 1;
    }

    vector<int> old_nodes{start,(start+1)%4,(start+2)%4,(start+3)%4};
    vector<int> fl_nodes{nnod,nnod+1,nnod+2,nnod+3};

    // Code to change the nodes assigned
    if(allocated_nodes.find(make_pair(nodes(old_nodes[1]),nodes(old_nodes[2]))) != allocated_nodes.end()){
      fl_nodes[0] = allocated_nodes[make_pair(nodes(old_nodes[1]),nodes(old_nodes[2]))].first;
      fl_nodes[1] = allocated_nodes[make_pair(nodes(old_nodes[1]),nodes(old_nodes[2]))].second;
    }
    else{
      allocated_nodes[make_pair(nodes(old_nodes[1]),nodes(old_nodes[2]))] = make_pair(fl_nodes[0],fl_nodes[1]);
      allocated_nodes[make_pair(nodes(old_nodes[2]),nodes(old_nodes[1]))] = make_pair(fl_nodes[1],fl_nodes[0]);
    }

    if(allocated_nodes.find(make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))) != allocated_nodes.end()){
      fl_nodes[2] = allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))].first;
      fl_nodes[3] = allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))].second;
    }
    else{
      allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))] = make_pair(fl_nodes[2],fl_nodes[3]);
      allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[3]))] = make_pair(fl_nodes[3],fl_nodes[2]);
    }


    x(fl_nodes[0],all) = point_interpolation(xv(old_nodes[1],all),xv(old_nodes[2],all),elem.edge[old_nodes[1]]-eps);
    x(fl_nodes[1],all) = point_interpolation(xv(old_nodes[1],all),xv(old_nodes[2],all),elem.edge[old_nodes[1]]+eps);//x(fl_nodes[0],all);
    un1(fl_nodes[0]*2) = disp_interpolation(ux(old_nodes[1]),ux(old_nodes[2]),elem.edge[old_nodes[1]]);
    un1(fl_nodes[0]*2+1) = disp_interpolation(uy(old_nodes[1]),uy(old_nodes[2]),elem.edge[old_nodes[1]]);
    un1(fl_nodes[1]*2) = un1(fl_nodes[0]*2);
    un1(fl_nodes[1]*2+1) = un1(fl_nodes[0]*2+1);

    x(fl_nodes[2],all) = point_interpolation(xv(old_nodes[3],all),xv(old_nodes[0],all),elem.edge[old_nodes[3]]-eps);
    x(fl_nodes[3],all) = point_interpolation(xv(old_nodes[3],all),xv(old_nodes[0],all),elem.edge[old_nodes[3]]+eps);//x(fl_nodes[2],all);
    un1(fl_nodes[2]*2) = disp_interpolation(ux(old_nodes[3]),ux(old_nodes[0]),elem.edge[old_nodes[3]]);
    un1(fl_nodes[2]*2+1) = disp_interpolation(uy((old_nodes[3])%4),uy(old_nodes[0]),elem.edge[old_nodes[3]]);
    un1(fl_nodes[3]*2) = un1(fl_nodes[2]*2);
    un1(fl_nodes[3]*2+1) = un1(fl_nodes[2]*2+1);

    if((start+1)%4 == crack_tip_edge){
      if(distance(xv(old_nodes[0],all),x(fl_nodes[0],all)) < distance(x(fl_nodes[3],all),xv(old_nodes[1],all))){
        elem.conn.push_back({nodes(old_nodes[0]), nodes(old_nodes[1]), fl_nodes[0]});
        elem.conn.push_back({fl_nodes[0], fl_nodes[3], nodes(old_nodes[0])});
      }
      else{
        elem.conn.push_back({nodes(old_nodes[0]), nodes(old_nodes[1]), fl_nodes[3]});
        elem.conn.push_back({nodes(old_nodes[1]), fl_nodes[0], fl_nodes[3]});
      }

      if(distance(x(fl_nodes[2],all),xv(old_nodes[2],all)) < distance(xv(old_nodes[3],all),x(fl_nodes[0],all))){
        elem.conn.push_back({fl_nodes[2], fl_nodes[0], nodes(old_nodes[2])});
        elem.conn.push_back({nodes(old_nodes[2]), nodes(old_nodes[3]), fl_nodes[2]});
      }
      else{
        elem.conn.push_back({fl_nodes[2], fl_nodes[0], nodes(old_nodes[3])});
        elem.conn.push_back({fl_nodes[0], nodes(old_nodes[2]), nodes(old_nodes[3])});
      }
      allocated_nodes[make_pair(nodes(old_nodes[1]),nodes(old_nodes[2]))].second = -allocated_nodes[make_pair(nodes(old_nodes[1]),nodes(old_nodes[2]))].second;
      allocated_nodes[make_pair(nodes(old_nodes[2]),nodes(old_nodes[1]))].first = -allocated_nodes[make_pair(nodes(old_nodes[2]),nodes(old_nodes[1]))].first;
    }
    else{
      if(distance(xv(old_nodes[0],all),x(fl_nodes[0],all)) < distance(x(fl_nodes[2],all),xv(old_nodes[1],all))){
        elem.conn.push_back({nodes(old_nodes[0]), nodes(old_nodes[1]), fl_nodes[0]});
        elem.conn.push_back({fl_nodes[0], fl_nodes[2], nodes(old_nodes[0])});
      }
      else{
        elem.conn.push_back({nodes(old_nodes[0]), nodes(old_nodes[1]), fl_nodes[2]});
        elem.conn.push_back({nodes(old_nodes[1]), fl_nodes[0], fl_nodes[2]});
      }

      if(distance(x(fl_nodes[2],all),xv(old_nodes[2],all)) < distance(xv(old_nodes[3],all),x(fl_nodes[1],all))){
        elem.conn.push_back({fl_nodes[2], fl_nodes[1], nodes(old_nodes[2])});
        elem.conn.push_back({nodes(old_nodes[2]), nodes(old_nodes[3]), fl_nodes[2]});
      }
      else{
        elem.conn.push_back({fl_nodes[2], fl_nodes[1], nodes(old_nodes[3])});
        elem.conn.push_back({fl_nodes[1], nodes(old_nodes[2]), nodes(old_nodes[3])});
      }
      allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))].second = -allocated_nodes[make_pair(nodes(old_nodes[3]),nodes(old_nodes[0]))].second;
      allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[3]))].first = -allocated_nodes[make_pair(nodes(old_nodes[0]),nodes(old_nodes[3]))].first;
    }

    elem.active = {1,1,1,1};

    ndof+=8;
  }
}

void partition_transition(element &elem, int lmn, MatrixXi &conn, map <pair<int,int>,pair<int,int> > &allocated_nodes, MatrixXd &x, VectorXd &un1, int &ndof){
  VectorXi nodes = conn(lmn,all);
  int start = -1, floating_node = -1;
  for(start = 0; start < 4; ++start){
    if(allocated_nodes.find(make_pair(nodes(start),nodes((start+1)%4))) != allocated_nodes.end()){
      if(allocated_nodes[make_pair(nodes(start),nodes((start+1)%4))].first > 0)
        floating_node = allocated_nodes[make_pair(nodes(start),nodes((start+1)%4))].first;
      else
        floating_node = allocated_nodes[make_pair(nodes(start),nodes((start+1)%4))].second;
      allocated_nodes[make_pair(nodes(start),nodes((start+1)%4))].first = abs(allocated_nodes[make_pair(nodes(start),nodes((start+1)%4))].first);
      allocated_nodes[make_pair(nodes(start),nodes((start+1)%4))].second = abs(allocated_nodes[make_pair(nodes(start),nodes((start+1)%4))].second);
      allocated_nodes[make_pair(nodes((start+1)%4),nodes(start))].first = abs(allocated_nodes[make_pair(nodes((start+1)%4),nodes(start))].first);
      allocated_nodes[make_pair(nodes((start+1)%4),nodes(start))].second = abs(allocated_nodes[make_pair(nodes((start+1)%4),nodes(start))].second);
      break;
    }
  }
  vector <int> old_nodes = {nodes(start), nodes((start+1)%4), nodes((start+2)%4), nodes((start+3)%4)};

  elem.conn.push_back({old_nodes[0],floating_node,old_nodes[3]});
  elem.conn.push_back({old_nodes[3],floating_node,old_nodes[2]});
  elem.conn.push_back({old_nodes[2],floating_node,old_nodes[1]});

  elem.active = {1,1,1,1};
}

void floating_nodes(vector<int> &discont, map <int,element> &fn_elements, map <pair<int,int>,pair<int,int> > &allocated_nodes, MatrixXi &conn, MatrixXd &x, VectorXd &un1, int &ndof){

  int nelm = conn.rows();
  for (int i = 0; i < nelm; ++i){
    if(discont[i] == 1){
      partition(fn_elements[i], i, conn, allocated_nodes, x, un1, ndof);
      discont[i]=2;
    }
  }

  for(int i = 0; i < nelm; ++i){
    if(discont[i] == 3){
      partition_crack_tip(fn_elements[i], i, conn, allocated_nodes, x, un1, ndof);
      discont[i] = 5;
    }
  }

  for(int i = 0; i < nelm; ++i){
    if(discont[i] == 4){
      partition_transition(fn_elements[i], i, conn, allocated_nodes, x, un1, ndof);
      discont[i] = 6;
    }
  }
}

void remove_singular_elements(map <int,element> &fn_elements, MatrixXd &x){
  for(map<int,element>::iterator it = fn_elements.begin(); it!=fn_elements.end(); ++it){
    if(!it->second.processed){
      vector<vector<int> > lconn = it->second.conn;
      for (int j = 0; j < lconn.size(); ++j){
        vector<int> nodes = lconn[j];
        MatrixXd xv = x(nodes,all);
        double area_elem = -1;

        area_elem = 0.5*(xv(0,0)*xv(1,1)-xv(1,0)*xv(0,1) + xv(1,0)*xv(2,1)-xv(2,0)*xv(1,1) + xv(2,0)*xv(0,1)-xv(0,0)*xv(2,1));

        if(area_elem < ar_tol){
          it->second.active[j] = 0;
        }
      }
      it->second.processed = 1;
    }
  }
}
