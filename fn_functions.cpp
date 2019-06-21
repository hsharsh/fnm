#include "utilities.cpp"

struct element{
  vector<double> edge = {NAN,NAN,NAN,NAN};
  vector<vector<int> > conn;
  vector<bool> active;
};

inline int count_nan(vector<double> input){
  int count = 0;
  for(int i = 0; i < input.size(); i++){
    if (isnan(input[i]))
      count++;
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
  else{
    // Code for T-crack and intersecting cracks
  }
  return type;
}

inline MatrixXd point_interpolation(VectorXd x1, VectorXd x2, double e){
  MatrixXd xv(1,2);
  xv << x1(0)*(1-e) + x2(0)*e, x1(1)*(1-e) + x2(1)*e;
  return xv;
}

inline double distance(VectorXd x1, VectorXd x2){
  return sqrt(pow(x1(0)-x2(0),2)+pow(x1(1)-x2(1),2));
}
void partition(element &elem, int lmn, MatrixXi &conn, MatrixXd &x, int &ndof){
  VectorXi nodes = conn(lmn,all);
  MatrixXd xv = x(nodes,all);
  int nnod = ndof/2;
  int type = det_type(elem.edge);
  if (type == 1){
    int start = 0;
    for(start = 1; start < 4; start++){
      if(!isnan(elem.edge[start-1]) && !isnan(elem.edge[start]))
        break;
    }

    x(nnod,all) = point_interpolation(xv(start,all),xv((start+1)%4,all),elem.edge[start]);
    x(nnod+1,all) = x(nnod,all);

    x(nnod+2,all) = point_interpolation(xv((start+3)%4,all),xv(start,all),elem.edge[(start+3)%4]);
    x(nnod+3,all) = x(nnod+2,all);

    elem.conn.push_back({nodes(start), nnod , nnod+3});
    elem.conn.push_back({nnod+1, nodes((start+1)%4), nodes((start+2)%4)});
    elem.conn.push_back({nnod+1, nodes((start+2)%4), nnod+3});
    elem.conn.push_back({nnod+2, nodes((start+2)%4), nodes((start+3)%4)});

    elem.active = {1,1,1,1};

    ndof+=8;
  }
  else if(type = 2){
    int start = 0;
    if (isnan(elem.edge[1]) && isnan(elem.edge[3])){
      start = 1;
    }
    x(nnod,all) = point_interpolation(xv((start+1)%4,all),xv((start+2)%4,all),elem.edge[(start+1)%4]);
    x(nnod+1,all) = x(nnod,all);

    x(nnod+2,all) = point_interpolation(xv((start+3)%4,all),xv(start,all),elem.edge[(start+3)%4]);
    x(nnod+3,all) = x(nnod+2,all);
    if(distance(xv(start,all),x(nnod,all)) < distance(x(nnod+3,all),xv((start+1)%4,all))){
      elem.conn.push_back({nodes(start), nodes((start+1)%4), nnod});
      elem.conn.push_back({nnod, nnod+3, nodes(start)});
    }
    else{
      elem.conn.push_back({nodes(start), nodes((start+1)%4), nnod+3});
      elem.conn.push_back({nodes((start+1)%4), nnod, nnod+3});
    }

    if(distance(x(nnod+2,all),xv((start+2)%4,all)) < distance(xv((start+3)%4,all),x(nnod+1,all))){
      elem.conn.push_back({nnod+2, nnod+1, nodes((start+2)%4)});
      elem.conn.push_back({nodes((start+2)%4), nodes((start+3)%4), nnod+2});
    }
    else{
      elem.conn.push_back({nnod+2, nnod+1, nodes((start+3)%4)});
      elem.conn.push_back({nnod+1, nodes((start+2)%4), nodes((start+3)%4)});
    }
    elem.active = {1,1,1,1};

    ndof+=8;;
  }
  else{
    // Code for T-crack and intersecting cracks
  }
}


void floating_nodes(vector<int> &discont, map <int,element> &fn_elements, MatrixXi &conn, MatrixXd &x, int &ndof){

  int nelm = conn.rows();
  for (int i = 0; i < nelm; i++){
    if(discont[i]==1){
      partition(fn_elements[i], i, conn, x, ndof);
      discont[i]=2;
    }
  }
}
