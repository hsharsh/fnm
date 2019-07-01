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

inline double disp_interpolation(double u1, double u2, double e){
  return u1*(1-e) + u2*e;
}

inline double distance(VectorXd x1, VectorXd x2){
  return sqrt(pow(x1(0)-x2(0),2)+pow(x1(1)-x2(1),2));
}
void partition(element &elem, int lmn, MatrixXi &conn, MatrixXd &x, VectorXd &un1, int &ndof){
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
    x(nnod,all) = point_interpolation(xv(start,all),xv((start+1)%4,all),elem.edge[start]);
    x(nnod+1,all) = x(nnod,all);
    un1(ndof) = disp_interpolation(ux(start),ux((start+1)%4),elem.edge[start]);
    un1(ndof+1) = disp_interpolation(uy(start),uy((start+1)%4),elem.edge[start]);
    un1(ndof+2) = un1(ndof);
    un1(ndof+3) = un1(ndof+1);

    x(nnod+2,all) = point_interpolation(xv((start+3)%4,all),xv(start,all),elem.edge[(start+3)%4]);
    x(nnod+3,all) = x(nnod+2,all);

    un1(ndof+4) = disp_interpolation(ux((start+3)%4),ux(start),elem.edge[(start+3)%4]);
    un1(ndof+5) = disp_interpolation(uy((start+3)%4),uy(start),elem.edge[(start+3)%4]);
    un1(ndof+6) = un1(ndof+4);
    un1(ndof+7) = un1(ndof+5);

    elem.conn.push_back({nodes(start), nnod , nnod+3});
    elem.conn.push_back({nnod+1, nodes((start+1)%4), nodes((start+2)%4)});
    elem.conn.push_back({nnod+1, nodes((start+2)%4), nnod+2});
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
    un1(ndof) = disp_interpolation(ux((start+1)%4),ux((start+2)%4),elem.edge[(start+1)%4]);
    un1(ndof+1) = disp_interpolation(uy((start+1)%4),uy((start+2)%4),elem.edge[(start+1)%4]);
    un1(ndof+2) = un1(ndof);
    un1(ndof+3) = un1(ndof+1);

    x(nnod+2,all) = point_interpolation(xv((start+3)%4,all),xv(start,all),elem.edge[(start+3)%4]);
    x(nnod+3,all) = x(nnod+2,all);
    un1(ndof+4) = disp_interpolation(ux((start+3)%4),ux(start),elem.edge[(start+3)%4]);
    un1(ndof+5) = disp_interpolation(uy((start+3)%4),uy(start),elem.edge[(start+3)%4]);
    un1(ndof+6) = un1(ndof+4);
    un1(ndof+7) = un1(ndof+5);

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


void floating_nodes(vector<int> &discont, map <int,element> &fn_elements, MatrixXi &conn, MatrixXd &x, VectorXd &un1, int &ndof){

  int nelm = conn.rows();
  for (int i = 0; i < nelm; ++i){
    if(discont[i]==1){
      partition(fn_elements[i], i, conn, x, un1, ndof);
      discont[i]=2;
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
        double area_elem = 0.5*(xv(0,0)*xv(1,1)-xv(1,0)*xv(0,1) + xv(1,0)*xv(2,1)-xv(2,0)*xv(1,1) + xv(2,0)*xv(0,1)-xv(0,0)*xv(2,1));
        if(area_elem < active_tol){
          it->second.active[j] = 0;
        }
      }
      it->second.processed = 1;
    }
  }
}
