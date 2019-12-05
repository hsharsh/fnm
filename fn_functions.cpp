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
  else if(type == 2){
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

    ndof+=8;
  }
  else if(type == 3){
    // Code for T-crack
    int start = -1;
    for(start = 0; start < 4; ++start){
      if(elem.edge[start] > 1.0)
        break;
    }
    x(nnod,all) = point_interpolation(xv((start+1)%4,all),xv((start+2)%4,all),elem.edge[(start+1)%4]);
    x(nnod+1,all) = x(nnod,all);
    un1(ndof) = disp_interpolation(ux((start+1)%4),ux((start+2)%4),elem.edge[(start+1)%4]);
    un1(ndof+1) = disp_interpolation(uy((start+1)%4),uy((start+2)%4),elem.edge[(start+1)%4]);
    un1(ndof+2) = un1(ndof);
    un1(ndof+3) = un1(ndof+1);

    x(nnod+2,all) = point_interpolation(xv((start+2)%4,all),xv((start+3)%4,all),elem.edge[(start+2)%4]);
    x(nnod+3,all) = x(nnod+2,all);
    un1(ndof+4) = disp_interpolation(ux((start+2)%4),ux((start+3)%4),elem.edge[(start+2)%4]);
    un1(ndof+5) = disp_interpolation(uy((start+2)%4),uy((start+3)%4),elem.edge[(start+2)%4]);
    un1(ndof+6) = un1(ndof+4);
    un1(ndof+7) = un1(ndof+5);

    x(nnod+4,all) = point_interpolation(xv((start+3)%4,all),xv(start,all),elem.edge[(start+3)%4]);
    x(nnod+5,all) = x(nnod+4,all);
    un1(ndof+8) = disp_interpolation(ux((start+3)%4),ux(start),elem.edge[(start+3)%4]);
    un1(ndof+9) = disp_interpolation(uy((start+3)%4),uy(start),elem.edge[(start+3)%4]);
    un1(ndof+10) = un1(ndof+8);
    un1(ndof+11) = un1(ndof+9);

    MatrixXd xf(1,2);
    xf = point_interpolation(xv(start,all),xv((start+1)%4,all),elem.edge[start]-1);
    double x1 = x(nnod+4,0), y1 = x(nnod+4,1), x2 = x(nnod,0), y2 = x(nnod,1), x3 = x(nnod+2,0), y3 = x(nnod+2,1), x4 = xf(0,0), y4 = xf(0,1);
    double e = (x1*y3-x1*y4-x3*y1+x3*y4+x4*y1-x4*y3)/(x1*y3-x1*y4-x2*y3+x2*y4-x3*y1+x3*y2+x4*y1-x4*y2);
    // double uf = disp_interpolation(ux((start)%4),ux((start+1)%4),elem.edge[start]-1);
    // double vf = disp_interpolation(uy((start)%4),uy((start+1)%4),elem.edge[start]-1);

    x(nnod+6,all) = point_interpolation(x(nnod+4,all), x(nnod,all),e);
    x(nnod+7,all) = x(nnod+6,all);
    un1(ndof+12) = disp_interpolation(un1(ndof+8),un1(ndof),e);
    un1(ndof+13) = disp_interpolation(un1(ndof+9),un1(ndof+1),e);
    un1(ndof+14) = un1(ndof+12);
    un1(ndof+15) = un1(ndof+13);

    elem.conn.push_back({nodes(start), nodes((start+1)%4) ,nnod, nnod+5});
    elem.conn.push_back({nnod+7, nnod+1, nodes((start+2)%4), nnod+2});
    elem.conn.push_back({nnod+4, nnod+6, nnod+3, nodes((start+3)%4)});

    elem.active = {1,1,1};

    ndof+=16;
  }
  else if(type == 4){
    // Code for intersecting crack

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
        double area_elem = -1;
        if(lconn.size() == 3)
          area_elem = 0.5*(xv(0,0)*xv(1,1)-xv(1,0)*xv(0,1) + xv(1,0)*xv(2,1)-xv(2,0)*xv(1,1) + xv(2,0)*xv(3,1)-xv(3,0)*xv(2,1) + xv(3,0)*xv(0,1)-xv(0,0)*xv(3,1));
        else
          area_elem = 0.5*(xv(0,0)*xv(1,1)-xv(1,0)*xv(0,1) + xv(1,0)*xv(2,1)-xv(2,0)*xv(1,1) + xv(2,0)*xv(0,1)-xv(0,0)*xv(2,1));

        if(area_elem < ar_tol){
          it->second.active[j] = 0;
        }
      }
      it->second.processed = 1;
    }
  }
}
