#include "utilities.cpp"

// Note that with 0-indexed arrays, first dof is node*(#dofs), next is node*(#dofs)+1, and so on
void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  VectorXi left(4);
  left << 203,  204,  205, 1409;
  left = left.array()-1;
  bcx = 0;
  bcy = 0;
  for(int i=0; i < left.size();++i){
    int node_bc = left(i);
    int idof = node_bc*2;
    un(idof) = bcx;
    un(idof+1) = bcy;
    un1(idof) = bcx;
    un1(idof+1) = bcy;
  }
  VectorXi right(4);
  right <<  90,   91,   92, 1189;
  right = right.array()-1;
  bcx = 0;
  bcy = 0;
  for(int i=0; i < right.size();++i){
    int node_bc = right(i);
    int idof = node_bc*2;
    un(idof) = bcx;
    un(idof+1) = bcy;
    un1(idof) = bcx;
    un1(idof+1) = bcy;
  }

  VectorXi top(6);
  top << 495,  496,  497,  498, 2868, 2869;
  top = top.array()-1;
  bcy = -1e-3;
  for(int i=0; i < top.size();++i){
    int node_bc = top(i);
    int idof = node_bc*2;
    vn(idof+1) = bcy;
    vn1(idof+1) = bcy;
  }
}

void temporary_bc(VectorXd &vn, VectorXd &vn1){
  double bcx, bcy;
  VectorXi left = VectorXi::LinSpaced(5,1,169);
  left = left.array()-1;
  bcx = 0.00866025;
  bcy = 0.006;
  for(int i=0; i < left.size();++i){
    int node_bc = left(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof) = bcx;
    vn1(idof+1) = bcy;
  }
}
