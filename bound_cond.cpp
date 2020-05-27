#include "utilities.cpp"

// Note that with 0-indexed arrays, first dof is node*(#dofs), next is node*(#dofs)+1, and so on
void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  // VectorXi left(9);
  // left << 203,  204,  205, 1384, 1385, 1386, 1408, 1409, 1410;
  VectorXi left(7);
  left << 764,  765,  766, 9741, 9754, 9755, 9756;
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
  // VectorXi right(9);
  // right <<  90,   91,   92, 1180, 1181, 1188, 1189, 1196, 1197;
  VectorXi right(7);
  right <<  231,  232,  233, 4990, 5052, 5053, 5054;
  right = right.array()-1;
  bcy = 0;
  for(int i=0; i < right.size();++i){
    int node_bc = right(i);
    int idof = node_bc*2;
    un(idof+1) = bcy;
    un1(idof+1) = bcy;
  }

  VectorXi top(4);
  top << 1565,  1566,  1567, 16025;
  top = top.array()-1;
  bcy = -1e-4;
  for(int i=0; i < top.size();++i){
    int node_bc = top(i);
    int idof = node_bc*2;
    vn(idof+1) = bcy;
    vn1(idof+1) = bcy;
  }
}
