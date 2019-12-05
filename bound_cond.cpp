#include "utilities.cpp"

// Note that with 0-indexed arrays, first dof is node*(#dofs), next is node*(#dofs)+1, and so on
void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  VectorXi right(6);
  right << 9,  10,  11,  12, 112, 170;
  right = right.array()-1;
  bcx = 0;
  bcy = 0;
  for(int i=0; i < right.size();++i){
    int node_bc = right(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof) = bcx;
    vn1(idof+1) = bcy;
  }

  VectorXi leftbottom(3);
  leftbottom << 6,  7, 72;
  leftbottom = leftbottom.array()-1;
  bcx = 0;
  bcy = -0.001;
  for(int i=0; i < leftbottom.size();++i){
    int node_bc = leftbottom(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof) = bcx;
    vn1(idof+1) = bcy;
  }

  VectorXi lefttop(3);
  lefttop << 1,  4, 52;
  lefttop = lefttop.array()-1;
  bcx = 0;
  bcy = 0.001;
  for(int i=0; i < lefttop.size();++i){
    int node_bc = lefttop(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof) = bcx;
    vn1(idof+1) = bcy;
  }
}


void temporary_bc(VectorXd &vn, VectorXd &vn1){
  double bcx, bcy;
  VectorXi left = VectorXi::LinSpaced(5,1,169);
  left = left.array()-1;
  bcx = 0.00866025;
  bcy = 0.005;
  for(int i=0; i < left.size();++i){
    int node_bc = left(i);
    int idof = node_bc*2;
    vn(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof) = bcx;
    vn1(idof+1) = bcy;
  }
}
