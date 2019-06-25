#include "utilities.cpp"

// Note that with 0-indexed arrays, first dof is node*(#dofs), next is node*(#dofs)+1, and so on
void boundary_conditions(VectorXd &vn, VectorXd &vn1){
  double bc;
  // Syntax: LinSpaced(#elements,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXi left = VectorXi::LinSpaced(4,1,13);
  left = left.array()-1;
  bc = -0.01;
  for(int i=0; i < left.size();++i){
    long node_bc = left(i);
    long idof = node_bc*2;
    vn(idof) = bc;
    vn1(idof) = bc;
  }

  VectorXi right = VectorXi::LinSpaced(4,4,16);
  right = right.array()-1;
  bc = 0.01;
  for(int i=0; i < right.size();++i){
    long node_bc = right(i);
    long idof = node_bc*2;
    vn(idof) = bc;
    vn1(idof) = bc;
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
