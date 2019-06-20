#include "utilities.cpp"

// Note that with 0-indexed arrays, first dof is node*(#dofs), next is node*(#dofs)+1, and so on
void boundary_conditions(VectorXd &vn, VectorXd &vn1){
  double bc;
  // Syntax: LinSpaced(increment,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXd left = VectorXd::LinSpaced(21,0,420);
  bc = -0.01;
  for(int i=0; i < left.size();i++){
    int node_bc = left(i);
    int idof = node_bc*2;
    vn(idof) = bc;
    vn1(idof) = bc;
  }

  VectorXd right = VectorXd::LinSpaced(21,20,440);
  bc = 0.01;
  for(int i=0; i < right.size();i++){
    int node_bc = right(i);
    int idof = node_bc*2;
    vn(idof) = bc;
    vn1(idof) = bc;
  }
}

void temporary_bc(VectorXd &vn, VectorXd &vn1){
  double bc;
  VectorXd left = VectorXd::LinSpaced(21,0,420);
  bc = 0.01;
  for(int i=0; i < left.size();i++){
    int node_bc = left(i);
    int idof = node_bc*2;
    vn(idof) = bc;
    vn1(idof) = bc;
  }
}
