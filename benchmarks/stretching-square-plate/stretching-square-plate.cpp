// Boundary conditions for square plate with velocity on both sides
// Activate only the boundary_conditions function. No temporary bc is needed

void boundary_conditions(VectorXd &vn, VectorXd &vn1){
  double bc;
  // Syntax: LinSpaced(increment,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXd left = VectorXd::LinSpaced(21,1,421);
  left = left.array()-1;
  bc = -0.01;
  for(int i=0; i < left.size();i++){
    long node_bc = left(i);
    long idof = node_bc*2;
    vn(idof) = bc;
    vn1(idof) = bc;
  }

  VectorXd right = VectorXd::LinSpaced(21,21,441);
  right = right.array()-1;
  bc = 0.01;
  for(int i=0; i < right.size();i++){
    long node_bc = right(i);
    long idof = node_bc*2;
    vn(idof) = bc;
    vn1(idof) = bc;
  }
}
