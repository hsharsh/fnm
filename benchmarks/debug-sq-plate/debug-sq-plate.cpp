// Boundary conditions for simple case used for debugging
// Activate only the boundary_conditions function. No temporary bc is needed

void boundary_conditions(VectorXd &vn, VectorXd &vn1){
  double bc;
  // Syntax: LinSpaced(increment,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXd left = VectorXd::LinSpaced(4,1,13);
  left = left.array()-1;
  bc = -0.01;
  for(int i=0; i < left.size();i++){
    long node_bc = left(i);
    long idof = node_bc*2;
    vn(idof) = bc;
    vn1(idof) = bc;
  }

  VectorXd right = VectorXd::LinSpaced(4,4,16);
  right = right.array()-1;
  bc = 0.01;
  for(int i=0; i < right.size();i++){
    long node_bc = right(i);
    long idof = node_bc*2;
    vn(idof) = bc;
    vn1(idof) = bc;
  }
}
