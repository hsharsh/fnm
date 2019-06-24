// Boundary conditions for simple case used for debugging
// Activate only the boundary_conditions function. No temporary bc is needed

void boundary_conditions(VectorXd &vn, VectorXd &vn1){
  double bc;
  // Syntax: LinSpaced(#elements,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXi left = VectorXi::LinSpaced(4,1,13);
  left = left.array()-1;
  bc = -0.01;
  for(int i=0; i < left.size();i++){
    long node_bc = left(i);
    long idof = node_bc*2;
    vn(idof) = bc;
    vn1(idof) = bc;
  }

  VectorXi right = VectorXi::LinSpaced(4,4,16);
  right = right.array()-1;
  bc = 0.01;
  for(int i=0; i < right.size();i++){
    long node_bc = right(i);
    long idof = node_bc*2;
    vn(idof) = bc;
    vn1(idof) = bc;
  }
}

void crack_def(vector<int> &discont, map<int,element> &fn_elements){
  vector <int> cracked;
  for (int i = 2; i < 9; i+=3){
    cracked.push_back(i-1);
  }
  for (int i = 0; i < cracked.size(); i++){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {0.5, NAN, 0.5, NAN};
  }
}
