#include "utilities.cpp"

// Note that with 0-indexed arrays, first dof is node*(#dofs), next is node*(#dofs)+1, and so on
void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  // Syntax: LinSpaced(#elements,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXi bottom = VectorXi::LinSpaced(41,1,41);
  bottom = bottom.array()-1;
  bcx = 0;
  bcy = 0;
  for(int i=0; i < bottom.size();++i){
    long node_bc = bottom(i);
    long idof = node_bc*2;
    vn(idof) = bcx;
    vn1(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof+1) = bcy;
  }

  VectorXi top = VectorXi::LinSpaced(41,862,902);
  top = top.array()-1;
  bcx = 1e-4;
  bcy = 0;
  for(int i=0; i < top.size();++i){
    long node_bc = top(i);
    long idof = node_bc*2;
    vn(idof) = bcx;
    vn1(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof+1) = bcy;
  }
}
