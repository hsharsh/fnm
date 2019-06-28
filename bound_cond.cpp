#include "utilities.cpp"

// Note that with 0-indexed arrays, first dof is node*(#dofs), next is node*(#dofs)+1, and so on
void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  // Syntax: LinSpaced(#elements,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXi bottom(38);
  bottom <<     3,   4,   9,  16,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76, 214,
215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 351, 352, 353, 354, 355,
356, 357, 358, 359, 360, 361;
  bottom = bottom.array()-1;
  bcx = -0.01;
  bcy = -0.01;
  for(int i=0; i < bottom.size();++i){
    long node_bc = bottom(i);
    long idof = node_bc*2;
    vn(idof) = bcx;
    vn1(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof+1) = bcy;
  }

  VectorXi top(38);
  top <<   14,  15,  19,  21, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 464,
 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 533, 534, 535, 536, 537,
 538, 539, 540, 541, 542, 543;
  top = top.array()-1;
  bcx = 0.01;
  bcy = 0.01;
  for(int i=0; i < top.size();++i){
    long node_bc = top(i);
    long idof = node_bc*2;
    vn(idof) = bcx;
    vn1(idof) = bcx;
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
