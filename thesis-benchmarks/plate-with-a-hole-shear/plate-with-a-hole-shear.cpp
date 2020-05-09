// Boundary conditions for plate with a hole and streching on top and bottom
// Activate only the boundary_conditions function. No temporary bc is needed

void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  // Syntax: LinSpaced(#elements,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXi bottom(31);
  bottom << 3,   4,   9,  16,  63,  64,  65,  66,  67,  68,  69,  70, 195, 196, 197, 198,
199, 200, 201, 202, 203, 204, 205, 319, 320, 321, 322, 323, 324, 325, 326;
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

  VectorXi top(31);
  top <<   14,  15,  19,  21, 278, 279, 280, 281, 282, 283, 284, 285, 422, 423, 424, 425,
 426, 427, 428, 429, 430, 431, 432, 485, 486, 487, 488, 489, 490, 491, 492;
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

tmax = 10
dt = 0.005
E = 1
nu = 0
rho = 1
alpha = 0
sy = 6e8
ar_tol = 1e-8
rf = 1
tc = 0.1
srate = 20
nlyrs = 3
init_c = 0
j_tol = 100
