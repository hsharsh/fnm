// Boundary conditions for plate with a hole and streching on top and bottom
// Activate only the boundary_conditions function. No temporary bc is needed
void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  // Syntax: LinSpaced(#elements,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXi bottom(21);
  bottom <<    6,   7, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115,
              116, 117, 118, 119, 120;
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

  VectorXi top(21);
  top << 2,  3, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
        41, 42, 43, 44, 45;
  top = top.array()-1;
  bcx = 0;
  bcy = 1e-4;
  for(int i=0; i < top.size();++i){
    long node_bc = top(i);
    long idof = node_bc*2;
    vn(idof) = bcx;
    vn1(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof+1) = bcy;
  }

  VectorXi left(42);
  left << 3,   5,   6,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,
        59,  60,  61,  62,  63,  64,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,
        93,  94,  95,  96,  97,  98,  99, 100, 101, 869;
  left = left.array()-1;
  bcx = 0;
  for(int i=0; i < left.size();++i){
    long node_bc = left(i);
    long idof = node_bc*2;
    un(idof) = bcx;
    un1(idof) = bcx;
  }
}

tmax = 10
dt = 0.02
E = 1
nu = 0
rho = 1
alpha = 0
sy = 6e8
ar_tol = 1e-8
rf = 2
tc = 0.1
sampling_rate = 10
init_c = 0
