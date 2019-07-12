// Boundary conditions for plate with a hole and streching on top and bottom
// Activate only the boundary_conditions function. No temporary bc is needed
void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
  double bcx, bcy;
  // Syntax: LinSpaced(#elements,start,end). Decrease the start and end by 1 (if taken from a 1-indexed mesh), to reflect 0-indexing
  VectorXi bottom(74);
  bottom <<     7,   8,  15,  16, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212,
   213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 461, 462,
   463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478,
   479, 480, 481, 482, 483, 484, 485, 486, 628, 629, 630, 631, 632, 633, 634, 635,
   636, 637, 638, 639, 640, 641, 642, 643, 644, 645;
  bottom = bottom.array()-1;
  bcx = 0;
  bcy = -1e-4;
  for(int i=0; i < bottom.size();++i){
    long node_bc = bottom(i);
    long idof = node_bc*2;
    vn(idof) = bcx;
    vn1(idof) = bcx;
    vn(idof+1) = bcy;
    vn1(idof+1) = bcy;
  }

  VectorXi top(75);
  top <<      1,   4,  11,  12,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
   110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 331, 332,
   333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348,
   349, 350, 351, 352, 353, 354, 355, 356, 538, 539, 540, 541, 542, 543, 544, 545,
   546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556;
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
}

void crack_def(vector<int> &discont, map<int,element> &fn_elements){
  VectorXi cracked(10);
  cracked << 2934, 2935, 2936, 2937, 2938, 2939, 2940, 2941, 2942, 2943;
  cracked = cracked.array()-1;
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {NAN, 0.5, NAN, 0.5};
  }
}
