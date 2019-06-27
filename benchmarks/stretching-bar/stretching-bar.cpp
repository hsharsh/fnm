// Boundary conditions for square plate with velocity on both sides
// Activate the temporary_bc for 4 seconds. Initiate crack just before 15 seconds.

 void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1){
   double bcx, bcy;
   VectorXi left = VectorXi::LinSpaced(5,1,169);
   left = left.array()-1;
   bcx = 0;
   bcy = 0;
   for(int i=0; i < left.size();++i){
     int node_bc = left(i);
     int idof = node_bc*2;
     un(idof) = bcx;
     un(idof+1) = bcy;
     un1(idof) = bcx;
     un1(idof+1) = bcy;
   }
   VectorXi right = VectorXi::LinSpaced(5,42,210);
   right = right.array()-1;
   bcx = 0.00866025;
   bcy = 0.006;
   for(int i=0; i < right.size();++i){
     int node_bc = right(i);
     int idof = node_bc*2;
     vn(idof) = bcx;
     vn(idof+1) = bcy;
     vn1(idof) = bcx;
     vn1(idof+1) = bcy;
   }
 }

void crack_def(vector<int> &discont, map<int,element> &fn_elements){
  vector <int> cracked;
  for (int i = 21; i <= 144; i+=41){
    cracked.push_back(i-1);
  }
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {0.5, NAN, 0.5, NAN};
  }
}
