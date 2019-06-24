// Boundary conditions for square plate with velocity on both sides
// Activate the temporary_bc for 4 seconds. Initiate crack just before 15 seconds.

 void temporary_bc(VectorXd &vn, VectorXd &vn1){
   double bc;
   VectorXi left = VectorXi::LinSpaced(5,1,169);
   left = left.array()-1;
   bc = 0.01;
   for(int i=0; i < left.size();i++){
     int node_bc = left(i);
     int idof = node_bc*2;
     vn(idof) = bc;
     vn1(idof) = bc;
   }
 }

void crack_def(vector<int> &discont, map<int,element> &fn_elements){
  vector <int> cracked;
  for (int i = 21; i <= 144; i+=21){
    cracked.push_back(i-1);
  }
  for (int i = 0; i < cracked.size(); i++){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {0.5, NAN, 0.5, NAN};
  }
}
