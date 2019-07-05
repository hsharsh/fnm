// Boundary conditions for square plate with velocity on both sides
// Activate the temporary_bc for 4 seconds. Initiate crack just before 15 seconds.

 void boundary_conditions(VectorXd &un, VectorXd &un1, VectorXd &vn, VectorXd &vn1, VectorXd &fg){
   double bcx, bcy;
   VectorXi left(4);
   left << 203,  204,  205, 1409;
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
   VectorXi right(4);
   right <<  90,   91,   92, 1189;
   right = right.array()-1;
   bcx = 0;
   bcy = 0;
   for(int i=0; i < right.size();++i){
     int node_bc = right(i);
     int idof = node_bc*2;
     un(idof) = bcx;
     un(idof+1) = bcy;
     un1(idof) = bcx;
     un1(idof+1) = bcy;
   }

   VectorXi top(6);
   top << 495,  496,  497,  498, 2868, 2869;
   top = top.array()-1;
   bcy = -1e-4;
   for(int i=0; i < top.size();++i){
     int node_bc = top(i);
     int idof = node_bc*2;
     vn(idof+1) = bcy;
     vn1(idof+1) = bcy;
   }
 }

 void crack_def(vector<int> &discont, map<int,element> &fn_elements){
   VectorXi cracked(4);
   cracked << 716, 741, 766,  791;
   cracked = cracked.array()-1;
   for (int i = 0; i < cracked.size(); ++i){
     discont[cracked[i]] = 1;
     fn_elements[cracked[i]].edge = {0.5, NAN, 0.5, NAN};
   }
 }
