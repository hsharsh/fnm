#include "utilities.cpp"

void crack_def(vector<int> &discont, map<int,element> &fn_elements){
  vector <int> cracked;
  for (int i = 21; i <= 144; i+=41){
    cracked.push_back(i-1);
  }
  for (int i = 0; i < cracked.size(); ++i){
    discont[cracked[i]] = 1;
    fn_elements[cracked[i]].edge = {0.99, NAN, 0.01, NAN};
  }
}
