////////////////DEBUG//////////////////////
void print_v(vector<ArrayCoordinateWithHeight> vector) {
  for (uint i = 0; i<vector.size(); i++){
    printf("%d %d %.2f, ", vector[i].row, vector[i].col, vector[i].h);
  }
  printf("\n");
}
void print_v(deque<ArrayCoordinateWithHeight> vector) {
  for (uint i = 0; i<vector.size(); i++){
    printf("%d %d %.2f, ", vector[i].row, vector[i].col, vector[i].h);
  }
  printf("\n");
}
///////////////////////////////////////////