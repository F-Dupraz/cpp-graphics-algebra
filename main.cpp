#include <iostream>
#include <vector>
#include "Matrix.cpp"
#include "LinSolve.h"

/* Prints a matrix */
template <class T>
void printMatrix(Matrix<T>& printed_matrix) {
  int number_columns = printed_matrix.getColumns();
  int number_rows = printed_matrix.getRows();
  for(int i = 0; i < number_rows; i++) {
    for(int j = 0; j < number_columns; j++) {
      std::cout << std::fixed << printed_matrix.getElement(i, j) << " ";
    }
    std::cout << std::endl;
  }
}

int main() {
  std::cout << std::endl;
  std::cout << "This is the matrix we have to resolve: " << std::endl;
  float my_test_matrix_data[9] = {1, -2, 1, 0, 2, -8, -4, 5, 9};
  Matrix<float> MyTestMatrix(3, 3, my_test_matrix_data);
  printMatrix(MyTestMatrix);

  std::cout << std::endl;
  std::cout << "This is the vector with the results of each row: " << std::endl;
  std::vector<float> my_test_matrix_results {0, 8, -9};
  for(int i = 0; i < my_test_matrix_results.size(); i++) {
    std::cout << my_test_matrix_results[i] << std::endl;
  }

  std::vector<float> testResult = LinearSolve<float>(&MyTestMatrix, &my_test_matrix_results);
  std::cout << std::endl;
  std::cout << "The result is: " << std::endl;
  for(int i = 0; i < testResult.size(); i++) {
    std::cout << testResult[i] << std::endl;
  }

  return 0;
}
