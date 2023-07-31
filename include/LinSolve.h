#ifndef LINEARSOLVE_H
#define LINEARSOLVE_H

#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <math.h>
#include <vector>
#include "Matrix.h"

/*
 * NO HAY QUE ARREGLAR NADA!!! FUNCIONA AHORAAA!!!
*/

template <typename T>
std::vector<T> LinearSolve(const Matrix<T> *aMatrix, const std::vector<T> *bVector) {
  Matrix<T> inputMatrix = *aMatrix;
  int num_dims = bVector->size();
  T vector_value[num_dims];
  for (int i = 0; i < num_dims; i++)
    vector_value[i] = (*bVector)[i];
  Matrix<T> bMatrix(num_dims, 1, vector_value);
  inputMatrix.Join(bMatrix);
  Matrix<T> rowEchelonMatrix = inputMatrix.RowEchelon();

  // Create a vector to store the intermediate results
  std::vector<T> intermediate_results(num_dims);

  int num_rows = rowEchelonMatrix.getRows();
  int num_cols = rowEchelonMatrix.getColumns();
  int start_row = num_rows - 1;

  for (int i = start_row; i >= 0; i--) {
    T current_result = rowEchelonMatrix.getElement(i, num_cols - 1);
    T cumulative_sum = static_cast<T>(0.0);

    for (int j = i + 1; j < num_rows; j++)
      cumulative_sum += (rowEchelonMatrix.getElement(i, j) * intermediate_results[j]);

    T final_answer = (current_result - cumulative_sum) / rowEchelonMatrix.getElement(i, i);
    intermediate_results[i] = final_answer;
  }

  // Copy the intermediate results to the output vector
  std::vector<T> output(intermediate_results);

  return output;
}

#endif // LINEARSOLVE_H
