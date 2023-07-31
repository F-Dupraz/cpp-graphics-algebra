#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <math.h>
#include "Matrix.h"

/* Construct from basic */
template <class T>
Matrix<T>::Matrix() {
  m_nRows = 1;
  m_nCols = 1;
  m_nElements = 1;
  m_matrixData = new T[m_nElements];
  m_matrixData[0] = 0.0;
}

/* Construct from number of rows and number of columns */
template <class T>
Matrix<T>::Matrix(int rows, int cols) {
  m_nRows = rows;
  m_nCols = cols;
  m_nElements = rows * cols;
  m_matrixData = new T[m_nElements];
  for(int i = 0; i < m_nElements; i++)
    m_matrixData[i] = 0.0;
}

/* Construct from const linear array */
template <class T>
Matrix<T>::Matrix(int rows, int cols, const T *inputData) {
  m_nRows = rows;
  m_nCols = cols;
  m_nElements = rows * cols;
  m_matrixData = new T[m_nElements];
  for(int i = 0; i < m_nElements; i++)
    m_matrixData[i] = inputData[i];
}

template <class T>
Matrix<T>::~Matrix() {
  //dtor
  if(m_matrixData != nullptr)
    delete[] m_matrixData;
}

/* Resize the matrix from the given values */
template <class T>
bool Matrix<T>::resizeMatrix(int numRows, int numCols) {
  m_nRows = numRows;
  m_nCols = numCols;
  m_nElements = m_nRows * m_nCols;
  delete[] m_matrixData;
  m_matrixData = new T[m_nElements];

  if(m_matrixData != nullptr) {
    for (int i = 0; i < m_nElements; i++) {
      m_matrixData[i] = 0.0;
    }
    return true;
  } else {
    return false;
  }
}

/* Transform the given matrix into the identity matrix */
template <class T>
void Matrix<T>::setToIdentity() {
  if(!IsSquare())
    throw std::invalid_argument("Cannot form an identity matrix, that is not square!");
  for(int row = 0; row < m_nRows; row++) {
    for(int col = 0; col < m_nCols; col++) {
      if(col == row)
        m_matrixData[Sub2Ind(row, col)] = 1.0;
      else
        m_matrixData[Sub2Ind(row, col)] = 0.0;
    }
  }
}

/* Returns the element by its column and its row */
template <class T>
T Matrix<T>::getElement(int row, int col) {
  int linear_index = Sub2Ind(row, col);
  if(linear_index >= 0)
    return m_matrixData[linear_index];
  else
    return 0.0;
}

/* Sets an element in a specified place of the matrix */
template <class T>
bool Matrix<T>::setElement(int row, int col, T elementValue) {
  int linear_index = Sub2Ind(row, col);
  if(linear_index >= 0) {
    m_matrixData[linear_index] = elementValue;
    return true;
  } else {
    return false;
  }
}

/* Returns the total rows */
template <class T>
int Matrix<T>::getRows() {
  return m_nRows;
}

/* Returns the total columns */
template <class T>
int Matrix<T>::getColumns() {
  return m_nCols;
}

template <class T>
int Matrix<T>::Sub2Ind(int row, int col) {
  if((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0))
    return (row * m_nCols) + col;
  else
    return -1;
}

/* Swaps rows i and j */
template <class T>
void Matrix<T>::SwapRow(int i, int j) {
  T *tempRow = new T[m_nCols];
  for(int k = 0; k < m_nCols; k++)
    tempRow[k] = m_matrixData[Sub2Ind(i, k)];
  for(int k = 0; k < m_nCols; k++)
    m_matrixData[Sub2Ind(i, k)] = m_matrixData[Sub2Ind(j, k)];
  for(int k = 0; k < m_nCols; k++)
    m_matrixData[Sub2Ind(j, k)] = tempRow[k];
  delete[] tempRow;
}

/* Adds a multiple row j to row i */
template <class T>
void Matrix<T>::MultAdd(int i, int j, T multFactor) {
  for(int k = 0; k < m_nCols; k++)
    m_matrixData[Sub2Ind(i, k)] += (m_matrixData[Sub2Ind(j, k)] * multFactor);
}

/* Multiplies a row by a given value */
template <class T>
void Matrix<T>::MultRow(int i, T multFactor) {
  for(int k = 0; k < m_nCols; k++)
    m_matrixData[Sub2Ind(i, k)] *= multFactor;
}

/* Finds the row with the column given */
template <class T>
int Matrix<T>::FindRowWithMaxElement(int colNumber, int startingRow) {
  T tempVuale = m_matrixData[Sub2Ind(startingRow, colNumber)];
  int rowIndex = startingRow;
  for(int k = startingRow + 1; k < m_nRows; k++) {
    if(fabs(m_matrixData[Sub2Ind(k, colNumber)]) > fabs(tempVuale)) {
      rowIndex = k;
      tempVuale = m_matrixData[Sub2Ind(k, colNumber)];
    }
  }
  return rowIndex;
}

/* Checks if it's square */
template <class T>
bool Matrix<T>::IsSquare() {
  if(m_nRows == m_nCols)
    return true;
  else
    return false;
}

/* Tests if whether the marix is in row echelon form */
template <class T>
bool Matrix<T>::IsRowEchelon() {
  T cumulative_sum = static_cast<T>(0.0);
  for(int i = 1; i < m_nRows; i++) {
    for(int j = 0; j < i; j++) {
      cumulative_sum += m_matrixData[Sub2Ind(i, j)];
    }
  }
  return CloseEnough(cumulative_sum, 0.0);
}

/* --- Operator Methods --- */

template <class T>
bool Matrix<T>::CloseEnough(T f1, T f2) {
  return fabs(f1-f2) < 1e-9;
}

/* Equals operator */
template <class T>
bool Matrix<T>::operator== (const Matrix<T>& lhs) {
  if((this -> m_nRows != lhs.m_nRows) && (this -> m_nCols != lhs.m_nCols))
    return false;
  bool flag = true;
  for(int i = 0; i < this -> m_nElements; i++) {
    if(!CloseEnough(this->m_matrixData[i], lhs.m_matrixData[i]))
      flag = false;
  }
  return flag;
}

/* Compares the matrix1 to the given matrix whit a certain tolerance */
template <class T>
bool Matrix<T>::Compare(const Matrix<T>& matrix1, double tolerance) {
  int numRows1 = matrix1.m_nRows;
  int numCols1 = matrix1.m_nCols;
  if((numRows1 != m_nRows) || (numCols1 != m_nCols))
    return false;
  double cumulativeSum = 0.0;
  for(int i = 0; i < m_nElements; i++) {
    T element1 = matrix1.m_matrixData[i];
    T element2 = m_matrixData[i];
    cumulativeSum += ((element1 - element2) * (element1 - element2));
  }
  double finalValue = sqrt(cumulativeSum / ((numRows1 * numRows1) - 1));
  if(finalValue < tolerance)
    return true;
  else
    return false;
}

/* Matrix plus matrix */
template <class T>
Matrix<T> operator+ (const Matrix<T>& lhs, const Matrix<T>& rhs) {
  int numRows = lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i < numElements; i++)
    tempResult[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];
  Matrix<T> result(numRows, numCols);
  for(int i = 0; i < numRows; i++) {
    for(int j = 0; j < numCols; j++) {
      result.setElement(i, j, tempResult[(i * numCols) + j]);
    }
  }
  delete[] tempResult;
  return result;
}

/* Scalar plus matrix */
template <class T>
Matrix<T> operator+ (const T& lhs, const Matrix<T>& rhs) {
  int numRows = lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i < numElements; i++)
    tempResult[i] = lhs + rhs.m_matrixData[i];
  Matrix<T> result(numRows, numCols);
  for(int i = 0; i < numRows; i++) {
    for(int j = 0; j < numCols; j++) {
      result.setElement(i, j, tempResult[(i * numCols) + j]);
    }
  }
  delete[] tempResult;
  return result;
}

/* Matrix plus scalar */
template <class T>
Matrix<T> operator+ (const Matrix<T>& lhs, const T& rhs) {
  int numRows = lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i < numElements; i++)
    tempResult[i] = lhs.m_matrixData[i] + rhs;
  Matrix<T> result(numRows, numCols);
  for(int i = 0; i < numRows; i++) {
    for(int j = 0; j < numCols; j++) {
      result.setElement(i, j, tempResult[(i * numCols) + j]);
    }
  }
  delete[] tempResult;
  return result;
}

/* Matrix minus matrix */
template <class T>
Matrix<T> operator- (const Matrix<T>& lhs, const Matrix<T>& rhs) {
  int numRows = lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i < numElements; i++)
    tempResult[i] = lhs.m_matrixData[i] - rhs.m_matrixData[i];
  Matrix<T> result(numRows, numCols);
  for(int i = 0; i < numRows; i++) {
    for(int j = 0; j < numCols; j++) {
      result.setElement(i, j, tempResult[(i * numCols) + j]);
    }
  }
  delete[] tempResult;
  return result;
}

/* Scalar minus matrix */
template <class T>
Matrix<T> operator- (const T& lhs, const Matrix<T>& rhs) {
  int numRows = lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i < numElements; i++)
    tempResult[i] = lhs - rhs.m_matrixData[i];
  Matrix<T> result(numRows, numCols);
  for(int i = 0; i < numRows; i++) {
    for(int j = 0; j < numCols; j++) {
      result.setElement(i, j, tempResult[(i * numCols) + j]);
    }
  }
  delete[] tempResult;
  return result;
}

/* Matrix minus scalar */
template <class T>
Matrix<T> operator- (const Matrix<T>& lhs, const T& rhs) {
  int numRows = lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i < numElements; i++)
    tempResult[i] = lhs.m_matrixData[i] - rhs;
  Matrix<T> result(numRows, numCols);
  for(int i = 0; i < numRows; i++) {
    for(int j = 0; j < numCols; j++) {
      result.setElement(i, j, tempResult[(i * numCols) + j]);
    }
  }
  delete[] tempResult;
  return result;
}

/* Matrix times matrix */
template <class T>
Matrix<T> operator* (const Matrix<T>& lhs, const Matrix<T>& rhs) {
  int r_numRows = rhs.m_nRows;
  int r_numCols = rhs.m_nCols;
  int l_numRows = lhs.m_nRows;
  int l_numCols = lhs.m_nCols;

  if(l_numCols == r_numRows) {
    T *tempResult = new T[lhs.m_nRows * rhs.m_nCols];
    for(int lhsRow = 0; lhsRow < l_numRows; lhsRow++) {
      for(int rhsCol = 0; rhsCol < r_numCols; rhsCol++) {
        T elementResult = 0.0;
        for(int lhsCol = 0; lhsCol < l_numCols; lhsCol++) {
          int lhsLinearIndex = (lhsRow * l_numCols) + lhsCol;
          int rhsLinearIndex = (lhsCol * r_numCols) + rhsCol;
          elementResult += (lhs.m_matrixData[lhsLinearIndex] * rhs.m_matrixData[rhsLinearIndex]);
        }
        int resultLinearindex = (lhsRow * r_numCols) + rhsCol;
        tempResult[resultLinearindex] = elementResult;
      }
    }
    Matrix<T> result(l_numRows, r_numCols);
    for(int i = 0; i < l_numRows; i++) {
      for(int j = 0; j < l_numCols; j++) {
        result.setElement(i, j, tempResult[(i * r_numCols) + j]);
      }
    }
    delete[] tempResult;
    return result;
  } else {
    Matrix<T> result(1, 1);
    return result;
  }
}

/* Scalar times matrix */
template <class T>
Matrix<T> operator* (const T& lhs, const Matrix<T>& rhs) {
  int numRows = lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i < numElements; i++)
    tempResult[i] = lhs * rhs.m_matrixData[i];
  Matrix<T> result(numRows, numCols);
  for(int i = 0; i < numRows; i++) {
    for(int j = 0; j < numCols; j++) {
      result.setElement(i, j, tempResult[(i * numCols) + j]);
    }
  }
  delete[] tempResult;
  return result;
}

/* Matrix times scalar */
template <class T>
Matrix<T> operator* (const Matrix<T>& lhs, const T& rhs) {
  int numRows = lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i < numElements; i++)
    tempResult[i] = lhs.m_matrixData[i] * rhs;
  Matrix<T> result(numRows, numCols);
  for(int i = 0; i < numRows; i++) {
    for(int j = 0; j < numCols; j++) {
      result.setElement(i, j, tempResult[(i * numCols) + j]);
    }
  }
  delete[] tempResult;
  return result;
}

template <class T>
bool Matrix<T>::Separate(Matrix<T> *matrix1, Matrix<T> *matrix2, int colNum) {
  int numRows = m_nRows;
  int numCols1 = colNum;
  // int numCols2 = m_nCols - colNum;

  matrix1->resizeMatrix(numRows, numCols1);
  matrix2->resizeMatrix(numRows, numCols1);

  for(int row = 0; row < m_nRows; row++) {
    for(int col = 0; col < m_nCols; col++) {
      if(col < colNum) {
        matrix1->setElement(row, col, this->getElement(row, col));
      } else {
        matrix2->setElement(row, col - colNum, this->getElement(row, col));
      }
    }
  }
  return true;
}

template <class T>
bool Matrix<T>::Join(const Matrix<T>& matrix2) {
  int numRows1 = m_nRows;
  int numCols1 = m_nCols;
  int numRows2 = matrix2.m_nRows;
  int numCols2 = matrix2.m_nCols;

  if(numRows1 != numRows2)
    throw std::invalid_argument("Attempt to join matrices with different number of rows is invalid!");

  T* newMatrixData = new T[numRows1 * (numCols1 + numCols2)];
  int linearIndex, resultLinearIndex;
  for(int i = 0; i < numRows1; i++) {
    for(int j = 0; j < (numCols1 + numCols2); j++) {
      resultLinearIndex = (i * (numCols1 + numCols2)) + j;
      if(j < numCols1) {
        linearIndex = (i * numCols1) + j;
        newMatrixData[resultLinearIndex] = m_matrixData[linearIndex];
      } else {
        linearIndex = (i * numCols2) + (j - numCols1);
        newMatrixData[resultLinearIndex] = matrix2.m_matrixData[linearIndex];
      }
    }
  }

  m_nCols = numCols1 + numCols2;
  m_nElements = m_nCols * m_nRows;
  delete[] m_matrixData;
  m_matrixData = new T[m_nElements];
  for(int i = 0; i < m_nElements; i++) {
    m_matrixData[i] = newMatrixData[i];
  }

  return true;
}

/* Computes the matrix inverse */

template <class T>
bool Matrix<T>::Inverse() {
  if(!IsSquare())
    throw std::invalid_argument("Cannot Compute the inverse of a matrix that is not square!");

  Matrix<T> identityMatrix(m_nRows, m_nCols);
  identityMatrix.setToIdentity();

  int originalNumCols = m_nCols;
  Join(identityMatrix);

  int cCol, cRow;
  int max_count = 100;
  int my_count = 0;
  bool complete_flag = false;

  while((!complete_flag) && (my_count < max_count)) {
    for(int diagIndex = 0; diagIndex < m_nRows; diagIndex++) {
      cRow = diagIndex;
      cCol = diagIndex;
      int max_index = FindRowWithMaxElement(cCol, cRow);

      if(max_index != cRow)
        SwapRow(cRow, max_index);

      if(m_matrixData[Sub2Ind(cRow, cCol)] != 1.0) {
        T mult_factor = 1.0 / m_matrixData[Sub2Ind(cRow, cCol)];
        MultRow(cRow, mult_factor);
      }

      for(int row_index = cRow + 1; row_index < m_nRows; row_index++) {
        if(!CloseEnough(m_matrixData[Sub2Ind(row_index, cCol)], 0.0)) {
          int row_one_index = cCol;
          T current_element_value = m_matrixData[Sub2Ind(row_index, cCol)];
          T row_one_value = m_matrixData[Sub2Ind(row_one_index, cCol)];

          if(!CloseEnough(row_one_value, 0.0)) {
            T correction_factor = -(current_element_value / row_one_value);
            MultAdd(row_index, row_one_index, correction_factor);
          }
        }
      }
      for(int col_index = cCol + 1; col_index < originalNumCols; col_index++) {
        if(!CloseEnough(m_matrixData[Sub2Ind(cRow, col_index)], 0.0)) {
          int row_one_index = col_index;
          T current_element_value = m_matrixData[Sub2Ind(cRow, col_index)];
          T row_one_value = m_matrixData[Sub2Ind(row_one_index, col_index)];
          if(!CloseEnough(row_one_value, 0.0)) {
            T correction_factor = -(current_element_value / row_one_value);
            MultAdd(cRow, row_one_index, correction_factor);
          }
        }
      }
    }

    Matrix<T> leftHalf;
    Matrix<T> rightHalf;

    this->Separate(&leftHalf, &rightHalf, originalNumCols);

    if(leftHalf == identityMatrix) {
      complete_flag = true;
      m_nCols = originalNumCols;
      m_nElements = m_nCols * m_nRows;
      delete[] m_matrixData;
      m_matrixData = new T[m_nElements];
      for(int i = 0; i < m_nElements; i++)
        m_matrixData[i] = rightHalf.m_matrixData[i];
    }
    my_count++;
  }
  return complete_flag;
}

/* Converts to row echelon form using Gaussian elemination */
template <class T>
Matrix<T> Matrix<T>::RowEchelon() {
  if(m_nCols < m_nRows)
    throw std::invalid_argument("The matrix must have at least the same number of columns as rows");
  int cRow, cCol;
  int max_count = 100;
  int my_count = 0;
  bool complete_flag = false;
  while((!complete_flag) && (my_count < max_count)) {
    for(int diag_index = 0; diag_index < m_nRows; diag_index++) {
      cRow = diag_index;
      cCol = diag_index;
      for(int row_index = cRow + 1; row_index < m_nRows; row_index++) {
        if(!CloseEnough(m_matrixData[Sub2Ind(row_index, cCol)], 0.0)) {
          int row_one_index = cCol;
          T current_element_value = m_matrixData[Sub2Ind(row_index, cCol)];
          T row_one_value = m_matrixData[Sub2Ind(row_one_index, cCol)];
          if(!CloseEnough(row_one_value, 0.0)) {
            T correction_factor = -(current_element_value / row_one_value);
            MultAdd(row_index, row_one_index, correction_factor);
          }
        }
      }
    }
    complete_flag = this->IsRowEchelon();
    my_count++; 
  }
  Matrix<T> outputMatrix(m_nRows, m_nCols, m_matrixData);
  return outputMatrix;
}
