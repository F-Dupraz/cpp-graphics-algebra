#ifndef MATRIX_H
#define MATRIX_H

template <class T> class Matrix {
  public:
    /* Defines the constructor */
    Matrix();
    Matrix(int rows, int cols);
    Matrix(int rows, int cols, const T *inputData);

    /* Defines the destructor */
    ~Matrix();

    /* Defines a method to change the rows and cols */
    bool resizeMatrix(int numRows, int numCols);
    /* Transform the given matrix into the identity matrix */
    void setToIdentity();
    /* Gives the reverse matrix */
    bool Inverse();
    /* Converts to row echelon form */
    Matrix<T> RowEchelon();

    /* --- Methods to access to the matrix --- */

    /* Gets one element of the matrix */
    T getElement(int row, int col);
    /* Set one element in the matrix */
    bool setElement(int row, int col, T elementValue);
    /* Gets the rows of the matrix */
    int getRows();
    /* Gets the columns of the matrix */
    int getColumns();

    /* --- Methods for the operators --- */

    /* Overload == operator */
    bool operator== (const Matrix<T>& rhs);
    bool Compare (const Matrix<T>& matrix1, double tolerance);

    /* Overload + operator */
    template <class U> friend Matrix<U> operator+ (const Matrix<U>& lhs, const Matrix<U>& rhs);
    template <class U> friend Matrix<U> operator+ (const U& lhs, const Matrix<U>& rhs);
    template <class U> friend Matrix<U> operator+ (const Matrix<U>& lhs, const U& rhs);

    /* Overload - operator */
    template <class U> friend Matrix<U> operator- (const Matrix<U>& lhs, const Matrix<U>& rhs);
    template <class U> friend Matrix<U> operator- (const U& lhs, const Matrix<U>& rhs);
    template <class U> friend Matrix<U> operator- (const Matrix<U>& lhs, const U& rhs);

    /* Overload * operator */
    template <class U> friend Matrix<U> operator* (const Matrix<U>& lhs, const Matrix<U>& rhs);
    template <class U> friend Matrix<U> operator* (const U& lhs, const Matrix<U>& rhs);
    template <class U> friend Matrix<U> operator* (const Matrix<U>& lhs, const U& rhs);

    /* Separates the given matrix in two matrices */
    bool Separate(Matrix<T> *matrix1, Matrix<T> *matrix2, int colNum);

  public:
    int Sub2Ind(int row, int col);
    bool IsSquare();
    bool IsRowEchelon();
    static int Rank(const Matrix<T> &inputMatrix);
    bool CloseEnough(T f1, T f2);
    void SwapRow(int i, int j);
    void MultAdd(int i, int j,  T multFactor);
    void MultRow(int i, T multFactor);
    bool Join(const Matrix<T>& matrix2);
    int FindRowWithMaxElement(int colNumber, int startingRow);
    void printMatrix();

  private:
    T *m_matrixData;
    int m_nRows, m_nCols, m_nElements;
};

#endif // MATRIX_H
