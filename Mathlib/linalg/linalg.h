#include <bits/stdc++.h>
#include <errno.h>
#include <stddef.h>
#include <stdint.h>


// it does not have to be the most efficient implementation
// the compiler will take care of that 
// the performance aspects will be handled later 
// this class handles 2 dimensional matrices
// later we will introduce a system that can rival numpy

#ifndef LINALG
#define LINALG
#endif

#define MAT_AT(m, i, j, stride) (m)[(i) * (stride) + (j)]

template <class T>

class Matrix{
    private:
        T* data;
        size_t nRows;
        size_t nCols;
        size_t nElements;
        size_t allocatedMemory;

        double determinant;
        bool isInvertible;
        bool isSymmetric;
        bool isAntisymmetric;
        bool isInvertible;

    public:
        // constructors
        Matrix(); // default constructor
        Matrix(size_t rows, size_t cols); // allocate and set all elements to zero
        Matrix(size_t rows, size_t cols, const T* input); // standard constructor
        Matrix(Matrix<T>& m); // copy constructor
        Matrix(Matrix<T>* m); // another copy constructor

        ~Matrix(); // destructor

        bool operator== (const Matrix<T>& rhs); // declare the comparision operator and provide a value to test against, rhs

        template <class U> friend Matrix<U> operator+ (const Matrix<U>& lhs, const Matrix<U>& rhs); // add two matrices
        template <class U> friend Matrix<U> operator* (Matrix<U>& lhs, Matrix<U> &rhs); // add two matrices
        template <class U> friend Matrix<U> operator+ (const U& lhs, const Matrix<U> &rhs); // add a matrix to a scalar
        template <class U> friend Matrix<U> operator+ (const Matrix<U> &lhs, const U& rhs); // add a scalar to a matrix

        template <class U> friend Matrix<U> operator+ (const Matrix<U>& lhs, const Matrix<U>& rhs); // subtract two matrices
        template <class U> friend Matrix<U> operator+ (const U& lhs, const Matrix<U> &rhs); // subtract a matrix from a scalar
        template <class U> friend Matrix<U> operator+ (const Matrix<U> &lhs, const U& rhs); // subtract a scalar from a matrix

        template <class U> friend Matrix<U> operator+ (const Matrix<U>& lhs, const Matrix<U>& rhs); // multiply two matrices
        template <class U> friend Matrix<U> operator+ (const U& lhs, const Matrix<U> &rhs); // linearly scale a matrix
        template <class U> friend Matrix<U> operator+ (const Matrix<U> &lhs, const U& rhs); // ??
        
        template <class U> friend U& operator[] (const Matrix<U> &target, int row, int col); // array indexing operation 


        void setElement(int row, int col, T& value);
        void setElement(int index, T& value);
        T& getElement(int row, int col);
        T& getElement(int index);
        size_t getNumRows(){return this->nRows;}
        size_t getNumCols(){return this->nCols;}
        size_t getNumElements(){return this->nElements;}
        void swapRows(size_t r1, size_t r2);
        void swapCols(size_t c1, size_t c2);
        Matrix<T> invert();
        size_t findPivot(size_t r1, size_t r2, size_t loc);

        bool resize(size_t row, size_t col);
        bool resize(size_t row, size_t col, T& value);
        void setToIdentity();
        void setToZero();
        template <class U> friend Matrix<U> multiply(const Matrix<U> &a, const Matrix<U> &b);

        // perhaps we could also define an operator to provide indexing operations

};


/*##########################################################  IMPLEMENTATION  ##########################################################*/

// you might be wondering why the implementation is provided within the header file itself
// this is because anything that uses templating must be processed entirely in the include statements and thus must be in the header file 






template <class T>
Matrix<T>::Matrix(){
    // default constructor implementation 
    this->nRows = 1;
    this->nCols = 1;
    this->nElements = 1;
    this->data = new T(sizeof(T));
    data[0] = 0; // the code is meant to deal with primarily number types, so we assign a value of zero to 0d scalar
}

template <class T>
Matrix<T>::Matrix(size_t rows, size_t cols){
    this->nRows = rows;
    this->nCols = cols;
    this->nElements = rows*cols;
    this->data = new T(sizeof(T)*(this->nElements));
    for(int i = 0; i < this->nRows; i++){
        for(int j = 0; j < this->nCols; j++){
            MAT_AT(this->data, i, j, this->nCols) = 0;
        }
    }
}

template <class T>
Matrix<T>::Matrix(size_t rows, size_t cols, const T* inputvector){
    this->nRows = rows;
    this->nCols = cols;
    this->nElements = rows*cols;
    this->data = new T(sizeof(T)*(this->nElements));
    for(int i = 0; i < this->nElements; i++) (this->data)[i] = inputvector[i];
    return;
}

template <class T>
Matrix<T>::Matrix(Matrix<T>&m){
    // copy constructor
    this->nRows = m.getNumRows();
    this->nCols = m.getNumCols();
    this->nElements = m.getNumElements();
    this->data = new T(sizeof(T)*(this->nElements));
    for(int i = 0; i < this->nElements; i++) (this->data)[i] = m.getElement(i);
    return;
}

template <class T>
Matrix<T>::Matrix(Matrix<T>* m){
    // copy constructor using pointer/this
    this->nRows = m->getNumRows();
    this->nCols = m->getNumCols();
    this->nElements = m->getNumElements();
    this->data = new T(sizeof(T)*(this->nElements));
    for(int i = 0; i < this->nElements; i++) (this->data)[i] = m->getElement(i);
    return;
}

// Destructor //
template <class T>
Matrix<T>::~Matrix(){
    if(this->data) delete this->data; // relinquish allocated memory back to the OS
}

// configuration functions //

template <class T>
T& Matrix<T>::getElement(int row, int col){
    int correctedRow = row % this->nRows;
    int correctedCol = col % this->nCols;
    int linearIndex = correctedRow*(this->nCols) + correctedCol;
    return (this->data)[linearIndex];
}

template <class T> 
T& Matrix<T>::getElement(int index){
    int correctedIndex = index%(this->numElements);
    return (this->data)[correctedIndex];
}

template <class T>
void Matrix<T>::setElement(int row, int col, T& value){
    int correctedRow = row % this->nRows;
    int correctedCol = col % this->nCols;
    int linearIndex = correctedRow*(this->nCols) + correctedCol;
    (this->data)[linearIndex] = value;
    return;
}


template <class T>
bool Matrix<T>::resize(size_t rows, size_t cols){
    // if available memory is greater than the requested memory, we will simply not relinquish it 
    // please note resizing in this manner will destroy all elements that fall out of memory 
    if(rows*cols < this->allocatedMemory){
        for(int i = 0; i < (rows*cols); i++){
            (this->data)[i] = 0;
        }
    }
    // if available memory is lesser than the requested memory, we will first allocate a new chunk of memory, 
    // and copy our existing items into, and fill the remaining values with zeros
    else if(rows >= this->nRows && cols >= this->nCols){
        T* newlyAllocatedMemory = new T(sizeof(T)*rows*cols);
        int totalElements = rows*cols;
        for(int i = 0; i < totalElements; i++) newlyAllocatedMemory[i] = 0;
        for(int i = 0; i < this->nRows; i++){
            for(int j = 0; j < this->nCols; j++){
                newlyAllocatedMemory[i*cols + j] = (this->data)[i*nCols + j];
            }
        }
        delete this->data; // free previously allocated memory
        this->data = newlyAllocatedMemory;
        this->allocatedMemory = rows*cols;
    }
    else{
        T* newlyAllocatedMemory = new T(sizeof(T)*rows*cols);
        // destroy all data
        int totalElements = rows*cols;
        for(int i = 0; i < totalElements; i++){
            newlyAllocatedMemory[i] = 0;
        }
        delete this->data; // free previously allocated memory
        this->data = newlyAllocatedMemory;
        this->allocatedMemory = rows*cols;
    }
    this->nRows = rows;
    this->nCols = cols;
    this->nElements = rows*cols;
    return true;
}

template <class T>
bool Matrix<T>::resize(size_t rows, size_t cols, T& value){
    if(rows*cols < this->allocatedMemory){
        for(int i = 0; i < (rows*cols); i++){
            (this->data)[i] = value;
        }
    }
    else{
        T* newlyAllocatedMemory = new T(sizeof(T)*rows*cols);
        // overwrite all data
        int totalElements = rows*cols;
        for(int i = 0; i < totalElements; i++){
            newlyAllocatedMemory[i] = value;
        }
        delete this->data; // free previously allocated memory 
        this->data = newlyAllocatedMemory;
        this->allocatedMemory = rows*cols;
    }
    this->nRows = rows;
    this->nCols = cols;
    this->nElements = rows*cols;
    return true;
}

template <class T>
void Matrix<T>::setToZero(){
    int count = this->nElements;
    for(int i = 0; i < count; i++) (this->data)[i] = 0;
}

template <class T> 
void Matrix<T>::setToIdentity(){
    assert(this->nRows == this->nCols);
    this->setToZero();
    for(int i = 0; i < nRows; i++) (this->data)[i*nCols + i] = 1;
    return;
}

// OPERATORS //


template <class T>
bool Matrix<T>::operator== (const Matrix <T> &rhs){
    if(this->nRows != rhs.nRows) return false;
    if(this->nCols != rhs.nCols) return false;
    for(int i = 0; i < a.nElements; i++) if((this->data)[i] != (rhs->data)[i]) return false;
    return true;
}

template <class T>
Matrix<T> operator+ (Matrix <T> &lhs, Matrix <T> &rhs){
    // adds one matrix to another and returns the result in a new matrix
    // since this is a friend function, it can access the private members of the class;
    assert((a.nRows == b.nRows) && (a.nCols == b.nCols));
    Matrix result(lhs);
    for(int i = 0; i < a.nElements; i++) (result.data)[i] += (rhs.data)[i];
    return result;
}

// ADDITION OPERATOR + //


template <class T>
Matrix<T> operator+ (T& lhs, Matrix <T> &rhs){
    // adds a constant value to every element of a matrix(rhs) and returns the result in a new matrix;
    Matrix result(rhs);
    for(int i = 0; i < rhs.nElements; i++) (result.data)[i] += lhs;
    return result;
}

template <class T>
Matrix<T> operator+ (Matrix <T> &lhs, T& rhs){
    // adds 
    Matrix result(lhs);
    for(int i = 0; i < lhs.nElements; i++) (result.data)[i] += rhs;
    return result;
}

// SUBTRACTION OPERATOR - //

template <class T>
Matrix<T> operator- (Matrix <T> &lhs, Matrix <T> &rhs){
    assert((a.nRows == b.nRows) && (a.nCols == b.nCols));
    Matrix result(lhs);
    for(int i = 0; i < a.nElements; i++) (result.data)[i] -= (rhs.data)[i];
    return result;
}

template <class T>
Matrix<T> operator- (T& lhs, Matrix <T> &rhs){
    Matrix result(rhs);
    for(int i = 0; i < rhs.nElements; i++) (result.data)[i] = lhs - (result.data)[i];
    return result;
}

template <class T>
Matrix<T> operator- (Matrix <T> &lhs, T& rhs){
    Matrix result(lhs);
    for(int i = 0; i < lhs.nElements; i++) (result.data)[i] -= rhs;
    return result;
}

// MULTIPLICATION OPERATOR * //

template <class T>
Matrix<T> operator* (Matrix <T> &lhs, Matrix <T> &rhs){
    // multiplication between two matrices
    assert((lhs.nCols == rhs.nRows));
    Matrix result(lhs.nRows, rhs.nCols);
    size_t eye = lhs.nRows;
    size_t jay = rhs.nCols;
    size_t kay = lhs.nCols;
    T& accumulator = 0;
    for(int i = 0; i < eye; i++){
        for(int j = 0; j < jay; j++){
            accumulator = 0;
            for(int k = 0; k < kay; k++){
                accumulator += lhs.getElement(i, k)*rhs.getElement(k, j);
            }
            result.setElement(i, j, accumulator);
        }
    }
    return result;
}


template <class T> 
Matrix <T> operator* (T& lhs, Matrix<T>& rhs){
    Matrix result(rhs);
    for(int i = 0; i < rhs.nElements; i++){
        (result.data)[i] *= lhs;
    }
    return result;
}

template <class T>
Matrix<T> operator* (Matrix<T> &lhs, T& rhs){
    Matrix result(lhs);
    for(int i = 0; i < lhs.nElements; i++) (result.data) *= rhs;
    return result;
}


template <class T>
Matrix<T> mulitply(const Matrix<T> &a, const Matrix<T> &b){
    // performs element wise multiplication 
    assert((a.nRows == b.nRows) && (a.nCols == b.nCols));
    int elCount = a.nElements;
    Matrix result(a.nRows, a.nCols);
    for(int i = 0; i < elCount; i++){
        (result.data)[i] = (a.data)[i]*((b.data)[i]);
    }
} 


// GAUSSIAN ELIMINATION ///////////////////


template <class T>
void Matrix<T>::swapRows(size_t r1, size_t r2){
    // swaps rows, used in Gaussian elimination 
    assert(r1 < this->nRows && r2 < this->nRows);
    T temp;
    int cols = this->nCols;
    for(int i = 0; i < cols; i++){
        temp = (this->data)[r1*(cols) + i];
        (this->data)[r1*(cols) + i] = (this->data)[r2*(cols) + i];
        (this->data)[r2*(cols) + i] = temp;
    }
    return;

}

template <class T>
void Matrix<T>::swapCols(size_t c1, size_t c2){
    // swaps cols, used in Gaussian elimination 
    assert(c1 < this->nCols && c2 < this->nCols);
    T temp;
    int rows = this->nRows;
    int cols = this->nCols;
    for(int i = 0; i < rows; i++){
        temp = (this->data)[i*(this->nCols) + c1];
        (a->data)[i*(cols) + c1] = (a->data)[i*(cols) + c2];
        (a->data)[i*(cols) + c2] = temp;
    }
    return;
}

template <class T> 
size_t Matrix<T>::findPivot(size_t rs, size_t re, size_t loc){
    for(int i = rs; i < re; i++){
        if((this->data)[i*nCols + loc] != static_cast<T>(0.0)) return i;
    }
    return -1;
}


template <class T>
Matrix<T> Matrix<T>::invert(){
    // this is a simple function that performs a matrix inversion
    // if at any point, we have a breakdown of elimination, we terminate and raise an error
    Matrix<T> inverse(this->nRows, this->nCols);
    Matrix<T> temp(this->nRows, this->nCols);
    for(int i = 0; i < n; i++){
        temp.setElement();
    }
    inverse.setToIdentity();
    // the elimination steps;
    /*
    
        1. Whenever the pivot evaluates to zero, swap it out with one of the rows beneath that has a non zero pivot
        2. Reduce all rows underneath by subtracting multiples of the pivot. A simple way would be to simply perform:
            Rp *= p*(mij/p);
    
    */
    int rows = this->nRows;
    int cols = this->nCols;

    for(int i = 0; i < cols; i++){
        // external loop evaulates all pivots
        // evaluate the status
        T p = (this->data)[i*cols + i]; 
        if(p == static_cast<T>(0.0)){
            int index = findPivot(this, i + 1, this->nRows);
            if(index == -1){
                inverse.setToIdentity();
                return inverse;
            }
        }
        for(int j = i + 1; j < cols; j++){
            for(int k = )
        }
    }

}