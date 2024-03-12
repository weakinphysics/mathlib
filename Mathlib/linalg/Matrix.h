#ifndef ULTIUltimatrix_H
#define ULTIUltimatrix_H
#endif 

#include <bits/stdc++.h>
#include <errno.h>
#include <stddef.h>
#include <stdint.h>
#include <iostream>


// it does not have to be the most efficient implementation
// the compiler will take care of that 
// the performance aspects will be handled later 
// this class handles 2 dimensional matrices
// later we will introduce a system that can rival numpy


#define MAT_AT(m, i, j, stride) (m)[(i) * (stride) + (j)]

template <class T>

class Ultimatrix{
    private:
        T* data;
        size_t nRows;
        size_t nCols;
        size_t nElements;
        size_t allocatedMemory;

        double determinant;
        int isInvertible;
        bool isSymmetric;
        bool isAntisymmetric;

    public:
        // constructors
        Ultimatrix(); // default constructor
        Ultimatrix(size_t rows, size_t cols); // allocate and set all elements to zero
        Ultimatrix(size_t rows, size_t cols, const T* input); // standard constructor
        Ultimatrix(Ultimatrix<T>& m); // copy constructor
        Ultimatrix(Ultimatrix<T>* m); // another copy constructor

        ~Ultimatrix(); // destructor

        bool operator== (const Ultimatrix<T>& rhs); // declare the comparision operator and provide a value to test against, rhs
        T& operator[] (size_t index); // array indexing operation 

        template <class U> friend Ultimatrix<U> operator+ ( Ultimatrix<U>& lhs,  Ultimatrix<U>& rhs); // add two matrices
        template <class U> friend Ultimatrix<U> operator+ ( U lhs,  Ultimatrix<U> &rhs); // add a Ultimatrix to a scalar
        template <class U> friend Ultimatrix<U> operator+ ( Ultimatrix<U> &lhs,  U rhs); // add a scalar to a Ultimatrix

        template <class U> friend Ultimatrix<U> operator- ( Ultimatrix<U>& lhs,  Ultimatrix<U>& rhs); // subtract two matrices
        template <class U> friend Ultimatrix<U> operator- ( U lhs,  Ultimatrix<U> &rhs); // subtract a Ultimatrix from a scalar
        template <class U> friend Ultimatrix<U> operator- ( Ultimatrix<U> &lhs,  U rhs); // subtract a scalar from a Ultimatrix

        template <class U> friend Ultimatrix<U> operator* ( Ultimatrix<U>& lhs,  Ultimatrix<U>& rhs); // multiply two matrices
        template <class U> friend Ultimatrix<U> operator* ( U lhs,  Ultimatrix<U> &rhs); // linearly scale a Ultimatrix
        template <class U> friend Ultimatrix<U> operator* ( Ultimatrix<U> &lhs,  U rhs); // ??
        
        
        // template <class U> friend std::ostream& operator<<(ostream& os,  Ultimatrix<U> &output); // print Ultimatrix

        template <class U> friend Ultimatrix<U> multiply(Ultimatrix<U> &a, Ultimatrix<U> &b); // element wise multiplication

        void setElement(int row, int col, T value);
        void setElement(int index, T value);
        T& getElement(int row, int col);
        T& getElement(int index);
        size_t getNumRows(){return this->nRows;}
        size_t getNumCols(){return this->nCols;}
        size_t getNumElements(){return this->nElements;}
        void swapRows(size_t r1, size_t r2);
        void swapCols(size_t c1, size_t c2);
        int getInvertible();
        Ultimatrix<T> invert();
        size_t findPivot(size_t r1, size_t r2, size_t loc);
        void printMatrix();
        void scaleRow(size_t r, T& mult);
        bool resize(size_t row, size_t col);
        bool resize(size_t row, size_t col, T& value);
        bool resize(size_t row, size_t col, T value);
        void setToIdentity();
        void setToZero();
        

        // perhaps we could also define an operator to provide indexing operations

};


/*##########################################################  IMPLEMENTATION  ##########################################################*/

// you might be wondering why the implementation is provided within the header file itself
// this is because anything that uses templating must be processed entirely in the include statements and thus must be in the header file 


template <class T>
Ultimatrix<T>::Ultimatrix(){
    // default constructor implementation 
    this->nRows = 1;
    this->nCols = 1;
    this->nElements = 1;
    this->data = new T;
    this->isInvertible = -1;
    data[0] = 0; // the code is meant to deal with primarily number types, so we assign a value of zero to 0d scalar
}

template <class T>
Ultimatrix<T>::Ultimatrix(size_t rows, size_t cols){
    this->nRows = rows;
    this->nCols = cols;
    this->nElements = rows*cols;
    // std::cout<<sizeof(T)<<std::endl;
    this->data = new T[rows*cols]; // we dont need to use sizeof(T)*rows*cols as new figures it out automatically, unlike malloc
    this->isInvertible = -1;
    for(int i = 0; i < this->nRows; i++){
        for(int j = 0; j < this->nCols; j++){
            // MAT_AT(this->data, i, j, this->nCols) = 0;
            (this->data)[i*cols + j] = 0;
        }
    }
}

template <class T>
Ultimatrix<T>::Ultimatrix(size_t rows, size_t cols, const T* inputvector){
    // please ensure no data leaks are present when using an inputvector
    // for(int i = 0; i < rows*cols; i++) std::cout<<inputvector[i]<<" ";
    // std::cout<<std::endl;
    this->nRows = rows;
    this->nCols = cols;
    this->nElements = rows*cols;
    this->data = new T[rows*cols];
    this->isInvertible = -1;
    for(int i = 0; i < (rows*cols); i++) (this->data)[i] = inputvector[i];
}

template <class T>
Ultimatrix<T>::Ultimatrix(Ultimatrix<T>&m){
    // copy constructor
    this->nRows = m.getNumRows();
    this->nCols = m.getNumCols();
    this->nElements = m.getNumElements();
    this->data = new T[m.getNumElements()];
    this->isInvertible = -1;
    for(int i = 0; i < this->nElements; i++) (this->data)[i] = m.getElement(i);
    return;
}

template <class T>
Ultimatrix<T>::Ultimatrix(Ultimatrix<T>* m){
    // copy constructor using pointer/this
    this->nRows = m->getNumRows();
    this->nCols = m->getNumCols();
    this->nElements = m->getNumElements();
    this->data = new T[this->nElements];
    this->isInvertible = -1;
    for(int i = 0; i < this->nElements; i++) (this->data)[i] = m->getElement(i);
    return;
}

// Destructor //
template <class T>
Ultimatrix<T>::~Ultimatrix(){
    if(this->data) delete this->data; // relinquish allocated memory back to the OS
}

// configuration functions //

template <class T>
T& Ultimatrix<T>::getElement(int row, int col){
    int correctedRow = row % this->nRows;
    int correctedCol = col % this->nCols;
    int linearIndex = correctedRow*(this->nCols) + correctedCol;
    return (this->data)[linearIndex];
}

template <class T> 
T& Ultimatrix<T>::getElement(int index){
    int correctedIndex = index%(this->nElements);
    return (this->data)[correctedIndex];
}

template <class T>
void Ultimatrix<T>::setElement(int row, int col, T value){
    int correctedRow = row % this->nRows;
    int correctedCol = col % this->nCols;
    int linearIndex = correctedRow*(this->nCols) + correctedCol;
    (this->data)[linearIndex] = value;
    return;
}

template <class T>
int Ultimatrix<T>::getInvertible(){
    if(this->isInvertible != -1) return this->isInvertible;
    Ultimatrix<T> garbage = this->invert();
    return this->isInvertible;     
}


template <class T>
void Ultimatrix<T>::printMatrix(){
    std::cout<<std::endl;
    int rows = this->nRows;
    int cols = this->nCols;
    std::cout<<"Ultimatrix Object of dimension "<<rows<<" by "<<cols<<std::endl;
    std::cout<<"Matrix Data: "<<std::endl;
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            std::cout<<(this->data)[i*cols + j]<<" ";
        }
        std::cout<<std::endl;
    }
    return;
}


template <class T>
bool Ultimatrix<T>::resize(size_t rows, size_t cols){
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
        T* newlyAllocatedMemory = new T[rows*cols];
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
        T* newlyAllocatedMemory = new T[rows*cols];
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
bool Ultimatrix<T>::resize(size_t rows, size_t cols, T& value){
    if(rows*cols < this->allocatedMemory){
        for(int i = 0; i < (rows*cols); i++){
            (this->data)[i] = value;
        }
    }
    else{
        T* newlyAllocatedMemory = new T[rows*cols];
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
bool Ultimatrix<T>::resize(size_t rows, size_t cols, T value){
    if(rows*cols < this->allocatedMemory){
        for(int i = 0; i < (rows*cols); i++){
            (this->data)[i] = value;
        }
    }
    else{
        T* newlyAllocatedMemory = new T[rows*cols];
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
void Ultimatrix<T>::setToZero(){
    int count = this->nElements;
    for(int i = 0; i < count; i++) (this->data)[i] = 0;
}

template <class T> 
void Ultimatrix<T>::setToIdentity(){
    assert(this->nRows == this->nCols);
    this->setToZero();
    for(int i = 0; i < nRows; i++) (this->data)[i*nCols + i] = 1;
    return;
}



// OPERATORS //


template <class T>
bool Ultimatrix<T>::operator== (const Ultimatrix <T> &rhs){
    if(this->nRows != rhs.nRows) return false;
    if(this->nCols != rhs.nCols) return false;
    for(int i = 0; i < rhs.nElements; i++) if((this->data)[i] != (rhs.data)[i]) return false;
    return true;
}

// template <class T>
// T& Ultimatrix<T>::operator[] (size_t index){
//     return
// }


// template <class T>
// std::ostream& operator<< (const std::ostream& os, const Ultimatrix<T>& m){
//     m.printMatrix();
//     return os;
// }

// ADDITION OPERATOR + //

template <class T>
Ultimatrix<T> operator+(Ultimatrix <T> &lhs, Ultimatrix <T> &rhs){
    // adds one Ultimatrix to another and returns the result in a new Ultimatrix
    // since this is a friend function, it can access the private members of the class;
    assert((lhs.nRows == rhs.nRows) && (lhs.nCols == rhs.nCols));
    Ultimatrix result(lhs);
    for(int i = 0; i < lhs.nElements; i++) (result.data)[i] += (rhs.data)[i];
    return result;
}



template <class T>
Ultimatrix<T> operator+( T lhs,  Ultimatrix <T> &rhs){
    // adds a ant value to every element of a Ultimatrix(rhs) and returns the result in a new Ultimatrix;
    Ultimatrix result(rhs);
    for(int i = 0; i < rhs.nElements; i++) (result.data)[i] += lhs;
    return result;
}

template <class T>
Ultimatrix<T> operator+( Ultimatrix <T> &lhs,  T rhs){
    // adds 
    Ultimatrix result(lhs);
    for(int i = 0; i < lhs.nElements; i++) (result.data)[i] += rhs;
    return result;
}

// SUBTRACTION OPERATOR - //

template <class T>
Ultimatrix<T> operator- ( Ultimatrix <T> &lhs,  Ultimatrix <T> &rhs){
    assert((lhs.nRows == rhs.nRows) && (lhs.nCols == rhs.nCols));
    Ultimatrix result(lhs);
    for(int i = 0; i < lhs.nElements; i++) (result.data)[i] -= (rhs.data)[i];
    return result;
}

template <class T>
Ultimatrix<T> operator- ( T lhs,  Ultimatrix <T> &rhs){
    Ultimatrix result(rhs);
    for(int i = 0; i < rhs.nElements; i++) (result.data)[i] = lhs - (result.data)[i];
    return result;
}

template <class T>
Ultimatrix<T> operator- ( Ultimatrix <T> &lhs,  T rhs){
    Ultimatrix result(lhs);
    for(int i = 0; i < lhs.nElements; i++) (result.data)[i] -= rhs;
    return result;
}

// MULTIPLICATION OPERATOR * //

template <class T>
Ultimatrix<T> operator* ( Ultimatrix <T> &lhs,  Ultimatrix <T> &rhs){
    // multiplication between two matrices
    assert((lhs.nCols == rhs.nRows));
    Ultimatrix<T> result(lhs.nRows, rhs.nCols);
    size_t eye = lhs.nRows;
    size_t jay = rhs.nCols;
    size_t kay = lhs.nCols;
    T accumulator = static_cast<T>(0.0);
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
Ultimatrix <T> operator* ( T lhs,  Ultimatrix<T>& rhs){
    Ultimatrix result(rhs);
    for(int i = 0; i < rhs.nElements; i++){
        (result.data)[i] *= lhs;
    }
    return result;
}

template <class T>
Ultimatrix<T> operator* ( Ultimatrix<T> &lhs,  T rhs){
    Ultimatrix result(lhs);
    for(int i = 0; i < lhs.nElements; i++) (result.data) *= rhs;
    return result;
}


template <class T>
Ultimatrix<T> multiply( Ultimatrix<T> &a,  Ultimatrix<T> &b){
    // performs element wise multiplication 
    assert((a.nRows == b.nRows) && (a.nCols == b.nCols));
    int elCount = a.nElements;
    Ultimatrix <T> result(a.nRows, a.nCols);
    for(int i = 0; i < elCount; i++){
        (result.data)[i] = (a.data)[i]*((b.data)[i]);
    }
    return result;
} 
 

// GAUSSIAN ELIMINATION ///////////////////


template <class T>
void Ultimatrix<T>::swapRows(size_t r1, size_t r2){
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
void Ultimatrix<T>::swapCols(size_t c1, size_t c2){
    // swaps cols, used in Gaussian elimination 
    assert(c1 < this->nCols && c2 < this->nCols);
    T temp;
    int rows = this->nRows;
    int cols = this->nCols;
    for(int i = 0; i < rows; i++){
        temp = (this->data)[i*(this->nCols) + c1];
        (this->data)[i*(cols) + c1] = (this->data)[i*(cols) + c2];
        (this->data)[i*(cols) + c2] = temp;
    }
    return;
}

template <class T> 
size_t Ultimatrix<T>::findPivot(size_t rs, size_t re, size_t loc){
    for(int i = rs; i < re; i++){
        if((this->data)[i*nCols + loc] != static_cast<T>(0.0)) return i;
    }
    return -1;
}


template <class T>
Ultimatrix<T> Ultimatrix<T>::invert(){

    // this is a simple function that performs a Ultimatrix inversion
    // if at any point, we have a breakdown of elimination, we terminate and raise an error
    Ultimatrix<T> inverse(this->nRows, this->nCols);
    if(this->nRows != this->nCols){
        this->isInvertible = false;
        return inverse;
    }

    Ultimatrix<T> temp(this);
    Ultimatrix<T> L(this);
    inverse.setToIdentity();
    L.setToIdentity();
    // the elimination steps;
    /*
    
        1. Whenever the pivot evaluates to zero, swap it out with one of the rows beneath that has a non zero pivot
        2. Reduce all rows underneath by subtracting multiples of the pivot. A simple way would be to simply perform:
            Rp *= p*(mij/p);
    
    */
    int rows = this->nRows;
    int cols = this->nCols;

    for(int i = 0; i < rows; i++){
        // external loop evaulates all pivots
        // evaluate the status
        T p = (temp.data)[i*cols + i]; 
        if(p == static_cast<T>(0.0)){
            int index = temp.findPivot(i + 1, rows, i);
            if(index == -1){
                inverse.setToIdentity();
                this->isInvertible = false;
                return inverse;
            }
            temp.swapRows(i, index);
            inverse.swapRows(i, index);
            p = (temp.data)[i*cols + i]; 
        }
        for(int j = i + 1; j < rows; j++){
            T multiplier = ((temp.data)[j*cols + i])/p;
            for(int k = 0; k <= i; k++){
                // create the Linverse
                (L.data)[j*cols + k] -= multiplier*((L.data)[i*cols + k]);
            } 
            for(int k = 0; k < cols; k++){
                (temp.data)[j*cols + k] -= ((temp.data)[i*cols + k])*multiplier;
                if(abs((temp.data)[j*cols+ k]) < static_cast<T>(1e-9)) (temp.data)[j*cols + k] = static_cast<T>(0.0); // resolve floating point errors
                (inverse.data)[j*cols + k] -= ((inverse.data)[i*cols + k])*multiplier;
                if(abs((inverse.data)[j*cols+ k]) < static_cast<T>(1e-9)) (inverse.data)[j*cols + k] = static_cast<T>(0.0); // resolve floating point errors
            }
        }
        
        L.printUltimatrix();
        // At this point standard Gaussian elimination is complete, and the resultant temp Ultimatrix can be used for solution of linear systems using back subst
    }

    // Now let us reduce the upper triangular Ultimatrix to the identity Ultimatrix 
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < i; j++){
            (L.data)[i*cols + j] *= static_cast<T>(-1.00);
        }
    }

    // the LU decomposition is now complete
    Ultimatrix<T> U(temp);
    for(int i = rows-1; i >= 0; i--){
        T p = (temp.data)[i*cols + i];
        // we don't need to worry about zero pivots, forward elimination would have spotted singularities
        for(int j = i - 1; j >= 0; j--){
            T multiplier = (temp.data)[j*cols + i]/p;
            (temp.data)[j*cols + i] = static_cast<T>(0.0);
            for(int k = cols-1; k >= 0; k--){
                (inverse.data)[j*cols + k] -= ((inverse.data)[i*cols + k])*multiplier;
                if(abs((inverse.data)[j*cols+ k]) < static_cast<T>(1e-9)) {
                    (inverse.data)[j*cols + k] = static_cast<T>(0.0); 
                }// resolve floating point errors
            }
        }
    }

    // Normalize pivots

    for(int i = 0; i < rows; i++){
        (inverse.data)[i*cols + i] /= (temp.data)[i*cols + i];
        (temp.data)[i*cols + i] = 1;
    }

    // At this point the Ultimatrix has been successfully inverted.
    L.printUltimatrix();
    U.printUltimatrix();

    return inverse;
}

template <class T>
void Ultimatrix<T>::scaleRow(size_t row, T& mult){
    // function to scale a row by a certain amount 
    int cols = this->nCols;
    for(int i = 0; i < cols; i++){
        (this->data)[row*cols + i] *= mult;
    }
    return;
}










