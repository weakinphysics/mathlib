#include "./Ultimatrix.h"
#include <bits/stdc++.h>
#include <math.h>

#ifndef xlr8
#define xlr8
#endif

/*

    Library: Vector 
    Author: Unknown
    Created: Some time in the early 21st century
    The name of the associated matrix library is ultimatrix, and thus I found it fitting to name this vector addendum as xlr8. Why? Vectors are commonly
    associated with locations, bearings/attitudes and velocities.
    The name of the "speed" alien in Ben10(Ultimatrix much?) is xlr8. Now you do the math(both figuratively and very very literally).

*/
template <class T>
class Vector: private Ultimatrix{
    // LEMUR
    Vector(); // default constructor
    Vector(size_t dim);
    Vector(size_t dim, T value); 
    Vector(size_t dim, T* inputvector);
    Vector(size_t d1, size_t d2); // constructor blocks illegal calls to base class
    Vector(size_t d1, size_t d2, T value);
    Vector(size_t d1, size_t d2, T* inputvector);
    Vector(Vector<T> &v);
    Vector(Vector<T>* v); // dont feed matrices;

    ~Vector(); // destructor

    // addition, subtraction and multiplication are covered in matrix mult
    // element wise multiplication is also covered by the ultimatrix implementation


    template<class U> friend Vector<U> operator=(Vector)
    template<class U> friend U operator*(Vector<U> &v1, Vector<U>& v2); // evaluate vector dot product 
    template<class U> friend U operator*(U factor, Vector<U> &v); // scale by a constant factor
    template<class U> friend Vector<U> cross(Vector<U> &lhs, Vector <U> &rhs); // evaluate vector cross product
    template<class U> friend Vector<U> goldsteinRotation(Vector<U> &axis, Vector <U> &vector, U theta); // evaluate rotation along a given axis
    template<class U> friend Vector<U> project(Vector<U> &target, Vector<U> &space);
    
    
    T magnitude();
    void scale(T value);
};

/*
    IMPLEMENTATION
*/
template<class T>
Vector<T>::Vector(){
    std::cout<<"Vector constructor called"<<std::endl;
    this->nRows = 1;
    this->nCols = 1;
    this->allocatedMemory = 1;
    this->isVector = true;
    this->data = new T[1];
    (this->data)[0] = static_cast<T>(0.0);
}

template<class T>
Vector<T>::Vector(size_t dim){
    std::cout<<"Vector constructor called"<<std::endl;
    this->nRows = dim;
    this->nCols = 1;
    this->nElements = dim;
    this->allocatedMemory = dim;
    this->isVector = true;
    this->data = new T[dim];
    for(int i = 0; i < dim; i++) (this->data)[i] = static_cast<T>(0.0);
}

template<class T>
Vector<T>::Vector(size_t dim, T value){
    std::cout<<"Vector constructor called"<<std::endl;
    this->nRows = dim;
    this->nCols = 1;
    this->nElements = dim;
    this->allocatedMemory = dim;
    this->isVector = true;
    this->data = new T[dim];
    for(int i = 0; i < dim; i++) (this->data)[i] = value;
}

template<class T>
Vector<T>::Vector(size_t dim, T* input){
    std::cout<<"Vector constructor called"<<std::endl;
    this->nRows = dim;
    this->nCols = 1;
    this->nElements = dim;
    this->allocatedMemory = dim;
    this->isVector = true;
    this->data = new T[dim];
    for(int i = 0; i < dim; i++) (this->data)[i] = input[i];
}

template<class T>
Vector<T>::Vector(Vector<T>& v){
    std::cout<<"Vector constructor called"<<std::endl;
    this->nRows = v.nRows;
    this->nCols = 1;
    this->nElements = v.nRows;
    this->allocatedMemory = v.nRows;
    this->isVector = true;
    this->data = new T[v.nRows];
    for(int i = 0; i < v.nRows; i++) (this->data)[i] = (v.data)[i];
}

template<class T>
Vector<T>::Vector(Vector<T>* v){
    std::cout<<"Vector constructor called"<<std::endl;
    this->nRows = v->nRows;
    this->nCols = 1;
    this->nElements = v->nRows;
    this->allocatedMemory = v->nRows;
    this->isVector = true;
    this->data = new T[v->nRows];
    for(int i = 0; i < v->nRows; i++) (this->data)[i] = (v->data)[i];
}


template <class T>
T operator*(Vector<T> &V, Vector <T> &Y){
    //evaluates v.y
    assert(V.nRows == Y.nCols);
    T acc = static_cast<T>(0.0);
    for(size_t i = 0; i < V.nRows; i++){
        acc += ((V.data)[i])*((Y.data)[i]);
    }
    return acc;
}

template <class T>
Vector <T> cross(Vector <T> &a, Vector <T> &b){
    assert(a.nRows == b.nRows);
    assert(a.nRows == 3);
    T c1 = ((a.data)[1])*((b.data)[2]) - ((a.data)[2])*((b.data)[1]);
    T c2 = ((a.data)[2])*((b.data)[0]) - ((a.data)[0])*((b.data)[2]);
    T c3 = ((a.data)[0])*((b.data)[1]) - ((a.data)[1])*((b.data)[0]);
    Vector <T> product(a.nRows);
    (product.data)[0] = c1;
    (product.data)[1] = c2;
    (product.data)[2] = c3;    
    return product;
}


template <class T>
Vector<T> project(Vector<T> &v, Vector<T> &y){
    // uses (vTy/yTy)*a
    // projects v onto y
    T mag = y.magnitude();
    T dot = v*y;
    T coefficient = mag*dot;
    Vector<T> projection(y.nRows); 
    for(size_t i = 0; i < y.nRows; i++){
        (projection.data)[i] = coefficient*((y.data)[i]);
    }
}


template <class T>
T Vector<T>::magnitude(){
    T acc = static_cast<T> (0.0);
    for(size_t i = 0; i < this->nRows; i++){
        acc += ((this->data)[i])*((this->data)[i]);
    }
    return acc;
}

template <class T>
void Vector<T>::scale(T value){
    for(size_t i = 0; i < this->nRows; i++){
        (this->data)[i] *= value;
    }
    return;
}


template <class T>
Vector <T> goldsteinRotation(Vector <T>& axis, Vector <T>& vector, T theta){
    // make sure the axis and vector are of the same dimension 
    assert(axis.nRows == vector.nRows);
    size_t rowCount = axis.nRows;
    T coeff1 = math.cos(theta);  
    T coeff2 = math.sin(theta); 
    T coeff3 = ((axis*vector)/axis.magnitude());
    Vector <T> auxilliary = project(vector, axis);
    Vector <T> error = vector - auxilliary;
    Vector <T> crossProduct = cross(auxilliary, error);
    error.scale(coeff1);
    crossProduct.scale(coeff2);
    Vector <T> ans = error + crossProduct;
    return ans;
}







// unfamiliar with my own matrix library!
// the face inside is right beneath my skin 