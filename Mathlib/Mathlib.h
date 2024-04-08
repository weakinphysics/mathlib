#ifndef MTLIB_C
#define MTLIB_C

#include <math.h>
#include <complex.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <pthread.h>
#include <bits/stdc++.h>
#include <stdlib.h>
#include <unistd.h>
#include <stddef.h>
#include <stdint.h>
#include <chrono.h>
#include <random.h>
#include <mmintrin.h>

#define MATRIX_ELEMENT(m, i, j) (*((m)->data + i*((m)->cols) + j))
#define PI 3.1415

typedef struct{
    double x;
} Number;


// class Matrix{
//     // This is for matrices. Not for other data structures such as tensors. That will come later.
//     public:     
//         double *data;   
//         size_t size;
//         size_t rows;
//         size_t cols;
//         Matrix(size_t r, size_t c){
//             assert(r != 0 && c != 0);
//             this->rows = r;
//             this->cols = c;
//             size_t memory = rows*cols; // please do not request absurd amounts of memory. this program is not protected against attempts to allocate ridiculously heavy memory slices
//             this->data = calloc(sizeof(double)*memory);
//             return;
//         }

//         Matrix(size_t r, size_t c, double *d){
//             assert(r != 0 && c != 0);
//             this->rows = r;
//             this->cols = c;
//             size_t memory = rows*cols; // please do not request absurd amounts of memory. this program is not protected against attempts to allocate ridiculously heavy memory slices
//             this->data = calloc(sizeof(double)*memory);
//             for(size_t i = 0; i < r; i++){
//                 for(size_t j = 0; j < c; j++) MATRIX_ELEMENT(this, i, j) = *(d + i*cols + j);
//             }
//         }
        
//         initializeMatrixToRandomValues(double randMin, double randMax);

//         initializeMatrix(double *input);
        
//         static Matrix multiplyMatrices(Matrix *a, Matrix *b){
//             assert(a->cols == b->rows);
//             size_t new_rows = a->rows;
//             size_t new_cols = b->cols;
//             Matrix response = new Matrix(new_rows, new_cols);
//             for(size_t i = 0; i < new_rows; i++){
//                 for(size_t j = 0; j < new_cols; j++){
//                     for(size_t k = 0; k < a->cols; k++){
//                         MATRIX_ELEMENT(response, i, j) += (MATRIX_ELEMENT(a, i, k)*MATRIX_ELEMENT(b, k, j))
//                     }
//                 }
//             }
//             return response; // this will create a new object of the type returned, which vexes me. The good news is the size of the newly generated object is minimal 

//         }

//         static Matrix rref(Matrix *a){
//             // write code for reduced row echelon form
//         }

//         static Matrix ldu(Matrix *a){
//             // write code for LDU decomposition
//         }

//         static Matrix nullspace(Matrix *a){
//             // code evaluates and returns the nullspace of a matrix
//         }

//         static Matrix leftNullspace(Matrix *a){
//             // code evaluates and returns the left nullspace of a matrix
//         }

// }




