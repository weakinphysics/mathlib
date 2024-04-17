#ifndef LINALG
#define LINALG
#endif


#include "./Ultimatrix.h"
#include "./Vector.h"

template <class T>
class LinearSystem{

    // designed with both real and complex systems in mind

    Ultimatrix <T> matrix;
    Vector <T> output;
    Ultimatrix<T> columnSpaceBasis;
    Ultimatrix<T> nullSpaceBasis;
    Ultimatrix<T> leftNullSpaceBasis;
    Ultimatrix<T> rowSpaceBasis;
    Ultimatrix<T> reduced;
    Ultimatrix<T> inverse;
    Ultimatrix<T> orthonormalBasis;

    // constructors

    LinearSystem();
    LinearSystem(size_t rows, size_t cols);
    LinearSystem(size_t rows, size_t cols, T* data, T* op);
    LinearSystem(Ultimatrix<T> &associated, Vector <T>& op);
    LinearSystem(Ultimatrix<T> associated, Vector<T> op);
};

