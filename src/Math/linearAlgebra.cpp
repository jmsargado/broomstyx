/*
  This file is part of the BROOMStyx project.

  BROOMStyx is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  BROOMStyx is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with BROOMStyx.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the AUTHORS file
  for the list of copyright holders.
*/

#include "linearAlgebra.hpp"
#include <cstdio>
#include <stdexcept>

#include <config.h>

#ifdef HAVE_MKL
    #include "mkl_cblas.h"
    #include "mkl_lapacke.h"
#else
    #include "cblas.h"
    #include "lapacke.h"
#endif

namespace broomstyx
{
    // Matrix addition
    // ---------------------------------------------------------------------------
    RealMatrix operator+( RealMatrix&& A, RealMatrix&& B ) 
    {
        int A_dim1 = A.dim1();
        int A_dim2 = A.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();

#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error( "\nSize mismatch in operands for matrix addition!\n\tdim(A) = [ "
                    + std::to_string( A_dim1 ) + " x " + std::to_string( A_dim2 ) + " ], dim(B) = [ "
                    + std::to_string( B_dim1 ) + " x " + std::to_string( B_dim2 ) + " ]" );
#endif
        cblas_daxpy( A_dim1 * A_dim2, 1.0, B.ptr(), 1, A.ptr(), 1 );
        return (RealMatrix&&)A;
    }

    RealMatrix operator+( RealMatrix&& A, const RealMatrix& B ) 
    {
        int A_dim1 = A.dim1();
        int A_dim2 = A.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error( "\nSize mismatch in operands for matrix addition!\n\tdim(A) = [ "
                    + std::to_string( A_dim1 ) + " x " + std::to_string( A_dim2 ) + " ], dim(B) = [ "
                    + std::to_string( B_dim1 ) + " x " + std::to_string( B_dim2 ) + " ]" );
#endif
        cblas_daxpy( A_dim1 * A_dim2, 1.0, B.ptr(), 1, A.ptr(), 1 );
        return (RealMatrix&&)A;
    }

    RealMatrix operator+( const RealMatrix& A, RealMatrix&& B ) 
    {
        int A_dim1 = A.dim1();
        int A_dim2 = A.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error( "\nSize mismatch in operands for matrix addition!\n\tdim(A) = [ "
                    + std::to_string( A_dim1 ) + " x " + std::to_string( A_dim2 ) + " ], dim(B) = [ "
                    + std::to_string( B_dim1 ) + " x " + std::to_string( B_dim2 ) + " ]" );
#endif
        cblas_daxpy( B_dim1 * B_dim2, 1.0, A.ptr(), 1, B.ptr(), 1 );

        return (RealMatrix&&)B;
    }

    RealMatrix operator+( const RealMatrix& A, const RealMatrix& B )
    {
        int A_dim1 = A.dim1();
        int A_dim2 = A.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error( "\nSize mismatch in operands for matrix addition!\n\tdim(A) = [ "
                    + std::to_string( A_dim1 ) + " x " + std::to_string( A_dim2 ) + " ], dim(B) = [ "
                    + std::to_string( B_dim1 ) + " x " + std::to_string( B_dim2 ) + " ]" );
#endif
        RealMatrix C = A;
        cblas_daxpy( A_dim1 * A_dim2, 1.0, B.ptr(), 1, C.ptr(), 1 );
        
        return C;
    }

    // Vector addition
    // ---------------------------------------------------------------------------
    RealVector operator+( RealVector&& A, RealVector&& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error( "\nSize mismatch in operands for vector addition!\n\tdim(A) = "
                + std::to_string( A_dim ) + ", dim(B) = " + std::to_string( B_dim ) );
#endif
        cblas_daxpy( A_dim, 1.0, B.ptr(), 1, A.ptr(), 1 );

        return (RealVector&&)A;
    }

    RealVector operator+( RealVector&& A, const RealVector& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error( "\nSize mismatch in operands for vector addition!\n\tdim(A) = "
                + std::to_string( A_dim ) + ", dim(B) = " + std::to_string( B_dim ) );
#endif
        cblas_daxpy( A_dim, 1.0, B.ptr(), 1, A.ptr(), 1 );

        return (RealVector&&)A;
    }

    RealVector operator+( const RealVector& A, RealVector&& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error( "\nSize mismatch in operands for vector addition!\n\tdim(A) = "
                + std::to_string( A_dim ) + ", dim(B) = " + std::to_string( B_dim ) );
#endif
        cblas_daxpy( B_dim, 1.0, A.ptr(), 1, B.ptr(), 1 );

        return (RealVector&&)B;
    }

    RealVector operator+( const RealVector& A, const RealVector& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error( "\nSize mismatch in operands for vector addition!\n\tdim(A) = "
                + std::to_string( A_dim ) + ", dim(B) = " + std::to_string( B_dim ) );
#endif
        RealVector C = A;
        cblas_daxpy( A_dim, 1.0, B.ptr(), 1, C.ptr(), 1 );

        return C;
    }

    // Matrix subtraction
    // ---------------------------------------------------------------------------
    RealMatrix operator-( RealMatrix&& A, RealMatrix&& B )
    {
        int A_dim1 = A.dim1();
        int A_dim2 = A.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error( "\nSize mismatch in operands for matrix subctraction!\n\tdim(A) = [ "
                    + std::to_string( A_dim1 ) + " x " + std::to_string( A_dim2 ) + " ], dim(B) = [ "
                    + std::to_string( B_dim1 ) + " x " + std::to_string( B_dim2 ) + " ]" );
#endif
        cblas_daxpy( A_dim1 * A_dim2, -1.0, B.ptr(), 1, A.ptr(), 1 );

        return (RealMatrix&&)A;
    }

    RealMatrix operator-( RealMatrix&& A, const RealMatrix& B )
    {
        int A_dim1 = A.dim1();
        int A_dim2 = A.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error( "\nSize mismatch in operands for matrix subctraction!\n\tdim(A) = [ "
                    + std::to_string( A_dim1 ) + " x " + std::to_string( A_dim2 ) + " ], dim(B) = [ "
                    + std::to_string( B_dim1 ) + " x " + std::to_string( B_dim2 ) + " ]" );
#endif
        cblas_daxpy( A_dim1 * A_dim2, -1.0, B.ptr(), 1, A.ptr(), 1 );

        return (RealMatrix&&)A;
    }

    RealMatrix operator-( const RealMatrix& A, RealMatrix&& B )
    {
        RealMatrix C = A;

        int A_dim1 = C.dim1();
        int A_dim2 = C.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error( "\nSize mismatch in operands for matrix subctraction!\n\tdim(A) = [ "
                    + std::to_string( A_dim1 ) + " x " + std::to_string( A_dim2 ) + " ], dim(B) = [ "
                    + std::to_string( B_dim1 ) + " x " + std::to_string( B_dim2 ) + " ]" );
#endif
        cblas_daxpy( A_dim1 * A_dim2, -1.0, B.ptr(), 1, C.ptr(), 1 );

        return C;
    }

    RealMatrix operator-( const RealMatrix& A, const RealMatrix& B )
    {
        RealMatrix C = A;

        int A_dim1 = C.dim1();
        int A_dim2 = C.dim2();
        int B_dim1 = B.dim1();
        int B_dim2 = B.dim2();
#ifndef NDEBUG
        if ( A_dim1 != B_dim1 || A_dim2 != B_dim2 )
            throw std::runtime_error( "\nSize mismatch in operands for matrix subctraction!\n\tdim(A) = [ "
                    + std::to_string( A_dim1 ) + " x " + std::to_string( A_dim2 ) + " ], dim(B) = [ "
                    + std::to_string( B_dim1 ) + " x " + std::to_string( B_dim2 ) + " ]" );
#endif
            cblas_daxpy( A_dim1 * A_dim2, -1.0, B.ptr(), 1, C.ptr(), 1 );

        return C;
    }

    // Vector subtraction
    // ---------------------------------------------------------------------------
    RealVector operator-( RealVector&& A, RealVector&& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error( "\nSize mismatch in operands for vector subtraction!\n\tdim(A) = "
                + std::to_string( A_dim ) + ", dim(B) = " + std::to_string( B_dim ) );
#endif
        cblas_daxpy( A_dim, -1.0, B.ptr(), 1, A.ptr(), 1 );

        return (RealVector&&)A;
    }

    RealVector operator-( RealVector&& A, const RealVector& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error( "\nSize mismatch in operands for vector subtraction!\n\tdim(A) = "
                + std::to_string( A_dim ) + ", dim(B) = " + std::to_string( B_dim ) );
#endif
        cblas_daxpy( A_dim, -1.0, B.ptr(), 1, A.ptr(), 1 );

        return (RealVector&&)A;
    }

    RealVector operator-( const RealVector& A, RealVector&& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error( "\nSize mismatch in operands for vector subtraction!\n\tdim(A) = "
                + std::to_string( A_dim ) + ", dim(B) = " + std::to_string( B_dim ) );
#endif
        RealVector C = A;
        cblas_daxpy( A_dim, -1.0, B.ptr(), 1, C.ptr(), 1 );

        return C;
    }

    RealVector operator-( const RealVector& A, const RealVector& B )
    {
        int A_dim = A.dim();
        int B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim != B_dim )
            throw std::runtime_error( "\nSize mismatch in operands for vector addition!\n\tdim(A) = "
                + std::to_string( A_dim ) + ", dim(B) = " + std::to_string( B_dim ) );
#endif
        RealVector C = A;
        cblas_daxpy( A_dim, -1.0, B.ptr(), 1, C.ptr(), 1 );

        return C;
    }

    // Scalar multiplication for matrices
    // ---------------------------------------------------------------------------
    RealMatrix operator*( RealMatrix&& A, double b )
    {
        int dim1 = A.dim1();
        int dim2 = A.dim2();
        double* ptr = A.ptr();

        cblas_dscal( dim1 * dim2, b, ptr, 1 );

        return (RealMatrix&&)A;
    }

    RealMatrix operator*( const RealMatrix& A, double b )
    {
        RealMatrix C = A;
        int dim1 = C.dim1();
        int dim2 = C.dim2();
        double* ptr = C.ptr();

        cblas_dscal( dim1 * dim2, b, ptr, 1 );

        return C;
    }

    RealMatrix operator*( double a, RealMatrix&& B )
    {
        int dim1 = B.dim1();
        int dim2 = B.dim2();
        double* ptr = B.ptr();

        cblas_dscal( dim1 * dim2, a, ptr, 1 );

        return (RealMatrix&&)B;
    }

    RealMatrix operator*( double a, const RealMatrix& B )
    {
        RealMatrix C = B;
        int dim1 = C.dim1();
        int dim2 = C.dim2();
        double* ptr = C.ptr();

        cblas_dscal( dim1 * dim2, a, ptr, 1 );

        return C;
    }

    // Scalar division for matrices
    // ---------------------------------------------------------------------------
    RealMatrix operator/( RealMatrix&& A, double b )
    {
        int dim1 = A.dim1();
        int dim2 = A.dim2();
        double* ptr = A.ptr();

        cblas_dscal( dim1 * dim2, 1./b, ptr, 1 );

        return (RealMatrix&&)A;
    }

    RealMatrix operator/( const RealMatrix& A, double b )
    {
        RealMatrix C = A;
        int dim1 = C.dim1();
        int dim2 = C.dim2();
        double* ptr = C.ptr();

        cblas_dscal( dim1 * dim2, 1./b, ptr, 1 );

        return C;
    }

    // Scalar multiplication for vectors
    // ---------------------------------------------------------------------------
    RealVector operator*( RealVector&& A, double b )
    {
        int dim = A.dim();
        double* ptr = A.ptr();

        cblas_dscal( dim, b, ptr, 1 );

        return (RealVector&&)A;
    }

    RealVector operator*( const RealVector& A, double b )
    {
        RealVector C = A;
        int dim = C.dim();
        double* ptr = C.ptr();

        cblas_dscal( dim, b, ptr, 1 );

        return C;
    }

    RealVector operator*( double a, RealVector&& B )
    {
        int dim = B.dim();
        double* ptr = B.ptr();

        cblas_dscal( dim, a, ptr, 1 );

        return (RealVector&&)B;
    }

    RealVector operator*( double a, const RealVector& B )
    {
        RealVector C = B;
        int dim = C.dim();
        double* ptr = C.ptr();

        cblas_dscal( dim, a, ptr, 1 );

        return C;
    }

    // Scalar division for vectors
    // ---------------------------------------------------------------------------
    RealVector operator/( RealVector&& A, double b )
    {
        int dim = A.dim();
        double* ptr = A.ptr();

        cblas_dscal( dim, 1./b, ptr, 1 );

        return (RealVector&&)A;
    }

    RealVector operator/( const RealVector& A, double b )
    {
        RealVector C = A;
        int dim = C.dim();
        double* ptr = C.ptr();

        cblas_dscal( dim, 1./b, ptr, 1 );

        return C;
    }

    // Matrix multiplication
    // ---------------------------------------------------------------------------
    RealMatrix operator*( RealMatrix&& A, RealMatrix&& B )
    {
        int opA_dim1, opA_dim2, opB_dim1, opB_dim2;

        int ldA = A.dim1();
        int ldB = B.dim1();

        opA_dim1 = A.dim1();
        opA_dim2 = A.dim2();

        opB_dim1 = B.dim1();
        opB_dim2 = B.dim2();

#ifndef NDEBUG
        if ( opA_dim2 != opB_dim1 )
            throw std::runtime_error( "\nSize mismatch in operands for matrix multiplication!\n\tdim(opA) = [ "
                + std::to_string(opA_dim1) + " x " + std::to_string(opA_dim2) + " ], dim(opB) = [ " + std::to_string(opB_dim1)
                + " x " + std::to_string(opB_dim2) + " ]" );
#endif
        RealMatrix C( opA_dim1, opB_dim2 );
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, opA_dim1, opB_dim2, opA_dim2, 1.0, A.ptr(), ldA, B.ptr(), ldB, 0.0, C.ptr(), opA_dim1 );

        return C;
    }

    RealMatrix operator*( RealMatrix&& A, const RealMatrix& B )
    {
        int opA_dim1, opA_dim2, opB_dim1, opB_dim2;

        int ldA = A.dim1();
        int ldB = B.dim1();

        opA_dim1 = A.dim1();
        opA_dim2 = A.dim2();

        opB_dim1 = B.dim1();
        opB_dim2 = B.dim2();

#ifndef NDEBUG
        if ( opA_dim2 != opB_dim1 )
            throw std::runtime_error( "\nSize mismatch in operands for matrix multiplication!\n\tdim(opA) = [ "
                + std::to_string(opA_dim1) + " x " + std::to_string(opA_dim2) + " ], dim(opB) = [ " + std::to_string(opB_dim1)
                + " x " + std::to_string(opB_dim2) + " ]" );
#endif
        RealMatrix C( opA_dim1, opB_dim2 );
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, opA_dim1, opB_dim2, opA_dim2, 1.0, A.ptr(), ldA, B.ptr(), ldB, 0.0, C.ptr(), opA_dim1 );

        return C;
    }

    RealMatrix operator*( const RealMatrix& A, RealMatrix&& B )
    {
        int opA_dim1, opA_dim2, opB_dim1, opB_dim2;

        int ldA = A.dim1();
        int ldB = B.dim1();

        opA_dim1 = A.dim1();
        opA_dim2 = A.dim2();

        opB_dim1 = B.dim1();
        opB_dim2 = B.dim2();

#ifndef NDEBUG
        if ( opA_dim2 != opB_dim1 )
            throw std::runtime_error( "\nSize mismatch in operands for matrix multiplication!\n\tdim(opA) = [ "
                + std::to_string(opA_dim1) + " x " + std::to_string(opA_dim2) + " ], dim(opB) = [ " + std::to_string(opB_dim1)
                + " x " + std::to_string(opB_dim2) + " ]" );
#endif
        RealMatrix C( opA_dim1, opB_dim2 );
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, opA_dim1, opB_dim2, opA_dim2, 1.0, A.ptr(), ldA, B.ptr(), ldB, 0.0, C.ptr(), opA_dim1 );

        return C;
    }

    RealMatrix operator*( const RealMatrix& A, const RealMatrix& B )
    {
        int opA_dim1, opA_dim2, opB_dim1, opB_dim2;

        int ldA = A.dim1();
        int ldB = B.dim1();

        opA_dim1 = A.dim1();
        opA_dim2 = A.dim2();

        opB_dim1 = B.dim1();
        opB_dim2 = B.dim2();

#ifndef NDEBUG
        if ( opA_dim2 != opB_dim1 )
            throw std::runtime_error( "\nSize mismatch in operands for matrix multiplication!\n\tdim(opA) = [ "
                + std::to_string(opA_dim1) + " x " + std::to_string(opA_dim2) + " ], dim(opB) = [ " + std::to_string(opB_dim1)
                + " x " + std::to_string(opB_dim2) + " ]" );
#endif
        RealMatrix C( opA_dim1, opB_dim2 );
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, opA_dim1, opB_dim2, opA_dim2, 1.0, A.ptr(), ldA, B.ptr(), ldB, 0.0, C.ptr(), opA_dim1 );

        return C;
    }

    // Matrix-vector multiplication
    // ---------------------------------------------------------------------------
    RealVector operator*( RealMatrix&& A, RealVector&& B )
    {
        int A_dim1, A_dim2, B_dim;

        A_dim1 = A.dim1();
        A_dim2 = A.dim2();

        int ldA = A_dim1;
        B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim2 != B_dim )
            throw std::runtime_error( "\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(A) = [ "
                + std::to_string( A_dim1 ) + " x " + std::to_string( A_dim2 ) + " ], dim(B) =" + std::to_string( B_dim ) );
#endif
        RealVector C( A_dim1 );
        cblas_dgemv( CblasColMajor, CblasNoTrans, A_dim1, A_dim2, 1.0, A.ptr(), ldA, B.ptr(), 1, 1.0, C.ptr(), 1 );

        return C;
    }

    RealVector operator*( RealMatrix&& A, const RealVector& B )
    {
        int A_dim1, A_dim2, B_dim;

        A_dim1 = A.dim1();
        A_dim2 = A.dim2();

        int ldA = A_dim1;
        B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim2 != B_dim )
            throw std::runtime_error( "\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(A) = [ "
                + std::to_string( A_dim1 ) + " x " + std::to_string( A_dim2 ) + " ], dim(B) =" + std::to_string( B_dim ) );
#endif
        RealVector C( A_dim1 );
        cblas_dgemv( CblasColMajor, CblasNoTrans, A_dim1, A_dim2, 1.0, A.ptr(), ldA, B.ptr(), 1, 1.0, C.ptr(), 1 );

        return C;
    }

    RealVector operator*( const RealMatrix& A, RealVector&& B )
    {
        int A_dim1, A_dim2, B_dim;

        A_dim1 = A.dim1();
        A_dim2 = A.dim2();

        int ldA = A_dim1;
        B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim2 != B_dim )
            throw std::runtime_error( "\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(A) = [ "
                + std::to_string( A_dim1 ) + " x " + std::to_string( A_dim2 ) + " ], dim(B) =" + std::to_string( B_dim ) );
#endif
        RealVector C( A_dim1 );
        cblas_dgemv( CblasColMajor, CblasNoTrans, A_dim1, A_dim2, 1.0, A.ptr(), ldA, B.ptr(), 1, 1.0, C.ptr(), 1 );

        return C;
    }

    RealVector operator*( const RealMatrix& A, const RealVector& B )
    {
        int A_dim1, A_dim2, B_dim;

        A_dim1 = A.dim1();
        A_dim2 = A.dim2();

        int ldA = A_dim1;
        B_dim = B.dim();
#ifndef NDEBUG
        if ( A_dim2 != B_dim )
            throw std::runtime_error( "\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(A) = [ "
                + std::to_string( A_dim1 ) + " x " + std::to_string( A_dim2 ) + " ], dim(B) =" + std::to_string( B_dim ) );
#endif
        RealVector C( A_dim1 );
        cblas_dgemv( CblasColMajor, CblasNoTrans, A_dim1, A_dim2, 1.0, A.ptr(), ldA, B.ptr(), 1, 1.0, C.ptr(), 1 );

        return C;
    }

    // Vector-matrix multiplication
    // ---------------------------------------------------------------------------
    RealVector operator*( RealVector&& A, RealMatrix&& B )
    {
        int B_dim1, B_dim2, A_dim;

        B_dim1 = B.dim1();
        B_dim2 = B.dim2();

        int ldB = B_dim1;
        A_dim = A.dim();
#ifndef NDEBUG
        if ( B_dim1 != A_dim )
            throw std::runtime_error( "\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(B) = [ "
                + std::to_string( B_dim1 ) + " x " + std::to_string( B_dim2 ) + " ], dim(A) = " + std::to_string( A_dim ) );
#endif
        RealVector C( B_dim2 );
        cblas_dgemv( CblasColMajor, CblasTrans, B_dim1, B_dim2, 1.0, B.ptr(), ldB, A.ptr(), 1, 1.0, C.ptr(), 1 );

        return C;
    }

    RealVector operator*( RealVector&& A, const RealMatrix& B )
    {
        int B_dim1, B_dim2, A_dim;

        B_dim1 = B.dim1();
        B_dim2 = B.dim2();

        int ldB = B_dim1;
        A_dim = A.dim();
#ifndef NDEBUG
        if ( B_dim1 != A_dim )
            throw std::runtime_error( "\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(B) = [ "
                + std::to_string( B_dim1 ) + " x " + std::to_string( B_dim2 ) + " ], dim(A) = " + std::to_string( A_dim ) );
#endif
        RealVector C( B_dim2 );
        cblas_dgemv( CblasColMajor, CblasTrans, B_dim1, B_dim2, 1.0, B.ptr(), ldB, A.ptr(), 1, 1.0, C.ptr(), 1 );

        return C;
    }

    RealVector operator*( const RealVector& A, RealMatrix&& B )
    {
        int B_dim1, B_dim2, A_dim;

        B_dim1 = B.dim1();
        B_dim2 = B.dim2();

        int ldB = B_dim1;
        A_dim = A.dim();
#ifndef NDEBUG
        if ( B_dim1 != A_dim )
            throw std::runtime_error( "\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(B) = [ "
                + std::to_string( B_dim1 ) + " x " + std::to_string( B_dim2 ) + " ], dim(A) = " + std::to_string( A_dim ) );
#endif
        RealVector C( B_dim2 );
        cblas_dgemv( CblasColMajor, CblasTrans, B_dim1, B_dim2, 1.0, B.ptr(), ldB, A.ptr(), 1, 1.0, C.ptr(), 1 );

        return C;
    }

    RealVector operator*( const RealVector& A, const RealMatrix& B )
    {
        int B_dim1, B_dim2, A_dim;

        B_dim1 = B.dim1();
        B_dim2 = B.dim2();

        int ldB = B_dim1;
        A_dim = A.dim();
#ifndef NDEBUG
        if ( B_dim1 != A_dim )
            throw std::runtime_error( "\nSize mismatch in operands for matrix-vector multiplication!\n\tdim(B) = [ "
                + std::to_string( B_dim1 ) + " x " + std::to_string( B_dim2 ) + " ], dim(A) = " + std::to_string( A_dim ) );
#endif
        RealVector C( B_dim2 );
        cblas_dgemv( CblasColMajor, CblasTrans, B_dim1, B_dim2, 1.0, B.ptr(), ldB, A.ptr(), 1, 1.0, C.ptr(), 1 );

        return C;
    }

    // Matrix transpose
    // ---------------------------------------------------------------------------
    RealMatrix trp( RealMatrix&& A )
    {
        return A.trp();
    }

    RealMatrix trp( const RealMatrix& A )
    {
        return A.trp();
    }

    // Matrix inverse
    // ------------------------------------------------------------------------
    RealMatrix inv( RealMatrix&& A )
    {
        int info;
        lapack_int* ipiv;

        int nrows = A.dim1();
        int ncols = A.dim2();

        if ( ncols != nrows )
            throw std::runtime_error( "Cannot invert non-square matrix!\ndim(A) = [ " + std::to_string( nrows )
                + " x " + std::to_string( ncols ) + " ]" );

        ipiv = new lapack_int[ nrows ]();
        double *locA = A.ptr();

        info = LAPACKE_dgetrf( LAPACK_COL_MAJOR, nrows, nrows, locA, nrows, ipiv );
        if ( info < 0 )
        {
            delete[] ipiv;
            throw std::runtime_error( "Error in matrix inversion: illegal value encountered during factorization!" );
        }
        else if ( info > 0 )
        {
            delete[] ipiv;
            throw std::runtime_error( "Cannot invert singular matrix!" );
        }

        info = LAPACKE_dgetri( LAPACK_COL_MAJOR, nrows, locA, nrows, ipiv );

        delete[] ipiv;

        if ( info < 0 )
            throw std::runtime_error( "Error in matrix inversion: illegal value encountered during factorization!" );
        else if ( info > 0 )
            throw std::runtime_error( "Cannot invert singular matrix!" );

        return (RealMatrix&&)A;
    }

    RealMatrix inv( const RealMatrix& B )
    {
        RealMatrix A = B;

        int info;
        lapack_int* ipiv;
        int nrows = A.dim1();
        int ncols = A.dim2();

        if ( ncols != nrows )
            throw std::runtime_error( "Cannot invert non-square matrix!\ndim(A) = [ " + std::to_string( nrows )
                + " x " + std::to_string( ncols ) + " ]" );

        ipiv = new lapack_int[ nrows ]();
        double *locA = A.ptr();

        info = LAPACKE_dgetrf( LAPACK_COL_MAJOR, nrows, nrows, locA, nrows, ipiv );
        if ( info < 0 )
        {
            delete[] ipiv;
            throw std::runtime_error( "Error in matrix inversion: illegal value encountered during factorization!" );
        }
        else if ( info > 0 )
        {
            delete[] ipiv;
            throw std::runtime_error( "Cannot invert singular matrix!" );
        }

        info = LAPACKE_dgetri( LAPACK_COL_MAJOR, nrows, locA, nrows, ipiv );

        delete[] ipiv;

        if ( info < 0 )
            throw std::runtime_error( "Error in matrix inversion: illegal value encountered during factorization!" );
        else if ( info > 0 )
            throw std::runtime_error( "Cannot invert singular matrix!" );

        return A;
    }
}
