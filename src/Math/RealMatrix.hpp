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

#ifndef REALMATRIX_HPP
#define	REALMATRIX_HPP

#include <config.h>

#include <cstdio>
#include <cassert>
#include <stdexcept>
#include <initializer_list>

#ifdef HAVE_MKL
    #include "mkl_cblas.h"
#else
    #include "cblas.h"
#endif

namespace broomstyx
{
    class RealMatrix final {
        friend RealMatrix operator*( RealMatrix&& A, double b );
        friend RealMatrix operator*( double a, RealMatrix&& B );
        friend RealMatrix operator/( RealMatrix&& A, double b );
        friend RealMatrix trp( RealMatrix&& A );
        
    public:
        // Default constructor
        RealMatrix()
            : _dim1( 0 )
            , _dim2( 0 )
            , _ptr( nullptr )
        {}

        // Constructor with initial memory allocation
        RealMatrix( int dim1, int dim2 )
        {
            assert( dim1 >= 1 );
            assert( dim2 >= 1 );

            _dim1 = dim1;
            _dim2 = dim2;
            _ptr = new double[ _dim1 * _dim2 ]();
        }
        
        // Constructor with initializer list
        RealMatrix( std::initializer_list< std::initializer_list<double>> initList )
        {
            _dim1 = (int)initList.size();
            _dim2 = (int)( *( initList.begin() ) ).size();
#ifndef NDEBUG
            for ( int i = 1; i < _dim1; i++ )
                assert( (int)( *(initList.begin() + i) ).size() == _dim2 );
#endif
            _ptr = new double[ _dim1 * _dim2 ];
            
            for ( int i = 0; i < _dim1; i++ )
                for ( int j = 0; j < _dim2; j++ )
                    _ptr[ j*_dim1 + i ] = *( ( *(initList.begin() + i) ).begin() + j );
        }
        
        // Copy constructor
        RealMatrix( const RealMatrix& source )
        {
            if ( source._dim1 > 0 || source._dim2 > 0 )
            {
                _dim1 = source._dim1;
                _dim2 = source._dim2;
                _ptr = new double[ _dim1 * _dim2 ];
#ifdef VERBOSE_REALMATRIX_CONSTRUCTION
                std::printf( "...RealMatrix Copy constructor called.\n" );
                std::fflush(stdout);
#endif
                std::copy( source._ptr, source._ptr + _dim1 * _dim2, _ptr );
            }
            else
            {
                _dim1 = 0;
                _dim2 = 0;
                _ptr = nullptr;
            }
        }

        // Move constructor
        RealMatrix( RealMatrix&& source )
            : _dim1( source._dim1 )
            , _dim2( source._dim2 )
            , _ptr( source._ptr )
        {
#ifdef VERBOSE_REALMATRIX_CONSTRUCTION
            std::printf("...RealMatrix Move constructor called.\n");
            std::fflush(stdout);
#endif
            source._dim1 = 0;
            source._dim2 = 0;
            source._ptr = nullptr;
        }

        // Destructor
        virtual ~RealMatrix()
        {
            if ( _ptr )
                delete[] _ptr;
        }

        // Assignment with initializer list
        RealMatrix& operator= ( std::initializer_list< std::initializer_list<double>> initList )
        {
            _dim1 = (int)initList.size();
            _dim2 = (int)( *(initList.begin()) ).size();
#ifndef NDEBUG
            for ( int i = 1; i < _dim1; i++ )
                assert( (int)(*(initList.begin() + i)).size() == _dim2 );
#endif
            if ( _ptr )
                delete[] _ptr;
                
            _ptr = new double[ _dim1 * _dim2 ];
            
            for ( int i = 0; i < _dim1; i++ )
                for ( int j = 0; j < _dim2; j++ )
                    _ptr[ j*_dim1 + i ] = *( ( *(initList.begin() + i) ).begin() + j );

            return *this;
        }
        
        // Copy assignment operator
        RealMatrix& operator=( const RealMatrix& source )
        {
#ifdef VERBOSE_REALMATRIX_CONSTRUCTION
            std::printf("...RealMatrix Copy assignment called.\n");
            std::fflush(stdout);
#endif
            // Check for self-assignment
            if ( this == &source )
                return *this;

            if ( _ptr )
                delete[] _ptr;

            _dim1 = source._dim1;
            _dim2 = source._dim2;
            _ptr = new double[ _dim1 * _dim2 ];
            std::copy( source._ptr, source._ptr + _dim1 * _dim2, _ptr );
            
            return *this;
        }

        // Move assignment operator
        RealMatrix& operator=( RealMatrix&& source )
        {
#ifdef VERBOSE_REALMATRIX_CONSTRUCTION
            std::printf("...RealMatrix Move assignment called.\n");
            std::fflush(stdout);
#endif
            int temp_dim1 = _dim1;
            int temp_dim2 = _dim2;
            double* temp_ptr = _ptr;

            _dim1 = source._dim1;
            _dim2 = source._dim2;
            _ptr = source.ptr();

            source._dim1 = temp_dim1;
            source._dim2 = temp_dim2;
            source._ptr = temp_ptr;
            
            return *this;
        }
        
        // Addition to self
        RealMatrix& operator+=( const RealMatrix& source )
        {
#ifndef NDEBUG
            if ( source._dim1 != _dim1 || source._dim2 != _dim2 )
                throw std::runtime_error("\nSize mismatch in operands for operator '+='!\n\tdim(A) = [ " + std::to_string(_dim1) + " x " 
                    + std::to_string(_dim2) + " ], dim(B) = [ " + std::to_string(source._dim1) + " x " + std::to_string(source._dim2) + " ]");
#endif
            cblas_daxpy(_dim1*_dim2, 1., source._ptr, 1, _ptr, 1);

            return *this;
        }
        
        RealMatrix& operator+=( RealMatrix&& source )
        {
#ifndef NDEBUG
            if ( source._dim1 != _dim1 || source._dim2 != _dim2 )
                throw std::runtime_error("\nSize mismatch in operands for operator '+='!\n\tdim(A) = [ " + std::to_string(_dim1) + " x " 
                    + std::to_string(_dim2) + " ], dim(B) = [ " + std::to_string(source._dim1) + " x " + std::to_string(source._dim2) + " ]");
#endif
            cblas_daxpy(_dim1*_dim2, 1., source._ptr, 1, _ptr, 1);

            return *this;
        }
        
        // Subtraction from self
        RealMatrix& operator-=( const RealMatrix& source )
        {
#ifndef NDEBUG
            if ( source._dim1 != _dim1 || source._dim2 != _dim2 )
                throw std::runtime_error("\nSize mismatch in operands for operator '+='!\n\tdim(A) = [ " + std::to_string(_dim1) + " x " 
                    + std::to_string(_dim2) + " ], dim(B) = [ " + std::to_string(source._dim1) + " x " + std::to_string(source._dim2) + " ]");
#endif
            cblas_daxpy(_dim1*_dim2, -1., source._ptr, 1, _ptr, 1);

            return *this;
        }
        
        RealMatrix& operator-=( RealMatrix&& source )
        {
#ifndef NDEBUG
            if ( source._dim1 != _dim1 || source._dim2 != _dim2 )
                throw std::runtime_error("\nSize mismatch in operands for operator '+='!\n\tdim(A) = [ " + std::to_string(_dim1) + " x " 
                    + std::to_string(_dim2) + " ], dim(B) = [ " + std::to_string(source._dim1) + " x " + std::to_string(source._dim2) + " ]");
#endif
            cblas_daxpy(_dim1*_dim2, -1., source._ptr, 1, _ptr, 1);

            return *this;
        }

        // In-place scalar multiplication
        RealMatrix& operator*=( double factor )
        {
            cblas_dscal( _dim1 * _dim2, factor, _ptr, 1 );
            return *this;
        }

        // In-place scalar division
        RealMatrix& operator/=( double factor )
        {
            cblas_dscal( _dim1 * _dim2, 1./factor, _ptr, 1 );
            return *this;
        }

        // NOTE: Storage of matrix components uses column-major format
        // Matrix component address access as variableName(i,j)
        double& operator()( int idx1, int idx2 )
        {
#ifndef NDEBUG
            if ( !_ptr )
                throw std::runtime_error("\nCannot access RealMatrix component (" + std::to_string(idx1) + "," + std::to_string(idx2) + ") -- "
                        + "\nmatrix is not initialized.");

            if ( idx1 < 0 || idx1 >= _dim1 || idx2 < 0 || idx2 >= _dim2 )
                throw std::runtime_error("\nCannot access RealMatrix component (" + std::to_string(idx1) + "," + std::to_string(idx2) + ")! "
                        + "Valid range is (0-" + std::to_string(_dim1 - 1) + ",0-" + std::to_string(_dim2 - 1) + ").");
#endif
            return _ptr[ idx2 * _dim1 + idx1 ];
        }

        // Matrix component value address access as variableName.elem(i,j)
        double operator()( int idx1, int idx2 ) const
        {
#ifndef NDEBUG
            if ( !_ptr )
                throw std::runtime_error("\nCannot access RealMatrix component (" + std::to_string(idx1) + "," + std::to_string(idx2) + ") -- "
                        + "\nmatrix is not initialized.");
            
            if ( idx1 < 0 || idx1 >= _dim1 || idx2 < 0 || idx2 >= _dim2 )
                throw std::runtime_error("\nCannot access RealMatrix component (" + std::to_string(idx1) + "," + std::to_string(idx2) + ")! "
                        + "Valid range is (0-" + std::to_string(_dim1 - 1) + ",0-" + std::to_string(_dim2 - 1) + ").");
#endif
            return _ptr[ idx2 * _dim1 + idx1 ];
        }

        // Dimensions of matrix
        int dim1() const { return _dim1; }
        int dim2() const { return _dim2; }

        // Erase contents
        void erase()
        {
            if ( _ptr )
                delete[] _ptr;
                
            _ptr = nullptr;
            _dim1 = 0;
            _dim2 = 0;
        }

        // Allocate/reallocate space for matrix of size (dim1 x dim2), with all
        // values set to zero
        void init( int dim1, int dim2 )
        {
            assert( dim1 >= 1 );
            assert( dim2 >= 1 );

            if ( _ptr)
                delete[] _ptr;

            _dim1 = dim1;
            _dim2 = dim2;
            _ptr = new double[ _dim1 * _dim2 ]();
        }

        // Show contents of matrix in scientific precision
        void print( const char *s, int n ) const 
        {
            std::printf( "\nRealMatrix %s:\n\n", s );

            if ( !_ptr )
                std::printf( "...is empty\n" );
            else
            {
                std::printf( "...size = %d x %d\n\n", _dim1, _dim2 );
                for ( int i = 0; i < _dim1; i++ )
                {
                    for ( int j = 0; j < _dim2; j++ )
                        std::printf( "%*.*e", n+10, n, _ptr[ j * _dim1 + i ] );
                    std::printf( "\n" );
                }
            }
            std::printf( "\n" );
        }

        // Give pointer to matrix components
        double* ptr() const { return _ptr; }
        
        // Transpose matrix
        RealMatrix trp() const
        {
            RealMatrix A( _dim2, _dim1 );
            double* t_ptr = A._ptr;

            for ( int i = 0; i < _dim1; i++ )
                for ( int j = 0; j < _dim2; j++ )
                    t_ptr[ i*_dim2 + j ] = _ptr[ j * _dim1 + i ];

            return A;
        }

    private:
        int     _dim1;
        int     _dim2;
        double* _ptr;
    };
}

#endif	/* REALMATRIX_HPP */
