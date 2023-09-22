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

#ifndef REALVECTOR_HPP
#define	REALVECTOR_HPP

#include <config.h>

#include <cstdio>
#include <stdexcept>
#include <initializer_list>
#include "RealMatrix.hpp"

#ifdef HAVE_MKL
    #include "mkl_cblas.h"
#else
    #include "cblas.h"
#endif

namespace broomstyx
{
    class RealVector final
    {
        friend RealVector operator*( RealVector&& A, double b );
        friend RealVector operator*( double a, RealVector&& B );
        friend RealVector operator/( RealVector&& A, double b );
        
    public:
        // Default constructor
        RealVector()
            : _dim( 0 )
            , _ptr( nullptr )
        {}

        // Constructor with initial memory allocation
        RealVector( int dim )
        {
#ifndef NDEBUG
            if ( dim < 1 )
                throw std::runtime_error("\nCannot construct RealVector with dim = " + std::to_string(dim));
#endif
            _dim = dim;
            _ptr = new double[ _dim ]();
        }
        
        // Constructor with initializer list
        RealVector( std::initializer_list<double> initList )
        {
            _dim = (int)initList.size();
            _ptr = new double[ _dim ]();

            std::copy( initList.begin(), initList.end(), _ptr );
        }

        // Copy constructor
        RealVector( const RealVector& source )
        {
            if ( source._dim > 0 )
            {
                _dim = source._dim;
#ifdef VERBOSE_REALVECTOR_CONSTRUCTION
                std::printf( "...RealVector Copy constructor called.\n" );
                std::fflush( stdout );
#endif
                _ptr = new double[ _dim ];
                std::copy( source._ptr, source._ptr + _dim, _ptr );
            }
            else
            {
                _dim = 0;
                _ptr = nullptr;
            }
        }

        // Move constructor
        RealVector( RealVector&& source )
            : _dim( source._dim )
            , _ptr( source._ptr )
        {
#ifdef VERBOSE_REALVECTOR_CONSTRUCTION
            std::printf( "...RealVector Move constructor called.\n" );
            std::fflush( stdout );
#endif
            source._ptr = nullptr;
            source._dim = 0;
        }

        // Destructor
        virtual ~RealVector()
        {
            if ( _ptr )
                delete[] _ptr;
        }

        // Assignment with initializer list
        RealVector& operator= ( std::initializer_list<double> initList )
        {
            if ( _ptr && _dim != (int)initList.size() )
            {
                delete[] _ptr;
                _ptr = nullptr;
            }
            
            _dim = (int)initList.size();
            if ( !_ptr )
                _ptr = new double[ _dim ];

            std::copy( initList.begin(), initList.end(), _ptr );

            return *this;
        }
        
        // Copy assignment operator
        RealVector& operator= ( const RealVector& source )
        {
#ifdef VERBOSE_REALVECTOR_CONSTRUCTION
            std::printf("...RealVector Copy assignment called.\n");
            std::fflush(stdout);
#endif
            // Check for self-assignment
            if ( this == &source )
                return *this;
            
            if ( _ptr )
                delete[] _ptr;

            _dim = source._dim;
            _ptr = new double[ _dim ];
            std::copy( source._ptr, source._ptr + _dim, _ptr );
            
            return *this;
        }

        // Move assignment operator
        RealVector& operator= ( RealVector&& source )
        {
#ifdef VERBOSE_REALVECTOR_CONSTRUCTION
            std::printf("...RealVector Move assignment called.\n");
            std::fflush(stdout);
#endif
            if ( _ptr )
                delete[] _ptr;
                
            _dim = source._dim;
            _ptr = source._ptr;

            source._dim = 0;
            source._ptr = nullptr;

            return *this;
        }

        // Addition to self
        RealVector& operator+= ( const RealVector& source )
        {
#ifndef NDEBUG
            if ( source._dim != _dim )
                throw std::runtime_error( "\nSize mismatch in operands for operator '+='!\n\tdim(A) = "
                    + std::to_string( _dim ) + ", dim(B) = " + std::to_string( source._dim ) );
#endif
            cblas_daxpy( _dim, 1.0, source._ptr, 1, _ptr, 1 );

            return *this;
        }
        
        RealVector& operator+= ( RealVector&& source )
        {
#ifndef NDEBUG
            if ( source._dim != _dim )
                throw std::runtime_error( "\nSize mismatch in operands for operator '+='!\n\tdim(A) = "
                    + std::to_string( _dim ) + ", dim(B) = " + std::to_string( source._dim ) );
#endif
            cblas_daxpy( _dim, 1.0, source._ptr, 1, _ptr, 1 );

            return *this;
        }
        
        // Subtraction from self
        RealVector& operator-= ( const RealVector& source )
        {
#ifndef NDEBUG
            if ( source._dim != _dim )
                throw std::runtime_error( "\nSize mismatch in operands for operator '-='!\n\tdim(A) = "
                    + std::to_string( _dim ) + ", dim(B) = " + std::to_string( source._dim ) );
#endif
            cblas_daxpy( _dim, -1.0, source._ptr, 1, _ptr, 1 );

            return *this;
        }
        
        RealVector& operator-= ( RealVector&& source )
        {
#ifndef NDEBUG
            if ( source._dim != _dim )
                throw std::runtime_error( "\nSize mismatch in operands for operator '-='!\n\tdim(A) = "
                    + std::to_string( _dim ) + ", dim(B) = " + std::to_string( source._dim ) );
#endif
            cblas_daxpy( _dim, -1.0, source._ptr, 1, _ptr, 1 );

            return *this;
        }

        // In-place scalar multiplication
        RealVector& operator*= ( double factor )
        {
            cblas_dscal( _dim, factor, _ptr, 1 );
            return *this;         
        }

        // In-place scalar division
        RealVector& operator/= ( double factor )
        {
            cblas_dscal( _dim, 1./factor, _ptr, 1 );
            return *this;         
        }
        
        // Vector component access as variableName(i)
        double& operator()( int idx )
        {
#ifndef NDEBUG
            if ( !_ptr )
                throw std::runtime_error( "\nCannot access RealVector component ("
                    + std::to_string( idx ) + ") -- vector is not initialized." );

            if ( idx < 0 || idx >= _dim )
                throw std::runtime_error( "\nCannot access RealVector component (" + std::to_string( idx )
                    + ")! Valid range is (0-" + std::to_string( _dim - 1 ) + ")." );
#endif
            return _ptr[ idx ];
        }

        double operator()( int idx ) const
        {
#ifndef NDEBUG
            if ( !_ptr )
                throw std::runtime_error( "\nCannot access RealVector component (" + std::to_string( idx )
                      + ") -- vector is not initialized." );

            if ( idx < 0 || idx >= _dim )
                throw std::runtime_error( "\nCannot access RealVector component (" + std::to_string( idx )
                      + ")! Valid range is (0-" + std::to_string( _dim - 1 ) + ")." );
#endif
            return _ptr[ idx ];
        }

        // Vector cross product
        RealVector cross( const RealVector& B )
        {
#ifndef NDEBUG
            if ( _dim != 3 || B._dim != 3 )
                throw std::runtime_error( "\nVector cross product only operates on vectors with dim = 3!" );
#endif
            double c0 = _ptr[ 1 ]*B._ptr[ 2 ] - _ptr[ 2 ]*B._ptr[ 1 ];
            double c1 = _ptr[ 2 ]*B._ptr[ 0 ] - _ptr[ 0 ]*B._ptr[ 2 ];
            double c2 = _ptr[ 0 ]*B._ptr[ 1 ] - _ptr[ 1 ]*B._ptr[ 0 ];

            return RealVector( { c0, c1, c2 } );
        }

        // Dimensions of vector
        int dim() const { return _dim; }

        // Vector dot product
        double dot( const RealVector& B )
        {
#ifndef NDEBUG
            if ( B._dim != _dim )
                throw std::runtime_error( "\nSize mismatch in vector dot product!\ndim(A) = " + std::to_string( _dim )
                      + ", dim(B) = " + std::to_string( B._dim ) );
#endif
            return cblas_ddot( _dim, _ptr, 1, B._ptr, 1 );
        }

        double dot( RealVector&& B )
        {
#ifndef NDEBUG
            if ( B._dim != _dim )
                throw std::runtime_error( "\nSize mismatch in vector dot product!\ndim(A) = " + std::to_string( _dim )
                      + ", dim(B) = " + std::to_string( B._dim ) );
#endif
            return cblas_ddot( _dim, _ptr, 1, B._ptr, 1 );
        }

        // Erase contents
        void erase()
        {
            if ( _ptr )
                delete[] _ptr;
                
            _ptr = nullptr;
            _dim = 0;
        }

        // Allocate/reallocate space for vector, initializing all values to zero
        void init( int dim )
        {
#ifndef NDEBUG
            if ( dim < 1 )
                throw std::runtime_error( "\nCannot initialize RealVector with dim = " + std::to_string( dim ) );
#endif
            if ( _ptr )
                delete[] _ptr;

            _dim = dim;
            _ptr = new double[ _dim ]();
        }
        
        // Show contents of matrix in scientific precision
        void print( const char *s, int n ) const
        {
            std::printf( "\nRealVector %s:\n\n", s );
            
            if ( !_ptr )
                std::printf( "...is empty\n" );
            else
            {
                std::printf( "\n" );
                for ( int i = 0; i < _dim; i++ )
                    std::printf( "%*.*e\n", n + 10, n, _ptr[ i ] );
            }
            std::printf( "\n" );
        }

        // Print contents to file
        void printTo( FILE* fp, int n ) const
        {
            if ( !_ptr )
                std::fprintf( fp, "...is empty\n" );
            else
                for ( int i = 0; i < _dim; i++ )
                    std::fprintf( fp, "%*.*e\n", n + 10, n, _ptr[ i ] );
        }
        
        // Give pointer to vector components
        double* ptr() const { return _ptr; }
        
        // Vector tensor product
        RealMatrix xMen( const RealVector& B ) const
        {
            int dimB = B.dim();
            RealMatrix C( _dim, dimB );
            
            cblas_dger( CblasColMajor, _dim, dimB, 1.0, _ptr, 1, B._ptr, 1, C.ptr(), _dim );

            return C;
        }

        RealMatrix xMen( RealVector&& B ) const
        {
            int dimB = B.dim();
            RealMatrix C( _dim, dimB );
            
            cblas_dger( CblasColMajor, _dim, dimB, 1.0, _ptr, 1, B._ptr, 1, C.ptr(), _dim );

            return C;
        }
        
    private:
        int     _dim;
        double* _ptr;
    };
}

#endif	/* REALVECTOR_HPP */
