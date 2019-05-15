/*
  Copyright (c) 2014 - 2019 University of Bergen

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

#ifndef ISTL_CPP
#define ISTL_CPP

#include <config.h>

#ifdef HAVE_DUNE_ISTL

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <string>
#include "../Core/ObjectFactory.hpp"
#include "../SparseMatrix/ISTLMat.hpp"
#include "../Util/readOperations.hpp"
#include "../Util/RealVector.hpp"

#include "ISTL.hpp"

using namespace broomstyx;

registerBroomstyxObject(LinearSolver, ISTL)

namespace detail
{

    template<class M, class X, class Y>
    class OpenMPMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
    {
    public:
        //! export types
        typedef M matrix_type;
        typedef X domain_type;
        typedef Y range_type;
        typedef typename X::field_type field_type;

        //! constructor: just store a reference to a matrix
        explicit OpenMPMatrixAdapter (const M& A) : _A_(A) {}

        //! apply operator to x:  \f$ y = A(x) \f$
        virtual void apply (const X& x, Y& y) const
        {
#ifdef _OPENMP
            const size_t nRows = _A_.N();

#pragma omp parallel for
            for( size_t i = 0; i < nRows; ++i )
            {
                const auto& row = _A_[ i ];

                y[ i ] = 0;
                const auto endj = row.end();
                for (auto j=row.begin(); j!=endj; ++j)
                {
                    (*j).umv( x[ j.index() ], y[ i ]);
                }
            }
#else
            _A_.mv(x,y);
#endif
        }

        //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
        virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
        {
#ifdef _OPENMP
            const size_t nRows = _A_.N();

#pragma omp parallel for
            for( size_t i = 0; i<nRows; ++i )
            {
              const auto& row = _A_[ i ];

              const auto endj = row.end();
              for (auto j=row.begin(); j!=endj; ++j)
              {
                (*j).usmv(alpha, x[ j.index() ], y[ i ]);
              }
            }
#else
            _A_.usmv(alpha,x,y);
#endif
        }

        //! get matrix via *
        virtual const M& getmat () const
        {
            return _A_;
        }

        //! Category of the solver (see SolverCategory::Category)
        virtual Dune::SolverCategory::Category category() const
        {
            return Dune::SolverCategory::sequential;
        }

    private:
        const M& _A_;
    };
}


ISTL::ISTL() {}

ISTL::~ISTL() {}

std::string ISTL::giveRequiredMatrixFormat()
{
    return std::string("ISTLMat");
}

bool ISTL::giveSymmetryOption()
{
    return false;
}

void ISTL::readDataFrom( FILE* fp )
{
	std::string src = "ISTL (LinearSolver)";

    _tol = getRealInputFrom(fp, "Failed to read relative tolerance for iterative linear solver from input file!", src);
    _maxIter = getIntegerInputFrom(fp, "Failed to read max. iterations for iterative linear solver from input file!", src);
    _relax = getRealInputFrom(fp, "Failed to read ILU0 preconditioner relaxation parameter from input file!", src);
}

void ISTL::setInitialGuessTo( RealVector& initGuess )
{
    _initGuess = initGuess;
}

RealVector ISTL::solve( SparseMatrix* coefMat, RealVector& rhs )
{
    ISTLMat* matrix = dynamic_cast< ISTLMat* > (coefMat);
    if( ! matrix )
    {
        std::cerr << "ERROR" << std::endl;
        std::abort();
    }

    typedef typename ISTLMat :: MatrixType M;

    M& mat = matrix->exportMatrix();

    typedef detail::OpenMPMatrixAdapter< M, BlockVectorType, BlockVectorType > AssembledOperatorType;
    AssembledOperatorType op( mat );
    
    typedef Dune::SeqILU< M, BlockVectorType, BlockVectorType > PreconditionerType;
    PreconditionerType precon( mat, _relax, true );

    Dune::BiCGSTABSolver< BlockVectorType > solver( op, precon, _tol, int(_maxIter), int(0));

    const int dim = rhs.dim();
    BlockVectorType B( rhs.dim() );
    BlockVectorType X( rhs.dim() );
    for( int i=0; i<dim; ++i )
    {
        B[ i ][0] = rhs( i );
        X[ i ][0] = _initGuess( i );
    }

    Dune::InverseOperatorResult info;
    solver.apply( X, B, _tol, info );

    RealVector x( dim );
    for( int i=0; i<dim; ++i )
    {
        x( i ) = X[ i ][ 0 ];
    }

    return x;
}

#endif /* HAVE_DUNE_ISTL */
#endif /* ISTL_CPP */