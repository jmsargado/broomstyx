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

#include "Cell.hpp"

#include "AnalysisModel.hpp"
#include "DomainManager.hpp"
#include "Dof.hpp"
#include "Node.hpp"

using namespace broomstyx;

Cell::Cell( int id, int elType, int label )
    : _elType( elType )
    , _label( label )
    , _id( id )
    , _partition( 0 )
{}

Cell::~Cell() = default;

// Public methods
// ---------------------------------------------------------------------------------------
int Cell::dimension() const
{
    return analysisModel().domainManager().giveDimensionOfPhysicalEntity( _label );
}
// ---------------------------------------------------------------------------------------
int Cell::id() const
{
    return _id;
}
// ---------------------------------------------------------------------------------------
int Cell::label() const
{
    return _label;
}
// ---------------------------------------------------------------------------------------
void Cell::showInfo()
{
    int cellDim = analysisModel().domainManager().giveDimensionOfPhysicalEntity( _label );
    std::printf( "\n  Cell ID = %d, label = %d, dim = %d\n", _id, _label, cellDim );

    std::printf( "  Cell nodes: " );
    for ( auto& i : _node )
        std::printf( "%d ", i->id() );
    
    std::printf( "\n  Attached Cells\n" );
    for ( int dim : { 0, 1, 2, 3 } )
    {
        std::printf( "   dim = %d: ", dim );
        for ( auto& it : _attachedCell[ dim ] )
            std::printf( "%d ", it->id() );
        std::printf( "\n" );
    }

    std::printf( "\n  Neighbors: " );
    for ( auto& i : _neighbor )
        if ( i )
            std::printf( "%d ", i->id() );
        else
            std::printf( "none" );
    std::printf( "\n" );
}