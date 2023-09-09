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

#ifndef CELL_HPP
#define CELL_HPP

#include <cstdio>
#include <cstdlib>
#include <set>
#include <vector>
#include "Math/RealVector.hpp"

namespace broomstyx
{    
    class Dof;
    class Node;
    class NumericsStatus;

    class Cell
    {
        friend class DomainManager;
        friend class DofManager;

    public:
        // Constructor with arguments
        Cell( int elType, int label, int dim );
        
        // Destructor
        virtual ~Cell();
        
        // Disable default constructor, copy constructor and assignment operator
        Cell() = delete;
        Cell( const Cell& ) = delete;
        Cell& operator=( const Cell& ) = delete;

        RealVector cellData;
        std::vector<NumericsStatus*> numericsStatus;
        
        int  id();
        void showInfo();

    private:
        int _elType;
        int _label;
        int _dim;
        int _id;
        int _partition;

        std::vector<bool>  _hasNumericsAtStage;
        std::vector<Node*> _node;
        std::vector<Dof*>  _dof;
        std::set<Cell*>    _attachedCell[4];
        std::vector<Cell*> _neighbor;
    };
}

#endif /* CELL_HPP */