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

#ifndef DOMAINMANAGER_HPP
#define	DOMAINMANAGER_HPP

#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <list>

#include "TimeData.hpp"
#include "Math/RealVector.hpp"

namespace broomstyx
{
    class Cell;
    class Dof;
    class Material;
    class Node;
    class Numerics;

    class DomainManager final
    {
        friend class AnalysisModel;
        
    public:
        // Disable copy constructor and assignment operator
        DomainManager( const DomainManager& ) = delete;
        DomainManager& operator=( const DomainManager& ) = delete;
        
        struct PhysicalEntity
        {
            int dimension;
            int entityNumber;
            std::string name;
        };
        
        void                   createPhysicalEntity( int dim, int number, const std::string& label );
        int                    giveDimensionOfPhysicalEntity( int physEntNum ) const;
        PhysicalEntity         givePhysicalEntity( int n ) const;
        std::vector<Material*> giveMaterialSetForDomain( int label, int stage ) const;
        int                    giveNumberOfPhysicalNames();
        Numerics*              giveNumericsForDomain( int label, int stage ) const;
        std::string            giveNameOfPhysicalEntity( int physEntNum ) const;
        int                    givePhysicalEntityNumberFor( const std::string& name ) const;
        void                   readDomainAssignmentsFrom( FILE* fp );

        // Methods involving node access
        
        void                    countNodes();
        static std::set<Cell*>& giveCellsAttachedTo( Node* node, int dim );
        static RealVector       giveCoordinatesOf( Node* node );
        static Dof*             giveNodalDof( int dofNum, Node* node );
        static double           giveFieldValueAt( Node* node, int fieldNum );
        static int              giveIdOf( Node* node );
        Node*  giveNode( int nodeNum ) const;
        int    giveNumberOfNodes() const;
        void   makeNewNodeAt( RealVector& location );
        void   performNodalPostProcessing() const;
        void   readNumberOfFieldsPerNodeFrom( FILE* fp );
        static void   setCoordinatesOf( Node* targetNode, const RealVector& coor );
        static void   setFieldValueAt( Node* targetNode, int fieldNum, double val );
        
        // Methods involving cell access

        void  countCells();
        void  finalizeCellDataAt( const TimeData& time, int stage );
        void  findCellAttachments();
        void  findCellsAttachedTo( Cell* targetCell );
        void  formDomainPartitions();
        Cell* giveCell( int num, int dim );
        static Dof*  giveCellDof( int dofNum, Cell* targetCell );
//        Cell* giveDomainCellInPartition( int partNum, int cellNum );
        static int   giveElementTypeOf( Cell* targetCell );
        static int   giveIdOf( Cell *targetCell );
        static std::vector<Cell*> giveNeighborsOf( Cell* targetCell );
        static std::vector<Node*> giveNodesOf( Cell* targetCell );
        int   giveNumberOfCellsWithDimension( int dim );
        static int   giveNumberOfNodesOf( Cell* targetCell );
        Numerics* giveNumericsFor( Cell* targetCell, int stage ) const;
        void  initializeMaterialsAtCells();
        void  initializeNumericsAtCells();
        Cell* makeNewCell( int id, int elType, int label );
//        void  mustConstructFaces();
        void  readNumberOfFieldsPerCellFrom( FILE* fp );
        void  removeAllCellConstraints();
        void  reorderNodesOf( Cell* targetCell, std::vector<int>& reordering );
        void  reportDetailedStatus();
        void  reportStatus();
        static void setElementTypeOf( Cell* targetCell, int elemType );
        static void setHaloOf( Cell *targetCell, std::vector<int>& halo );
        void  setNodesOf( Cell* targetCell, std::vector<int>& cellNodes );
        static void setPartitionOf( Cell* targetCell, int partition );

    private:
        std::string _name;
        int _nStage;

        std::vector<PhysicalEntity> _physEnt;
        
        std::vector<std::map<int, Numerics*> > _numerics;
        std::vector<std::map<int, std::vector< Material*> > > _materialSet;
        
        int _fieldsPerNode;
        std::list<Node*>   _nodeList;
        std::vector<Node*> _node;
                
        int _fieldsPerCell;

        // Separate lists and vectors for cells of different dimensions
        std::list<Cell*>   _cellList[4];
        std::vector<Cell*> _cell[4];
        
        std::vector<std::vector<Cell*> > _partition;
        
        DomainManager();
        virtual ~DomainManager();
    };
}

#endif	/* DOMAINMANAGER_HPP */
