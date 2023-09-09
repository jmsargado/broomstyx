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
        
        void                   createPhysicalEntity( int dim, int number, std::string label );
        PhysicalEntity         givePhysicalEntity( int n );
        std::vector<Material*> giveMaterialSetForDomain( int label, int stage );
        int                    giveNumberOfPhysicalNames();
        Numerics*              giveNumericsForDomain( int label, int stage );
        std::string            givePhysicalEntityNameFor( int physEntNum );
        int                    givePhysicalEntityNumberFor( std::string name );
        void                   readDomainAssignmentsFrom( FILE* fp );
        void                   readNumberOfStagesFrom( FILE* fp );

        // Methods involving node access
        
        void            countNodes();
        std::set<Cell*> giveCellsAttachedTo( Node* node, int dim );
        RealVector      giveCoordinatesOf( Node* node );
        Dof*   giveNodalDof( int dofNum, Node* node );
        double giveFieldValueAt( Node* node, int fieldNum );
        int    giveIdOf( Node* node );
        Node*  giveNode( int nodeNum );
        int    giveNumberOfNodes();
        void   makeNewNodeAt( RealVector& location );
        void   performNodalPostProcessing();
        void   readNumberOfFieldsPerNodeFrom( FILE* fp );
        void   setCoordinatesOf( Node* targetNode, const RealVector& coor );
        void   setFieldValueAt( Node* targetNode, int fieldNum, double val );
        
        // Methods involving cell access

        void  countCells();
        void  finalizeCellDataAt( const TimeData& time, int stage );
        void  findCellAttachments();
        // void   findBoundaryAssociations();
        // void   findDomainCellsAssociatedWith( Cell* targetCell );
        void  findCellsAttachedTo( Cell *targetCell );
        // void   formDomainPartitions();
        // Cell*  giveBoundaryCell( int cellNum );
        Dof*  giveCellDof( int dofNum, Cell *targetCell );
        // Cell*  giveDomainCell( int cellNum );
        // std::vector<Cell*> 
        //        giveDomainCellsAssociatedWith( Cell *targetCell );
        Cell* giveDomainCellInPartition( int partNum, int cellNum );
        int   giveElementTypeOf( Cell* targetCell );
        int   giveIdOf( Cell *targetCell );
        int   giveLabelOf( Cell *targetCell );
        
        std::set<Cell*> giveNeighborsOf( Cell* targetCell, int dim );
        std::vector<Node*> giveNodesOf( Cell *targetCell );
        
        // int   giveNumberOfBoundaryCells();
        // int   giveNumberOfDomainCells();
        // int   giveNumberOfDomainCellsInPartition( int partNum );
        int   giveNumberOfNodesOf( Cell *targetCell );
        // int   giveNumberOfPartitions();
        Numerics* giveNumericsFor( Cell* targetCell, int stage );
        void  initializeMaterialsAtCells();
        void  initializeNumericsAtCells();
        Cell* makeNewCellWithLabel( int cellLabel, int dim );
        // void  makeNewFaceBetween( Cell* posCell, Cell* negCell, int posFaceNum );
        void  mustConstructFaces();
        void  readNumberOfFieldsPerCellFrom( FILE* fp );
        // void  readNumberOfFieldsPerFaceFrom( FILE* fp );
        // void  removeAllCellConstraints();
        void  reorderNodesOf( Cell* targetCell, std::vector<int>& reordering );
        void  reportDetailedStatus();
        void  reportStatus();
        void  setElementTypeOf( Cell* targetCell, int elemType );
        void  setNeighborsOf( Cell *targetCell, std::set<Cell*>& neighbors, int dim );
        void  setNodesOf( Cell *targetCell, std::vector<int>& cellNodes );

    private:
        std::string _name;

        std::vector<PhysicalEntity> _physEnt;
        
        int _nStages;
        std::vector<std::map<std::string, Numerics*> > _numerics;
        std::vector<std::map<std::string, std::vector< Material*> > > _materialSet;
        
        int _fieldsPerNode;
        std::list<Node*>   _nodeList;
        std::vector<Node*> _node;
                
        int _fieldsPerCell;

        // Separate lists and vectors for cells of different dimensions
        std::list<Cell*>   _cellList[4];
        std::vector<Cell*> _cell[4];

        // std::vector<std::vector<Cell*> > _domCell;

//         std::list<Cell*> _bndCellList;
//         std::vector<Cell*> _bndCell;
//
//         std::list<Cell*> _domCellList;
//         std::vector<Cell*> _domCell;
        
        // Cells representing faces of domain cells
        // int _fieldsPerFace;
        // bool _constructFaces;
        // std::list<Cell*> _faceList;
        // std::vector<Cell*> _face;
        
        std::vector<std::vector<Cell*> > _partition;
        
        DomainManager();
        virtual ~DomainManager();
    };
}

#endif	/* DOMAINMANAGER_HPP */
