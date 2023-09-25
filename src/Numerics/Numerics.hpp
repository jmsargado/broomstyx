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

#ifndef NUMERICS_HPP
#define	NUMERICS_HPP

#include <cstdio>
#include <map>
#include <vector>
#include <string>
#include <tuple>

#include "Core/BoundaryCondition.hpp"
#include "Core/Cell.hpp"
#include "Core/FieldCondition.hpp"
#include "Core/InitialCondition.hpp"
#include "Core/DofManager.hpp"
#include "Core/TimeData.hpp"

namespace broomstyx
{
    class AnalysisModel;
    class Material;
    
    class NumericsStatus
    {
    public:
        NumericsStatus();
        virtual ~NumericsStatus();
    };

    class Numerics
    {
        friend class NumericsManager;
        
    public:
        Numerics();
        virtual ~Numerics();
        
        // Disable copy constructor and assignment operator
        Numerics( const Numerics& ) = delete;
        Numerics& operator=( const Numerics& ) = delete;

        std::string giveName();
        
        std::tuple<RealVector,RealVector> giveCellFieldOutputAtEvaluationPointsOf( Cell* targetCell, int fieldNum );
        int  giveSpatialDimension() const;
        void readDataFrom( FILE* fp );
        int  requiredNumberOfDofsPerCell() const;
        int  requiredNumberOfMaterials() const;
        int  requiredNumberOfDofsPerNode() const;
        int  requiredNumberOfNodes() const;
        void setIdTo( int id );
        void setStageTo( int stage );

        virtual void initializeMaterialsAt( Cell* targetCell );
        virtual bool performAdditionalConvergenceCheckAt( Cell* targetCell );
        virtual void performPreIterationOperationsAt( int iterNum );
        virtual void printPostIterationMessage();
        virtual void readAdditionalDataFrom( FILE* fp );
        virtual void removeConstraintsOn( Cell* targetCell );

        virtual void finalizeDataAt( Cell* targetCell, const TimeData& time ) = 0;
        virtual void deleteNumericsAt( Cell* targetCell ) = 0;
        virtual void initializeNumericsAt( Cell* targetCell ) = 0;
        
        virtual std::tuple< std::vector<Dof*>
                          , std::vector<Dof*>
                          , RealVector >
            giveStaticCoefficientMatrixAt( Cell*           targetCell
                                         , int             subsys
                                         , const TimeData& time );

        virtual std::tuple< std::vector<Dof*>, RealVector >
            giveStaticLeftHandSideAt( Cell*           targetCell
                                    , int             subsys
                                    , const TimeData& time );

        virtual std::tuple< std::vector<Dof*>, RealVector >
            giveStaticRightHandSideAt( Cell*                    targetCell
                                     , int                      subsys
                                     , const BoundaryCondition& bndCond
                                     , const TimeData&          time );

        virtual std::tuple< std::vector<Dof*>, RealVector >
            giveStaticRightHandSideAt( Cell*                 targetCell
                                     , int                   subsys
                                     , const FieldCondition& fldCond
                                     , const TimeData&       time );

        virtual std::tuple< std::vector<Dof*>
                          , std::vector<Dof*>
                          , RealVector >
            giveTransientCoefficientMatrixAt( Cell*           targetCell
                                            , int             subsys
                                            , const TimeData& time );

        virtual std::tuple< std::vector<Dof*>, RealVector >
            giveTransientLeftHandSideAt( Cell*           targetCell
                                       , int             subsys
                                       , const TimeData& time
                                       , ValueType       valType );
        
        virtual void
            imposeConstraintAt( Cell*                    targetCell
                              , const BoundaryCondition& bndCond
                              , const TimeData&          time );
        
        // Error-generating virtual functions (must be implemented in derived class when called)
        virtual double 
            giveCellFieldValueAt( Cell* targetCell, int fieldNum );
        
        virtual RealVector 
            giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum );
        
        virtual std::vector<RealVector> 
            giveEvaluationPointsFor( Cell* targetCell );
        
        virtual std::tuple< RealVector,RealVector > 
            giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag );
        
        virtual RealVector
            giveNumericsParameter( const std::string& paramTag );
        
        virtual void
            imposeInitialConditionAt( Cell*                   targetCell
                                    , const InitialCondition& initCond );

        virtual void performPostprocessingAt( Cell* targetCell, std::string tag );
        virtual void performPrefinalizationCalculationsAt( Cell* targetCell );
        virtual void performPreprocessingAt( Cell* targetCell, std::string tag );
        virtual void setDofStagesAt( Cell* targetCell );
        
    protected:
        std::string _name;

        int  _id;
        int  _stage;
        bool _stageAssigned;

        int _dim;
        int _nDofsPerCell;
        int _nDofsPerNode;
        int _nMaterials;
        int _nNodes;
        int _nSubsystems;

        // Mapping from cell field number to field argument
        int _nCellFieldOutput;
        std::map<int,std::string> _cellFieldOutput;

        std::vector<int> _nodalDof;
        std::vector<int> _cellDof;
        std::vector<int> _subsystem;
        
        // Helper methods
        std::vector<Material*> giveMaterialSetFor( Cell* targetCell ) const;
        void error_unimplemented( const std::string& method );
    };
}

#endif	/* NUMERICS_HPP */
