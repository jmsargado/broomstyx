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

#ifndef LINF_LINF_STOPPINGCRITERION_HPP
#define	LINF_LINF_STOPPINGCRITERION_HPP

#include "ConvergenceChecker.hpp"
#include "Util/RealMatrix.hpp"

namespace broomstyx
{
    class LInf_LInf_StoppingCriterion : public ConvergenceChecker
    {
    public:
        LInf_LInf_StoppingCriterion();
        virtual ~LInf_LInf_StoppingCriterion();

        bool checkConvergenceOf( const std::vector<RealVector>& resid
                               , const std::vector<int>& subsysNumbers
                               , const std::vector<Dof*>& dof );
        void initialize( int nDofGroups );
        void processLocalResidualContribution( RealVector& contrib, std::vector<int>& dofGrp, int threadNum );
        void reportConvergenceStatus();
        void resetResidualCriteria();

    private:
        int _nThreads;
        RealVector _contribCount;
        RealVector _dofGrpCount;
        
        RealMatrix _contribCountPerThread;
        RealMatrix _corrCritPerThread;
        RealMatrix _corrNormPerThread;
        RealMatrix _dofGrpCountPerThread;
        RealMatrix _residCritPerThread;
        RealMatrix _residNormPerThread;
        
        int giveIndexForDofGroup( int dofGroupNum )
        {
            int idx = -1;
            for ( int i = 0; i < _nDofGroups; i++ )
                if ( dofGroupNum == _dofGrpNum[i] )
                    idx = i;
            
            return idx;
        }

        int giveIndexForSubsystem( int subsysNum, const std::vector<int>& subsysNumVector )
        {
            int idx = -1;
            for ( int i = 0; i < (int)subsysNumVector.size(); i++ )
                if ( subsysNum == subsysNumVector[i] )
                    idx = i;
            
            return idx;
        }
    };
}

#endif	/* LINF_LINF_STOPPINGCRITERION_HPP */