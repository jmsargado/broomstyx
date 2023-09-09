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

#ifndef MESHREADER_HPP
#define	MESHREADER_HPP

#include <cstdio>
#include <string>
#include <vector>

namespace broomstyx
{
    class AnalysisModel;

    class MeshReader
    {
    public:
        MeshReader() {}
        virtual ~MeshReader() {}
        
        // Disable copy constructor and assignment operator
        MeshReader( const MeshReader& ) = delete;
        MeshReader& operator=( const MeshReader& ) = delete;

        virtual void readMeshFile( std::string filename ) = 0;
        virtual std::vector<int> giveFaceNodeNumbersForElementType( int elType, int face ) = 0;
        virtual int  giveNumberOfFacesForElementType( int elType ) = 0;
    };
}

#endif	/* MESHREADER_HPP */