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

#ifndef OUTPUTWRITER_HPP
#define	OUTPUTWRITER_HPP

#include <cstdio>
#include <string>

namespace broomstyx
{
    class AnalysisModel;
    
    class OutputWriter
    {
    public:
        OutputWriter() {}
        virtual ~OutputWriter() {}
        
        // Disable copy constructor and assignment operator
        OutputWriter( const OutputWriter& ) = delete;
        OutputWriter& operator=( const OutputWriter& ) = delete;

        // Virtual public methods
        virtual void initialize() = 0;
        virtual void readDataFrom( FILE* fp ) = 0;
        virtual void writeOutput( double time ) = 0;

    protected:
        std::string _name;
    };
}

#endif	/* OUTPUTWRITER_HPP */

