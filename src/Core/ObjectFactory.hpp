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

#ifndef OBJECTFACTORY_HPP
#define OBJECTFACTORY_HPP

#include <map>
#include <string>

#define registerBroomstyxObject( baseClass, derivedClass ) \
    static bool dummy_ ## derivedClass __attribute__((unused)) = \
    objectFactory().register ## baseClass( # derivedClass, createDerivedObject<baseClass,derivedClass> );

namespace broomstyx
{
    class ConvergenceCriterion;
    class LinearSolver;
    class Material;
    class MeshReader;
    class Numerics;
    class OutputQuantity;
    class OutputWriter;
    class SolutionMethod;
    class SparseMatrix;
    class UserFunction;
    
    template<typename Base, typename Derived> Base* createDerivedObject()
    {
        return new Derived();
    }

    class ObjectFactory
    {
    public:
        ObjectFactory();
        virtual ~ObjectFactory();
        
        // Disable copy constructor and assignment operator
        ObjectFactory( const ObjectFactory& ) = delete;
        ObjectFactory& operator=( const ObjectFactory& ) = delete;
        
        bool hasError() const;
        
        ConvergenceCriterion* instantiateConvergenceCriterion( const std::string& name );
        LinearSolver*         instantiateLinearSolver( const std::string& name );
        Material*             instantiateMaterial( const std::string& name );
        MeshReader*           instantiateMeshReader( const std::string& name );
        Numerics*             instantiateNumerics( const std::string& name );
        OutputQuantity*       instantiateOutputQuantity( const std::string& name );
        OutputWriter*         instantiateOutputWriter( const std::string& name );
        SolutionMethod*       instantiateSolutionMethod( const std::string& name );
        SparseMatrix*         instantiateSparseMatrix( const std::string& name );
        UserFunction*         instantiateUserFunction( const std::string& name );
        
        bool registerConvergenceCriterion( const std::string& name, ConvergenceCriterion* (*)() );
        bool registerLinearSolver( const std::string& name, LinearSolver* (*)() );
        bool registerMaterial( const std::string& name, Material* (*)() );
        bool registerMeshReader( const std::string& name, MeshReader* (*)() );
        bool registerNumerics( const std::string& name, Numerics* (*)() );
        bool registerOutputQuantity( const std::string& name, OutputQuantity* (*)() );
        bool registerOutputWriter( const std::string& name, OutputWriter* (*)() );
        bool registerSolutionMethod( const std::string& name, SolutionMethod* (*)() );
        bool registerSparseMatrix( const std::string& name, SparseMatrix* (*)() );
        bool registerUserFunction( const std::string& name, UserFunction* (*)() );

    private:
        std::map< std::string, ConvergenceCriterion* (*)() > _convergenceCriterion;
        std::map< std::string, LinearSolver* (*)() >         _linearSolver;
        std::map< std::string, Material* (*)() >             _material;
        std::map< std::string, MeshReader* (*)() >           _meshReader;
        std::map< std::string, Numerics* (*)() >             _numerics;
        std::map< std::string, OutputQuantity* (*)() >       _outputQuantity;
        std::map< std::string, OutputWriter* (*)() >         _outputWriter;
        std::map< std::string, SolutionMethod* (*)() >       _solutionMethod;
        std::map< std::string, SparseMatrix* (*)() >         _sparseMatrix;
        std::map< std::string, UserFunction* (*)() >         _userFunction;
        
        bool _errorInRegistration;
    };

    ObjectFactory& objectFactory();
}

#endif /* OBJECTFACTORY_HPP */
