#include "Tester.hpp"

using namespace broomstyx;

Tester::Tester()
    : _testsPassed(0)
    , _testsFailed(0)
{}

Tester::Tester( std::string title )
    : _testsPassed(0)
    , _testsFailed(0)
{
    _title = title;
}

Tester::~Tester() {}


bool Tester::checkThat( bool expr )
{
    if ( expr )
    {
        _testsPassed += 1;
        return true;
    }
    else
    {
        _testsFailed += 1;
        return false;
    }
}

void Tester::reportResults()
{
    if ( _testsFailed == 0 )
    {
        std::printf( " ---------------------------------------------------------------------------\n" );
        std::printf( "  %-30s : %3d check(s) passed,    none failed\n", _title.c_str(), _testsPassed );
    }
    else
    {
        std::printf( " ---------------------------------------------------------------------------\n" );
        std::printf( "  %-30s : %3d check(s) passed, %3d check(s) failed  <--- ATTENTION!!! \n", _title.c_str(), _testsPassed, _testsFailed );
    }
}
