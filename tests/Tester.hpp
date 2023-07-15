#ifndef TESTER_HPP
#define TESTER_HPP

#include <string>

namespace broomstyx
{
    class Tester final
    {
    public:
        Tester();
        Tester( std::string title );
        virtual ~Tester();

        // Disable copy constructor and assignment operator
        Tester( const Tester& ) = delete;
        Tester& operator=( const Tester& ) = delete;

        bool checkThat( bool expr );
        void reportResults();

    private:
        std::string _title;
        int _testsPassed;
        int _testsFailed;
    };
}


#endif // TESTER_HPP
