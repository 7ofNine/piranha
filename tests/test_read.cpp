// test for reading a file

#include <string>
#include <fstream>
#include "piranha.h"


using namespace piranha;
using namespace piranha::manipulators;
using namespace std;

int main()
{
    // very simple just construct it
    //ofstream check("check.txt");
    //check << "check where is it" << endl;
    //check.close();

    std::string testfile("testxxxx.qps");
    qps test(testfile);
}