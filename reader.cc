#include <Sequence/samrecord.hpp>
#include <iostream>

using namespace std;
using namespace Sequence;

int main(int argc, char ** argv)
{
  samrecord r;
  while(!cin.eof())
    {
      cin >> r >> ws;
      samflag f = r.flag();
      cout << f << ' '
	   << f.query_unmapped << ' '
	   << f.mate_unmapped << ' '
	   << (!f.query_unmapped && !f.mate_unmapped) << '\n';
    }
}
