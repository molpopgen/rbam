#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cstdlib>
#include <cstdint>

using namespace std;

int main(int argc, char ** argv)
{
  using U32 = std::uint32_t;

  U32 bin=555,mapq=255,l_read_name=21;

  U32 bin_mq_nl = (l_read_name|bin<<16|mapq<<8);

  cout << bin << '\n'
       << mapq << '\n'
       << l_read_name << '\n'
       << bin_mq_nl << '\n';

  //U32 A,B,C;
  int32_t A;
  int32_t B,C;

  A = bin_mq_nl >> 16;
  B = (bin_mq_nl ^ (A<<16))>>8;
  C = bin_mq_nl ^ (A<<16) ^ (B<<8);
  cout << bin_mq_nl << '\n'
       << A << '\n'
       << B << '\n'
       << C << '\n';
}
