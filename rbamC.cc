#include <zlib.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <sstream>
#include <memory>
#include <cassert>
#include <algorithm>

#include <bamreader.hpp>

using namespace std;
using Sequence::bamreader;
using Sequence::bamrecord;
using Sequence::samflag;

int main( int argc, char ** argv )
{
  gzFile in = gzopen(argv[1],"rb");
  if(in == NULL) { exit(0);}
  
  bamreader reader( argv[1] );
  unsigned nread=0;
  std::vector<bamrecord> PAR,UL;
  //vector< pair<string,bamrecord > > vb2;
  PAR.reserve(50000000);
  UL.reserve(50000000);
  while(! reader.eof() && !reader.error() )
    {
      bamrecord b(reader.next_record());
      if(b.empty()) break;

      if( b.refid() != -1 && b.next_refid() != -1 ) //if both reads mapped
	{
	  if( b.refid() != b.next_refid() ) //both map to different scaffolds
	    {
	      UL.emplace_back(std::move(b));
	    }
	  else
	    {
	      samflag sf = b.flag();
	      assert(sf.is_paired);
	      if( sf.qstrand == sf.mstrand && b.pos() != b.next_pos())
		{
		  if( b.read_name() == "1:1:234846:0" ||
		      b.read_name() == "1:1:234846:1" )
		    {
		      cout << b.mapq() << ' ' << b.seq() << ' '
			   << b.aux() << '\n';
		    }
		  if( b.refid() ==  b.next_refid() ) {
		    PAR.emplace_back(std::move(b));
		  }
		}
	    }
	}
      ++nread;
    }

  //Sort the alignments by read name
  std::sort( PAR.begin(), PAR.end(),
	     [](const bamrecord & lhs, const bamrecord & rhs) {
	       return lhs.read_name() < rhs.read_name();
	     } );

  std::sort( UL.begin(), UL.end(),
	     [](const bamrecord & lhs, const bamrecord & rhs) {
	       return lhs.read_name() < rhs.read_name();
	     } );
  std::for_each(PAR.begin(),PAR.end(),
		[](const bamrecord & b) {
		  cout << b.read_name() << '\n';
		});
  cout << nread << " alignments processed\n";
}



