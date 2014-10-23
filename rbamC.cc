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
#include <algorithm>

#include <bamreader.hpp>

using namespace std;
using namespace Sequence;

int main( int argc, char ** argv )
{
  gzFile in = gzopen(argv[1],"rb");
  if(in == NULL) { exit(0);}
  
  bamreader reader( argv[1] );
  unsigned nread=0;
  std::vector<bamrecord> vb;
  vector< pair<string,bamrecord > > vb2;
  vb.reserve(50000000);
  while(! reader.eof() && !reader.error() )
    {
      bamrecord b(reader.next_record());
      if(b.empty()) break;
      //cout << b.aux() << '\n';
      //auto s = b.seq(),c=b.cigar();//,q=b.qual();
      //cout << s << ' ' << c << '\n';
	  //cout << s << '\n';
	  /*
	  auto s = b.seq(),q=b.qual();
	  cout << s << ' ';
	  for( auto __q : q ) cout << char(__q+32);
	  cout << '\n';
	  */
      //cout << b.pos() << ' ' << b.next_pos() << ' ' << b.refid() << ' ' << b.next_refid() << ' ' << b.tlen() << '\n';
      ++nread;
	  //if(vb.size()<50000000) vb.emplace_back(std::move(b)); //calls move constructor for bamrecord && (faster)
	  //vb2.push_back( make_pair(b.read_name(),b) );
	  /*
	  for(unsigned i=0;i<vb2.size();++i)
	    {
	      bamrecord b2 = vb2[i].second;
	      cout << b2.seq() << '\n';
	      //cout << vb2[i].first << ' ' << vb2[i].second.seq() << '\n';
	    }
	  */
	  //if(vb.size()<50000000) vb.push_back(b); //calls copy constructor for const bamrecord &
	  //else
	  //exit(1);

	  // for( auto rec:vb )
	  //   {
	  //     auto s = rec.seq(),
	  // 	c=rec.cigar(),
	  // 	q=rec.qual();
	  //     cout << s << ' ' << c << ' ' << q << '\n';
	  //     exit(1);
	  //   }
	  //if(nread<10000000) vb.emplace_back(std::move(b));
	  /*
	  else
	    {
	      for_each(vb.begin(),vb.end(),
		       [](const bamrecord & rec)
		       {
			 auto s = rec.seq(),
			   c=rec.cigar(),
			   q=rec.qual();
			 cout << s << ' ' << c << ' ' << q << '\n';
		       }
		       );
	      exit(1);
	    }
	  */
    }
  cout << nread << " alignments processed\n";
}



