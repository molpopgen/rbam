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
using Sequence::bamaux;
using Sequence::samflag;

using readbucket = map< string,pair<bamrecord,bamrecord> >;

void addRead( readbucket & rb, bamrecord & b )
{
  auto n = b.read_name();
  n.erase(n.end()-2,n.end());
  auto i = rb.find(n);
  if(i == rb.end())
    {
      rb.insert( make_pair(n, make_pair(std::move(b),bamrecord())) );
    }
  else
    {
      i->second.second = std::move(b);
      assert(!i->second.second.empty());
    }
}

unsigned countReads(const readbucket & rb)
{
  unsigned rv=0;
  for( auto i = rb.cbegin() ; i != rb.cend() ; ++i )
    {
      if(! i->second.second.empty()) ++rv;
    }
  return rv;
}

int main( int argc, char ** argv )
{
  gzFile in = gzopen(argv[1],"rb");
  if(in == NULL) { exit(0);}
  
  bamreader reader( argv[1] );
  unsigned nread=0;
  readbucket PAR,UL,DIV,UM;
  while(! reader.eof() && !reader.error() )
    {
      bamrecord b(reader.next_record());
      if(b.empty()) break;
       ++nread;
       samflag sf = b.flag();
      //if( b.refid() != -1 && b.next_refid() != -1 ) //if both reads mapped
      //RELEARNING OLD LESSONS HERE: SAMTOOLS AND/OR BWA GOOF
      if(!sf.query_unmapped  && !sf.mate_unmapped)
      	{
	  bamaux ba = b.aux("XT");
	  if(ba.value[0]=='U') //Read is flagged as uniquely-mapping
	    {
	      if( b.refid() != b.next_refid() ) //both map to different scaffolds
		{
		  addRead(UL,b);
		}
	      else
		{
		  if( sf.qstrand == sf.mstrand && b.pos() != b.next_pos()) //On same strand and chromo
		    {
		      addRead(PAR,b);
		    }
		  else if( (!sf.qstrand && b.pos() > b.next_pos()) ||
		       (!sf.mstrand && b.next_pos() > b.pos()) ) //DIV
		    {
		      addRead(DIV,b);
		    }
		}
	    }
	}
    }
       
  unsigned NPAR=countReads(PAR),
    NUL=countReads(UL),
    NDIV=countReads(DIV);
  /*
  cout << nread << " alignments processed\n"
       << NPAR << " PAR reads\n"
       << NUL << " UL reads\n"
       << NDIV << " DIV reads\n";
  */
  for(auto i = DIV.cbegin(); i!=DIV.cend();++i)
    {
      cout << i->first << '\n';
    }
}



