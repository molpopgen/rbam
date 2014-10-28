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
#include <chrono>

#include <bamreader.hpp>

using namespace std;
using Sequence::bamreader;
using Sequence::bamrecord;
using Sequence::bamaux;
using Sequence::samflag;

using readbucket = map< string,pair<bamrecord,bamrecord> >;
using Mbucket = map<string,z_off_t>;

void addRead( readbucket & rb, bamrecord & b, bool skipSecond = false )
{
  auto n = b.read_name();
  n.erase(n.end()-2,n.end());
  auto i = rb.find(n);
  if(i == rb.end())
    {
      rb.insert( make_pair(n, make_pair(std::move(b),bamrecord())) );
    }
  else if (!skipSecond)
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
  bamreader reader( argv[1] );
  //We have now processed all the header info, so let's record where we are
  auto pos = reader.tell();
  unsigned nread=0;
  readbucket PAR,UL,DIV;
  Mbucket Mpos; //locations of repetitive reads in the file
  while(! reader.eof() && !reader.error() )
  //while(1)
    {
      auto recstart = reader.tell();
      bamrecord b(reader.next_record());
      if(b.empty()) {  break; }
      ++nread;
      samflag sf(b.flag());
      if(!sf.query_unmapped  && !sf.mate_unmapped) //Both reads are mapped
      	{
	  bamaux ba = b.aux("XT");
	  if(ba.value[0]=='U') //Read is flagged as uniquely-mapping
	    {
	      if( b.refid() != b.next_refid() ) //both map to different scaffolds
	  	{
	  	  addRead(UL,b);
	  	}
	      else if ( b.pos() != b.next_pos()) //Don't map to same position
	  	{
	  	  if( sf.qstrand == sf.mstrand )
	  	    {
	  	      addRead(PAR,b);
	  	    }
	  	  else if( (sf.qstrand == 0 && b.pos() > b.next_pos()) ||
	  		   (sf.mstrand == 0 && b.next_pos() > b.pos() ) )
	  	    {
	  	      addRead(DIV,b);
	  	    }
	  	}
	    } 
	  else if (ba.value[0] == 'M') //Read is flagged as repetitively-mapping
	    {
	      //We pass true here so that we do not record multi/multi mappings
	      //addRead(UM,b,true);
	      auto n = b.read_name();
	      n.erase(n.end()-2,n.end());
	      auto i = Mpos.find(n);
	      if( i == Mpos.end() )
	  	{
	  	  //This is a rep. read that we've not seen before
	  	  Mpos.insert(make_pair(move(n),move(recstart)));
	  	}
	      else
	  	{
	  	  Mpos.erase(i);
	  	}
	    }
	}
    }
  unsigned NPAR=countReads(PAR),
    NUL=countReads(UL),
    NDIV=countReads(DIV);
  
  cout << nread << " alignments processed\n"
       << NPAR << " PAR reads\n"
       << NUL << " UL reads\n"
       << NDIV << " DIV reads\n";
  exit(0);
  //Pass2: find the U reads that go with any M reads
  //This is fucking slow.
  cerr << "beginning pass2\n";
  //reset reader to start of reads (see above)
  //std::cerr << reader.tell() << ' ';
  auto rv = reader.seek( pos, SEEK_SET );
  //std::cerr << reader.tell() << ' ' << rv << '\n';
  nread=0;
  unsigned NUM=0;
  //std::chrono::system_clock clock;
  //auto start = clock.now();
  //while(! reader.eof() && !reader.error() )
  while(1)
    {
      bamrecord b(reader.next_record());
      if(b.empty()) break;
      ++nread;
      /*
      if(nread && nread % 10000 == 0.) 
	{
	  auto now = clock.now();
	  auto elapsed = now-start;
	  start = now;
	  std::cout << elapsed.count() << ' '
	  << nread << ' ' << NUM << '\n';
	  }
      */
      samflag r(b.flag());
      if(!r.query_unmapped)
	{
	  bamaux ba = b.aux("XT");
	  if(ba.value[0]=='U') //Read is flagged as uniquely-mapping
	    {
	      auto n = b.read_name();
	      //cerr << n << '\n';
	      n.erase(n.end()-2,n.end());
	      auto i = Mpos.find(n);
	      if(i != Mpos.end()) //then the Unique reads redundant mate exists
		{
		  //get the multiple read now
		  bamrecord mate = reader.record_at_pos(i->second);
#ifndef NDEBUG
		  auto n2=mate.read_name();
		  n2.erase(n2.end()-2,n2.end());
		  assert(n == n2);
#endif
		  ++NUM;
		  Mpos.erase(i);
		}
	    }
	}
    }
  
  // unsigned NPAR=countReads(PAR),
  //   NUL=countReads(UL),
  //   NDIV=countReads(DIV);
  
  // cout << nread << " alignments processed\n"
  //      << NPAR << " PAR reads\n"
  //      << NUL << " UL reads\n"
  //      << NDIV << " DIV reads\n"
  //      << NUM << " UM reads\n";
  
  // for(auto i = DIV.cbegin(); i!=DIV.cend();++i)
  //   {
  //     if(!i->second.second.empty())
  // 	cout << i->first << '\n';
  //   }
}



