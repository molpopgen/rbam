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

using namespace std;

struct Xtra
{
  ostringstream data;
  char tag[2],val_type;
  //Xtra();
  int read(gzFile in,
	   const int & curr_pos,
	   const int & block_size);
};

/*
Xtra::Xtra() : data(nullptr)
{
}
*/
int Xtra::read(gzFile in,
	       const int & curr_pos,
	       const int & block_size)
{
  int rv=0;
  rv+=gzread(in,&tag[0],2);
  rv+=gzread(in,&val_type,1);

  //cout << tag << ' ' << val_type << ' ';
  if( val_type == 'A' ) {
    char temp;
    rv+=gzread(in,&temp,sizeof(char));
    data << temp;
  //‘cCsSiI’ in BAM, representing int8 t, uint8 t, int16 t, uint16 t, int32 t and
  //uint32 t, respectively. In SAM, all single integer types are mapped to int32 t.
  } else if ( val_type == 'c' ) {
    int8_t temp;
    rv+=gzread(in,&temp,sizeof(int8_t));
    data << temp;
  } else if ( val_type == 'C' ) {
    uint8_t temp;
    rv+=gzread(in,&temp,sizeof(uint8_t));
    data << temp;
  } else if ( val_type == 's' ) {
    int16_t temp;
    rv+=gzread(in,&temp,sizeof(int16_t));
    data << temp;
  } else if ( val_type == 'S' ) {
    uint16_t temp;
    rv+=gzread(in,&temp,sizeof(uint16_t));
    data << temp;
  } else if ( val_type == 'i' ) {
    int32_t temp;
    rv+=gzread(in,&temp,sizeof(int32_t));
    data << temp;
  } else if ( val_type == 'I' ) {
    uint32_t temp;
    rv+=gzread(in,&temp,sizeof(uint32_t));
    data << temp;
  } 
  //cout << data.str() << '\n';
  return rv;   
}

int main( int argc, char ** argv )
{
  gzFile in = gzopen(argv[1],"rb");
  if(in == NULL) { exit(0);}
  
  char magic[4];

  gzread( in, &magic[0], 4*sizeof(char) );

  //cout << magic << '\n';

  using U32 = std::uint32_t;
  using I32 = std::int32_t;

  I32 l_text,n_ref;

  gzread( in, &l_text, sizeof(I32) );

  //cout << l_text << '\n';
  char htext[l_text];

  gzread(in,&htext[0],l_text*sizeof(char));
  //cout << htext << '\n';

  gzread(in,&n_ref,sizeof(I32));

  //cout << n_ref << '\n';

  for( I32 i = 0 ; i < n_ref ; ++i )
    {
      I32 l_name,l_ref;
      gzread( in,&l_name,sizeof(I32) );
      char name[l_name];
      gzread( in,&name[0],l_name*sizeof(char));
      gzread( in,&l_ref,sizeof(I32) );
      //cout << l_name << ' ' << name << ' ' << l_ref << '\n';
    }
  unsigned nread=0;
  //while(!gzeof(in))
  do
    {

      I32 block_size,refID,pos,l_seq,next_refID,next_pos,tlen;
      U32 bin_mq_nl,flag_nc,BIN,MAPQ,l_read_name,flag,nc;
      
      int test = gzread(in,&block_size,sizeof(I32));
      if(!test || test == -1) break;//EOF or error
      ++nread;
      char block[block_size];
      //This works:
      /*
	char block[block_size];
	gzread(in,&block[0],block_size);
	cout << "block starting:\n";
	for( auto i=0;i<block_size;++i)
	{
	cout << block[i];
	}
  cout << "block ending\n";
  */
  int32_t read = 0;
  read+=gzread(in,&refID,sizeof(I32));
  read+=gzread(in,&pos,sizeof(I32));
  read+=gzread(in,&bin_mq_nl,sizeof(U32));
  read+=gzread(in,&flag_nc,sizeof(U32));
  read+=gzread(in,&l_seq,sizeof(I32));
  read+=gzread(in,&next_refID,sizeof(I32));
  read+=gzread(in,&next_pos,sizeof(I32));
  read+=gzread(in,&tlen,sizeof(I32));

  BIN = bin_mq_nl >> 16;
  MAPQ= (bin_mq_nl ^ (BIN<<16))>>8;
  l_read_name = bin_mq_nl ^ (BIN<<16) ^ (MAPQ<<8);
  
  flag = flag_nc >> 16;
  nc = flag_nc ^ (flag<<16);
  
  char read_name[l_read_name];
  read+=gzread(in,&read_name[0],l_read_name);

  //cout << read_name << '\n';

  U32 ncigop[nc];

  read+=gzread(in,&ncigop[0],nc*sizeof(U32));

  // for(auto i=0;i<nc;++i)
  //   {
  //     cout << ncigop[i]<<' ';
  //   }
  // cout << '\n';

  uint8_t seq[ (l_seq+1)/2 ];

  read+=gzread(in,&seq[0],((l_seq+1)/2)*sizeof(uint8_t));

  char qual[l_seq];

  read+=gzread(in,&qual[0],l_seq);

  // cout << flag << ' ' << seq << ' ';
  // for( auto q : qual )
  //   cout << char(q+32);
  // cout << ' ' << read << ' ' << block_size << '\n';
  //Read the xtra stuff

  if(read<block_size)
    {
      char temp[block_size-read];
      gzread(in,&temp[0],block_size-read);
    }
  // while(read < block_size)
  //   {
  //     Xtra x;
  //     read += x.read(in,read,block_size);
  //   }
    }
  while(!gzeof(in));
  cout << nread << " alignments processed\n";
}



