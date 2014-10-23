#include <bamrecord.hpp>
#include <bamutil.hpp>
#include <algorithm>

using namespace std;
using Sequence::samflag;

//Implementation class details
class bamrecordImpl
{
public:
  using U32 = std::uint32_t;
  using I32 = std::int32_t;
  using U8 = std::uint8_t;

  bool __empty;
  I32 __read,__block_size,__rdata1[6];
  U32 __rdata2[7];
  std::unique_ptr<char[]> __read_name,__qual,__rest;
  std::unique_ptr<U32[]> __ncigop;
  std::unique_ptr<U8[]> __seq;

  bamrecordImpl( gzFile in );
  bamrecordImpl( const bamrecordImpl & );
  bamrecordImpl& operator=( const bamrecordImpl &);
};

bamrecordImpl::bamrecordImpl( const bamrecordImpl & bI) {
  *this = bI;
}

bamrecordImpl& bamrecordImpl::operator=( const bamrecordImpl & bI)  
{
  using bamutil::U32s;
  using bamutil::I32s;
  __empty=bI.__empty;
  __read=bI.__read;
  __block_size=bI.__block_size;
  __qual=nullptr;
  __rest=nullptr;
  __ncigop=nullptr;
  __seq=nullptr;
  if (! __empty )
    {
      //Copy data into current object
      std::copy( &bI.__rdata1[0],
		 &bI.__rdata1[bamutil::lI32s],
		 &__rdata1[0] );
      std::copy( &bI.__rdata2[0],
		 &bI.__rdata2[bamutil::lU32s],
		 &__rdata2[0] );
      //Allocate our local pointers
      __read_name = std::unique_ptr<char[]>(new char[__rdata2[static_cast<size_t>(U32s::L_READ_NAME)]]);
      __ncigop = std::unique_ptr<U32[]>( new U32[__rdata2[static_cast<size_t>(U32s::NC)]] );
      __seq = std::unique_ptr<U8[]>(new U8[(__rdata1[static_cast<size_t>(I32s::L_SEQ)]+1)/2]);
      __qual = std::unique_ptr<char[]>(new char[__rdata1[static_cast<size_t>(I32s::L_SEQ)]]);
      __rest = std::unique_ptr<char[]>(new char[__block_size-__read]);

      //And copy
      std::copy( &bI.__read_name[0],
		 &bI.__read_name[__rdata2[static_cast<size_t>(U32s::L_READ_NAME)]],
		 &__read_name[0] );
      std::copy( &bI.__ncigop[0],
		 &bI.__ncigop[__rdata2[static_cast<size_t>(U32s::NC)]],
		 &__ncigop[0] );
      std::copy( &bI.__seq[0],
		 &bI.__seq[(__rdata1[static_cast<size_t>(I32s::L_SEQ)]+1)/2],
		 &__seq[0] );
      std::copy( &bI.__qual[0],
		 &bI.__qual[__rdata1[static_cast<size_t>(I32s::L_SEQ)]],
		 &__qual[0] );
      std::copy( &bI.__rest[0],
		 &bI.__rest[__block_size-__read],
		 &__rest[0] );
    }
  return *this;
}

bamrecordImpl::bamrecordImpl( gzFile in ) : __empty(false),
				    __read(0),
				    __read_name(nullptr),
				    __qual(nullptr),
				    __rest(nullptr),
				    __ncigop(nullptr),
				    __seq(nullptr)
{
  using bamutil::I32s;
  using bamutil::U32s;
  int test = gzread(in,&__block_size,sizeof(I32));
  if(!test || test == -1) { 
    __empty = true;
    return;
  }
  __read += gzread(in,&__rdata1[static_cast<size_t>(I32s::REFID)],2*sizeof(I32));
  __read += gzread(in,&__rdata2[static_cast<size_t>(U32s::BIN_MQ_NL)],2*sizeof(I32));
  __read += gzread(in,&__rdata1[static_cast<size_t>(I32s::L_SEQ)],4*sizeof(I32));

  __rdata2[static_cast<size_t>(U32s::BIN)] = __rdata2[static_cast<size_t>(U32s::BIN_MQ_NL)] >> 16;
  __rdata2[static_cast<size_t>(U32s::MAPQ)] = (__rdata2[static_cast<size_t>(U32s::BIN_MQ_NL)] ^ __rdata2[static_cast<size_t>(U32s::BIN)]<<16)>>8;
  __rdata2[static_cast<size_t>(U32s::L_READ_NAME)] = __rdata2[static_cast<size_t>(U32s::BIN_MQ_NL)] ^ __rdata2[static_cast<size_t>(U32s::BIN)]<<16 ^ __rdata2[static_cast<size_t>(U32s::MAPQ)]<<8;
  __rdata2[static_cast<size_t>(U32s::FLAG)] = __rdata2[static_cast<size_t>(U32s::FLAG_NC)]>>16;
  __rdata2[static_cast<size_t>(U32s::NC)] = __rdata2[static_cast<size_t>(U32s::FLAG_NC)] ^ (__rdata2[static_cast<size_t>(U32s::FLAG)]<<16);

  __read_name = std::unique_ptr<char[]>(new char[__rdata2[static_cast<size_t>(U32s::L_READ_NAME)]]);
  __read += gzread( in,__read_name.get(),__rdata2[static_cast<size_t>(U32s::L_READ_NAME)]*sizeof(char) );

  __ncigop = std::unique_ptr<U32[]>( new U32[__rdata2[static_cast<size_t>(U32s::NC)]] );
  __read += gzread( in, __ncigop.get(),__rdata2[static_cast<size_t>(U32s::NC)]*sizeof(U32) );

  __seq = std::unique_ptr<U8[]>(new U8[(__rdata1[static_cast<size_t>(I32s::L_SEQ)]+1)/2]);
  __read += gzread( in, __seq.get(), ((__rdata1[static_cast<size_t>(I32s::L_SEQ)]+1)/2)*sizeof(U8) );

  __qual = std::unique_ptr<char[]>(new char[__rdata1[static_cast<size_t>(I32s::L_SEQ)]]);
  __read += gzread( in, __qual.get(),__rdata1[static_cast<size_t>(I32s::L_SEQ)]*sizeof(char) );

  //HACK
  __rest = std::unique_ptr<char[]>(new char[__block_size-__read]);
  __read += gzread( in,__rest.get(), __block_size - __read );
}

bamrecord::bamrecord( gzFile in ) : __impl( std::unique_ptr<bamrecordImpl>(new bamrecordImpl(in)) ) {}

bamrecord::bamrecord( const bamrecord & rhs) : __impl(std::unique_ptr<bamrecordImpl>(new bamrecordImpl(*rhs.__impl))) {}

bool bamrecord::empty() const 
{
  return __impl->__empty;
}

bamrecord::bamrecord(bamrecord&&rhs) : __impl(std::move(rhs.__impl))
{
  rhs.__impl = nullptr;
}

bamrecord & bamrecord::operator=(bamrecord&&rhs) 
{
  bamrecordImpl  * t = this->__impl.release();
  delete t;
  //this->__impl = nullptr;
  std::swap(this->__impl,rhs.__impl);
  //rhs.__impl.release();
  /*
  this->__impl.release();
  this->__impl = std::move(rhs.__impl);
  rhs.__impl = nullptr;
  */
  return *this;
}

bamrecord & bamrecord::operator=(const bamrecord &rhs )
{
  bamrecordImpl * t = this->__impl.release();
  delete t;
  this->__impl = std::unique_ptr<bamrecordImpl>(new bamrecordImpl(*rhs.__impl));
  return *this;
}

bamrecord::~bamrecord() {

}

std::string bamrecord::read_name() const
{
  return std::string(__impl->__read_name.get());
}

std::string bamrecord::seq() const 
{
  using bamutil::I32s;
  std::unique_ptr<char[]> rv(new char[__impl->__rdata1[static_cast<size_t>(I32s::L_SEQ)]+1]);
  size_t __pos = 0;
  std::for_each(this->seq_cbeg(),
		this->seq_cend(),
		[&](const U8 & i) {
		  rv[__pos++]= bamutil::int2seq[ (((i) >> 4) & 0x0F) ];
		  rv[__pos++]= bamutil::int2seq[ (((i)) & 0x0F) ];
		});
  rv[__pos]='\0';
  return(string(rv.get()));
 }

const std::uint8_t * bamrecord::seq_cbeg() const
{
  return &__impl->__seq[0];
}

const std::uint8_t * bamrecord::seq_cend() const
{
  using bamutil::I32s;
  return &__impl->__seq[0] + (__impl->__rdata1[static_cast<size_t>(I32s::L_SEQ)]+1)/2;
}

std::string bamrecord::cigar() const 
{
  using bamutil::U32s;
  std::string rv;
  std::for_each( &__impl->__ncigop[0],
		 &__impl->__ncigop[0]+ __impl->__rdata2[static_cast<size_t>(U32s::NC)],
		 [&](const U32 & i )
		 {
		   U32 __op = (i>>4);
		   rv += to_string(__op) + bamutil::bamCig[i^(__op<<4)];
		 } );
  return rv;
}

std::string bamrecord::qual() const 
{
  return std::string(this->qual_cbeg(),
		     this->qual_cend());
}

const char * bamrecord::qual_cbeg() const 
{
  return &__impl->__qual[0];
}

const char * bamrecord::qual_cend() const 
{
  using bamutil::I32s;
  return &__impl->__qual[0]+__impl->__rdata1[static_cast<size_t>(I32s::L_SEQ)];
}

samflag bamrecord::flag() const
{
  using bamutil::U32s;
  return samflag( __impl->__rdata2[static_cast<size_t>(U32s::FLAG)]);
}
 
std::uint32_t bamrecord::mapq() const
{
  using bamutil::U32s;
  return __impl->__rdata2[static_cast<size_t>(U32s::MAPQ)];
}

std::int32_t bamrecord::pos() const {
  using bamutil::I32s;
  return __impl->__rdata1[static_cast<size_t>(I32s::POS)];
}

std::int32_t bamrecord::next_pos() const {
  using bamutil::I32s;
  return __impl->__rdata1[static_cast<size_t>(I32s::NEXT_POS)];
}

std::int32_t bamrecord::refid() const {
  using bamutil::I32s;
  return __impl->__rdata1[static_cast<size_t>(I32s::REFID)];
}

std::int32_t bamrecord::next_refid() const {
  using bamutil::I32s;
  return __impl->__rdata1[static_cast<size_t>(I32s::NEXT_REFID)];
}

std::int32_t bamrecord::tlen() const {
  using bamutil::I32s;
  return __impl->__rdata1[static_cast<size_t>(I32s::TLEN)];
}
