#include <bamrecord.hpp>
#include <bamutil.hpp>
#include <algorithm>
#include <cstring>

using namespace std;
using Sequence::samflag;

namespace {
  static const size_t lI32s = 6;
  static const size_t lU32s = 7;
  //! What stored in __rdata1 of a bamrecordImpl
  enum class I32s : size_t {REFID,POS,L_SEQ,NEXT_REFID,NEXT_POS,TLEN};
  //! What stored in __rdata2 of a bamrecordImpl
  enum class U32s : size_t {BIN_MQ_NL,FLAG_NC,BIN,MAPQ,L_READ_NAME,FLAG,NC};
  //14An integer may be stored as one of ‘cCsSiI’ in BAM, representing int8 t, uint8 t, int16 t, uint16 t, int32 t and
  //uint32 t, respectively. In SAM, all single integer types are mapped to int32 t
  static size_t sc = sizeof(int8_t),
    sC = sizeof(uint8_t),
    ss = sizeof(int16_t),
    sS = sizeof(uint16_t),
    si = sizeof(int32_t),
    sI = sizeof(uint32_t);
}

using U32 = std::uint32_t;
using I32 = std::int32_t;
using U16 = std::uint16_t;
using I16 = std::int16_t;
using U8 = std::uint8_t;
using I8 = std::int8_t;

namespace Sequence 
{
  //Implementation class details
  class bamrecordImpl
  {
  public:
    bool __empty;
    I32 __read,__block_size,__rdata1[6];
    U32 __rdata2[7];
    std::unique_ptr<char[]> __read_name,__qual,__aux;
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
    __empty=bI.__empty;
    __read=bI.__read;
    __block_size=bI.__block_size;
    __qual=nullptr;
    __aux=nullptr;
    __ncigop=nullptr;
    __seq=nullptr;
    if (! __empty )
      {
	//Copy data into current object
	std::copy( &bI.__rdata1[0],
		   &bI.__rdata1[lI32s],
		   &__rdata1[0] );
	std::copy( &bI.__rdata2[0],
		   &bI.__rdata2[lU32s],
		   &__rdata2[0] );
	//Allocate our local pointers
	__read_name = std::unique_ptr<char[]>(new char[__rdata2[static_cast<size_t>(U32s::L_READ_NAME)]]);
	__ncigop = std::unique_ptr<U32[]>( new U32[__rdata2[static_cast<size_t>(U32s::NC)]] );
	__seq = std::unique_ptr<U8[]>(new U8[(__rdata1[static_cast<size_t>(I32s::L_SEQ)]+1)/2]);
	__qual = std::unique_ptr<char[]>(new char[__rdata1[static_cast<size_t>(I32s::L_SEQ)]]);
	__aux = std::unique_ptr<char[]>(new char[strlen(bI.__aux.get())]);

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
	std::copy( &bI.__aux[0],
		   &bI.__aux[strlen(bI.__aux.get())],
		   &__aux[0] );
      }
    return *this;
  }

  bamrecordImpl::bamrecordImpl( gzFile in ) : __empty(false),
					      __read(0),
					      __read_name(nullptr),
					      __qual(nullptr),
					      __aux(nullptr),
					      __ncigop(nullptr),
					      __seq(nullptr)
  {
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
    __aux = std::unique_ptr<char[]>(new char[__block_size-__read]);
    __read += gzread( in,__aux.get(), __block_size - __read );
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
    std::swap(this->__impl,rhs.__impl);
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
    return &__impl->__seq[0] + (__impl->__rdata1[static_cast<size_t>(I32s::L_SEQ)]+1)/2;
  }

  std::string bamrecord::cigar() const 
  {
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
    return &__impl->__qual[0]+__impl->__rdata1[static_cast<size_t>(I32s::L_SEQ)];
  }

  samflag bamrecord::flag() const
  {
    return samflag( __impl->__rdata2[static_cast<size_t>(U32s::FLAG)]);
  }
 
  std::uint32_t bamrecord::mapq() const
  {
    return __impl->__rdata2[static_cast<size_t>(U32s::MAPQ)];
  }

  std::int32_t bamrecord::pos() const {
    return __impl->__rdata1[static_cast<size_t>(I32s::POS)];
  }

  std::int32_t bamrecord::next_pos() const {
    return __impl->__rdata1[static_cast<size_t>(I32s::NEXT_POS)];
  }

  std::int32_t bamrecord::refid() const {
    return __impl->__rdata1[static_cast<size_t>(I32s::REFID)];
  }

  std::int32_t bamrecord::next_refid() const {
    return __impl->__rdata1[static_cast<size_t>(I32s::NEXT_REFID)];
  }

  std::int32_t bamrecord::tlen() const {
    return __impl->__rdata1[static_cast<size_t>(I32s::TLEN)];
  }

  std::int32_t bamrecord::l_seq() const {
    return __impl->__rdata1[static_cast<size_t>(I32s::L_SEQ)];
  }

  const char * bamrecord::hasTag(const char * tag) {
    if(tag==NULL||tag==nullptr) return nullptr;
    const char * rv = strstr(__impl->__aux.get(),tag);
    return(rv==NULL) ? nullptr : rv;
  }

  const char * bamrecord::hasTag(const char * start, const char * tag) {
    if(start == NULL || start == nullptr || tag == NULL || tag == nullptr) return nullptr;
    const char * rv = strstr(start,tag);
    return(rv==NULL) ? nullptr : rv;
  }

  std::string bamrecord::aux() const
  {
    string rv;
    size_t start = 0,end=strlen(__impl->__aux.get());
    char tag[3];
    tag[2]='\0';
    while( start < end )
      {
	tag[0]=__impl->__aux[start++];
	tag[1]=__impl->__aux[start++];
	char val_type = __impl->__aux[start++];
	rv += string(tag) + ":" + val_type + ":";
	if ( val_type == 'A' )
	  {
	    rv += __impl->__aux[start++];
	  }
	else if ( val_type == 'c' )
	  {
	    I32 v = 0;
	    for( size_t i = 0 ; i < sizeof(int8_t) ; ++i )
	      {
		v |= __impl->__aux[start++];
	      }
	    rv += to_string(v);
	  }
	else if (val_type == 'C')
	  {
	    I32 v = 0;
	    for( size_t i = 0 ; i < sizeof(uint8_t) ; ++i )
	      {
		v |= __impl->__aux[start++];
	      }
	    rv += to_string(v);
	  }
	else if (val_type == 's')
	  {
	    I32 v = 0;
	    for( size_t i = 0 ; i < sizeof(int16_t) ; ++i )
	      {
		v |= __impl->__aux[start++];
	      }
	    rv += to_string(v);
	  }
	else if (val_type == 'S')
	  {
	    I32 v = 0;
	    for( size_t i = 0 ; i < sizeof(uint16_t) ; ++i )
	      {
		v |= __impl->__aux[start++];
	      }
	    rv += to_string(v);
	  }
	else if (val_type == 'i')
	  {
	    I32 v = 0;
	    for( size_t i = 0 ; i < sizeof(int32_t) ; ++i )
	      {
		v |= __impl->__aux[start++];
	      }
	    rv += to_string(v);
	  }
	else if (val_type == 'I')
	  {
	    I32 v = 0;
	    for( size_t i = 0 ; i < sizeof(uint32_t) ; ++i )
	      {
		v |= __impl->__aux[start++];
	      }
	    rv += to_string(v);
	  }
	else if (val_type == 'B')//not impl
	  {
	  }
	else if (val_type == 'Z')
	  {
	    char * __x = std::find(&__impl->__aux[start],&__impl->__aux[end],'\0');
	    rv += string( &__impl->__aux[start],__x );
	    start += (__x - &__impl->__aux[start])+1;
	  }
	rv += '\t';
      }
    return rv;
  }

}
