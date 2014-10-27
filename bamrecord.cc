#include <bamutil.hpp>
#include <bamrecord.hpp>
#include <algorithm>
#include <cstring>
#include <cassert>

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

  size_t auxTagSize(const char & c)
  {
    switch(c)
      {
      case 'c':
	return sizeof(int8_t);
	break;
      case 'C':
	return sizeof(uint8_t);
	break;
      case 's':
	return sizeof(int16_t);
	break;
      case 'S':
	return sizeof(uint16_t);
	break;
      case 'i':
	return sizeof(int32_t);
	break;
      case 'I':
	return sizeof(uint32_t);
	break;
      default:
	return 0;
	break;
      }
    return 0;
  }
  /*
  static size_t sc = sizeof(int8_t),
    sC = sizeof(uint8_t),
    ss = sizeof(int16_t),
    sS = sizeof(uint16_t),
    si = sizeof(int32_t),
    sI = sizeof(uint32_t);
  */
}

using U32 = std::uint32_t;
using I32 = std::int32_t;
using U16 = std::uint16_t;
using I16 = std::int16_t;
using U8 = std::uint8_t;
using I8 = std::int8_t;

namespace Sequence 
{
  using bamutil::int2seq;
  using bamutil::bamCig;

  bamaux::bamaux( ) : size(0),
		      value_type(char()),
		      tag(nullptr),
		      value(nullptr)
  {
  }

  bamaux::bamaux( size_t __size,
		  std::unique_ptr<char[]> & __tag,
		  char __value_type,
		  std::unique_ptr<unsigned char[]> & __value) : size(std::move(__size)),
							 value_type(std::move(__value_type)),
							 tag(std::move(__tag)),
								value(std::move(__value))
  {
  }

  bamaux::bamaux( bamaux && ba ) : size(std::move(ba.size)),
				   value_type(std::move(ba.value_type)),
				   tag(std::move(ba.tag)),
				   value(std::move(ba.value))
  {
    //value = 
  }

  rawbam::rawbam( size_t && __block_size,
		  std::unique_ptr<char[]> && __block) : block_size(move(__block_size)),
							block(move(__block))
  {
    //To be implemented
  }

  //Implementation class details
  class bamrecordImpl
  {
  public:
    bool __empty;
    I32 __read,__block_size,__rdata1[6];
    U32 __rdata2[7];
    size_t __aux_size;
    std::unique_ptr<char[]> __read_name,__qual,__aux;
    std::unique_ptr<U32[]> __ncigop;
    std::unique_ptr<U8[]> __seq;

    bamrecordImpl(  );
    bamrecordImpl( BGZF * in );
    bamrecordImpl( const bamrecordImpl & );
    bamrecordImpl& operator=( const bamrecordImpl &);
  };

  bamrecordImpl::bamrecordImpl() : __empty(true),
				   __aux_size(0),
				   __read_name(nullptr),
				   __qual(nullptr),
				   __aux(nullptr),
				   __ncigop(nullptr),
				   __seq(nullptr)
  {
  }
  
  bamrecordImpl::bamrecordImpl( const bamrecordImpl & bI) {
    *this = bI;
  }

  bamrecordImpl& bamrecordImpl::operator=( const bamrecordImpl & bI)  
  {
    __empty=bI.__empty;
    __read=bI.__read;
    __block_size=bI.__block_size;
    __aux_size = bI.__aux_size;
    __read_name=nullptr;
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
	if(bI.__aux_size)
	  __aux = std::unique_ptr<char[]>(new char[bI.__aux_size]);

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
	if(bI.__aux_size) {
	  std::copy( &bI.__aux[0],
		     &bI.__aux[0]+bI.__aux_size,
		     &__aux[0] );
	}
      }
    return *this;
  }

  bamrecordImpl::bamrecordImpl( BGZF * in ) : __empty(false),
					      __read(0),
					      __aux_size(0),
					      __read_name(nullptr),
					      __qual(nullptr),
					      __aux(nullptr),
					      __ncigop(nullptr),
					      __seq(nullptr)
  {
    int test = bgzf_read(in,&__block_size,sizeof(I32));
    if(!test || test == -1) { 
      __empty = true;
      return;
    }
    __read += bgzf_read(in,&__rdata1[static_cast<size_t>(I32s::REFID)],2*sizeof(I32));
    __read += bgzf_read(in,&__rdata2[static_cast<size_t>(U32s::BIN_MQ_NL)],2*sizeof(I32));
    __read += bgzf_read(in,&__rdata1[static_cast<size_t>(I32s::L_SEQ)],4*sizeof(I32));

    __rdata2[static_cast<size_t>(U32s::BIN)] = __rdata2[static_cast<size_t>(U32s::BIN_MQ_NL)] >> 16;
    __rdata2[static_cast<size_t>(U32s::MAPQ)] = (__rdata2[static_cast<size_t>(U32s::BIN_MQ_NL)] ^ __rdata2[static_cast<size_t>(U32s::BIN)]<<16)>>8;
    __rdata2[static_cast<size_t>(U32s::L_READ_NAME)] = __rdata2[static_cast<size_t>(U32s::BIN_MQ_NL)] ^ __rdata2[static_cast<size_t>(U32s::BIN)]<<16 ^ __rdata2[static_cast<size_t>(U32s::MAPQ)]<<8;
    __rdata2[static_cast<size_t>(U32s::FLAG)] = __rdata2[static_cast<size_t>(U32s::FLAG_NC)]>>16;
    __rdata2[static_cast<size_t>(U32s::NC)] = __rdata2[static_cast<size_t>(U32s::FLAG_NC)] ^ (__rdata2[static_cast<size_t>(U32s::FLAG)]<<16);

    __read_name = std::unique_ptr<char[]>(new char[__rdata2[static_cast<size_t>(U32s::L_READ_NAME)]]);
    __read += bgzf_read( in,__read_name.get(),__rdata2[static_cast<size_t>(U32s::L_READ_NAME)]*sizeof(char) );

    __ncigop = std::unique_ptr<U32[]>( new U32[__rdata2[static_cast<size_t>(U32s::NC)]] );
    __read += bgzf_read( in, __ncigop.get(),__rdata2[static_cast<size_t>(U32s::NC)]*sizeof(U32) );

    __seq = std::unique_ptr<U8[]>(new U8[(__rdata1[static_cast<size_t>(I32s::L_SEQ)]+1)/2]);
    __read += bgzf_read( in, __seq.get(), ((__rdata1[static_cast<size_t>(I32s::L_SEQ)]+1)/2)*sizeof(U8) );

    __qual = std::unique_ptr<char[]>(new char[__rdata1[static_cast<size_t>(I32s::L_SEQ)]]);
    __read += bgzf_read( in, __qual.get(),__rdata1[static_cast<size_t>(I32s::L_SEQ)]*sizeof(char) );
    if( __read < __block_size )
      {
	__aux_size = __block_size-__read;
	__aux = std::unique_ptr<char[]>(new char[__aux_size]);
	__read += bgzf_read( in,__aux.get(), __aux_size );
      }
    assert(__read == __block_size);
  }

  bamrecord::bamrecord(  ) : __impl( std::unique_ptr<bamrecordImpl>(new bamrecordImpl()) ) {}

  bamrecord::bamrecord( BGZF * in ) : __impl( std::unique_ptr<bamrecordImpl>(new bamrecordImpl(in)) ) {}

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
		    rv[__pos++]= int2seq[ (((i) >> 4) & 0x0F) ];
		    rv[__pos++]= int2seq[ (((i)) & 0x0F) ];
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
		     rv += to_string(__op) + bamCig[i^(__op<<4)];
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

  rawbam bamrecord::raw() const
  {
    size_t bsize = __impl->__block_size;
    std::unique_ptr<char[]> data;

    return rawbam(std::move(bsize),std::move(data));
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

  const char * bamrecord::hasTag(const char * tag) const {
    if(tag==NULL||tag==nullptr || !__impl->__aux_size) {
      return nullptr;
    }
    const char * rv = __impl->__aux.get(),
      * end = (__impl->__aux.get() + __impl->__aux_size);
    
    for( ; rv < end-1 ; rv++)
      {
	if( *rv == tag[0] &&
	    *(rv+1) == tag[1]) return rv;
      }
    return nullptr;
  }
  
  const char * bamrecord::hasTag(const char * start, const char * tag) const {
    if(start == NULL || start == nullptr || tag == NULL || tag == nullptr || !__impl->__aux_size) return nullptr;
    const char * rv = strstr(start,tag);
    return(rv==NULL) ? nullptr : rv;
  }

  bamaux bamrecord::aux(const char * tag) const {
    const char * tagspot = this->hasTag(tag);
    if( tagspot == nullptr ) return bamaux();
    std::unique_ptr<char[]> __tag(new char[3]);
    __tag[0]=*tagspot;
    __tag[1]=*(tagspot+1);
    __tag[2] = '\0';
    char __val_type = *(tagspot+2);
    //cout << "check: "<< tag << ' ' << __val_type << '\n';
    size_t valsize = (__val_type == 'A' ) ? 2*sizeof(char) : auxTagSize(__val_type);
    if(!valsize && !(__val_type=='Z'||__val_type=='B'))  return bamaux(); //something goofed
    std::unique_ptr<unsigned char[]> value;
    if(__val_type == 'Z')
      {
	const char * EOS = std::find( tagspot+2, reinterpret_cast<const char*>(&__impl->__aux[0]) + __impl->__aux_size, '\0' );
	valsize = EOS - (tagspot+2) + 1;
	value = std::unique_ptr<unsigned char[]>(new unsigned char[valsize]);
	copy( tagspot + 3, EOS+1, value.get() );
      }
    else if (__val_type == 'B')
      {
      }
    else
      {
	//Copy the tag value
	value = std::unique_ptr<unsigned char[]>(new unsigned char[valsize]);
	std::copy( tagspot+3, tagspot+3+valsize-1, value.get() );
	value[valsize-1]='\0';
      }
    return bamaux(valsize,__tag,__val_type,value);
  }

  std::string bamrecord::allaux() const
  {
    if(!__impl->__aux_size) return std::string();
    string rv;
    size_t start = 0,end=__impl->__aux_size;

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
