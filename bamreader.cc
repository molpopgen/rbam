#include <bamreader.hpp>
#include <htslib/bgzf.h>

using std::string;

namespace Sequence
{

  //!Impl class for bamreader
  class bamreaderImpl
  {
  public:
    using I32 = std::int32_t;
    BGZF * in;
    bool __errorstate;
    char __magic[4];
    I32 __l_text,__n_ref;
    std::unique_ptr<char[]> __htext;
    std::vector< std::pair<std::string,I32> > __refdata;

    bamreaderImpl(const char * bamfilename);
    ~bamreaderImpl();
  }; 

  bamreaderImpl::~bamreaderImpl() 
  {
    if(in != NULL) bgzf_close(in);
  }

  bamreaderImpl::bamreaderImpl(const char * bamfilename) :
    in((bamfilename != nullptr) ? bgzf_open(bamfilename,"rb") : NULL),
    __errorstate(false),
    __htext(nullptr),
    __refdata(std::vector< std::pair<std::string,I32> >())
  {
    if(gzopen != NULL)
      {
	bgzf_read( in, &__magic[0], 4*sizeof(char) );
	if(string({__magic[0],__magic[1],__magic[2]}) != string("BAM")) __errorstate = 1;
	if(!__errorstate)
	  {
	    bgzf_read( in, &__l_text, sizeof(I32) );
	    __htext = std::unique_ptr<char[]>( new char[__l_text] );
	    bgzf_read(in,__htext.get(),__l_text*sizeof(char));
	    bgzf_read(in,&__n_ref,sizeof(I32));
	    for(decltype(__n_ref) i = 0 ; i < __n_ref ; ++i )
	      {
		I32 l_name,l_ref;
		bgzf_read( in,&l_name,sizeof(I32) );
		char name[l_name];
		bgzf_read( in,&name[0],l_name*sizeof(char));
		bgzf_read( in,&l_ref,sizeof(I32) );
		__refdata.push_back(std::make_pair(std::string(name),l_ref));
	      }
	  }
      }
  }

  bamreader::bamreader( const char * bamfilename) :
    __impl( new bamreaderImpl(bamfilename) )
  {
  }

  bamreader::~bamreader(){}

  bamrecord bamreader::next_record() const
  {
    return bamrecord(__impl->in);
  }

  bamrecord bamreader::record_at_pos( std::int64_t offset ) const 
  {
    auto current = bgzf_tell(__impl->in);
    int check = bgzf_seek(__impl->in, offset, SEEK_SET );
    if(check==-1) return bamrecord();
    bamrecord b(__impl->in);
    //restore offset
    check = bgzf_seek(__impl->in, current, SEEK_SET );
    if(check == -1) return bamrecord();
    return b;
  }

  bool bamreader::eof() const
  {
    //return gzeof(__impl->in);
  }

  bool bamreader::error() const
  {
    return __impl->__errorstate;
  }

  int bamreader::rewind() 
  {
    return bgzf_seek(__impl->in,0L,SEEK_SET);
    //return gzrewind(__impl->in);
  }

  int bamreader::seek( z_off_t offset, int whence )
  {
    return bgzf_seek(__impl->in,std::move(offset),std::move(whence));
  }

  int bamreader::close()
  {
    return bgzf_close(__impl->in);
  }

  std::int64_t bamreader::tell() 
  {
    return bgzf_tell(__impl->in);
  }

  bamreader::refdataObj bamreader::operator[](const size_type & i) {
    return __impl->__refdata[i];
  }

  bamreader::refdata_citr bamreader::ref_cbegin() const
  {
    return __impl->__refdata.cbegin();
  }

  bamreader::refdata_citr bamreader::ref_cend() const
  {
    return __impl->__refdata.cend();
  }

  std::string bamreader::header() const
  {
    return std::string(__impl->__htext.get());
  }

  std::uint32_t bamreader::n_ref() const {
    return __impl->__n_ref;
  }

}
