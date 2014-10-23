#include <bamreader.hpp>
#include <zlib.h>

using std::string;

namespace Sequence
{

  //!Impl class for bamreader
  class bamreaderImpl
  {
  public:
    using I32 = std::int32_t;
    gzFile in;
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
    if(in != NULL) gzclose(in);
  }

  bamreaderImpl::bamreaderImpl(const char * bamfilename) :
    in((bamfilename != nullptr) ? gzopen(bamfilename,"rb") : NULL),
    __errorstate(false),
    __htext(nullptr),
    __refdata(std::vector< std::pair<std::string,I32> >())
  {
    if(gzopen != NULL)
      {
	gzread( in, &__magic[0], 4*sizeof(char) );
	if(string({__magic[0],__magic[1],__magic[2]}) != string("BAM")) __errorstate = 1;
	if(!__errorstate)
	  {
	    gzread( in, &__l_text, sizeof(I32) );
	    __htext = std::unique_ptr<char[]>( new char[__l_text] );
	    gzread(in,__htext.get(),__l_text*sizeof(char));
	    gzread(in,&__n_ref,sizeof(I32));
	    for(decltype(__n_ref) i = 0 ; i < __n_ref ; ++i )
	      {
		I32 l_name,l_ref;
		gzread( in,&l_name,sizeof(I32) );
		char name[l_name];
		gzread( in,&name[0],l_name*sizeof(char));
		gzread( in,&l_ref,sizeof(I32) );
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

  bool bamreader::eof() const
  {
    return gzeof(__impl->in);
  }

  bool bamreader::error() const
  {
    return __impl->__errorstate;
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
