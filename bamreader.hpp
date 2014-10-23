#ifndef __SEQUENCE__BAMREADER_HPP__
#define __SEQUENCE__BAMREADER_HPP__

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <bamrecord.hpp>

//fwd declaration
class bamreaderImpl;

class bamreader
{
public:
  using I32 = std::int32_t;
private:
  std::unique_ptr<bamreaderImpl> __impl;
public:
  bamreader( const char * bamfilename = nullptr);
  ~bamreader();

  bamrecord next_record() const;
  bool eof() const;
  bool error() const;

  using refdata_citr = std::vector< std::pair<std::string,I32> >::const_iterator;
  refdata_citr ref_cbeg() const;
  refdata_citr ref_cend() const;
};

#endif
