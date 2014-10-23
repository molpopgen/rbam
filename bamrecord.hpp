#ifndef __SEQUENCE__BAMRECORD_HPP__
#define __SEQUENCE__BAMRECORD_HPP__

#include <cstdint>
#include <memory>
#include <string>
#include <zlib.h>
#include <Sequence/samrecord.hpp>

//fwd declaration
class bamrecordImpl;
class bamrecord 
{
private:
  std::unique_ptr<bamrecordImpl> __impl;
public:
  using U32 = std::uint32_t;
  using I32 = std::int32_t;
  using U8 = std::uint8_t;
  //constructors
  bamrecord( gzFile in );
  bamrecord( const bamrecord & );
  //move constructors
  bamrecord(bamrecord&&);
  //destructor
  ~bamrecord();


  bamrecord & operator=(bamrecord&&);
  bamrecord & operator=(const bamrecord&);

  //member functions
  bool empty() const;

  //Member functions relating to stuff that users care about
  std::string read_name() const;
  std::string seq() const;
  std::string cigar() const;
  std::string qual() const;
  Sequence::samflag flag() const;
  std::uint32_t mapq() const;
  std::int32_t pos() const;
  std::int32_t next_pos() const;
  std::int32_t refid() const;
  std::int32_t next_refid() const;
  std::int32_t tlen() const;
  std::int32_t l_seq() const;
  std::string aux() const;
  //Iterator member functions
  //! Beginning of encoded sequence data
  const std::uint8_t * seq_cbeg() const;
  //! One past end of encoded sequence data
  const std::uint8_t * seq_cend() const;
  //! Beginning of quality data
  const char * qual_cbeg() const;
  //! One past end of quality data
  const char * qual_cend() const;

  const char * hasTag(const char * tag);
  const char * hasTag(const char * start, const char * tag);
};

#endif
