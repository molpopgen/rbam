#ifndef __SEQUENCE__BAMRECORD_HPP__
#define __SEQUENCE__BAMRECORD_HPP__

#include <cstdint>
#include <memory>
#include <string>
#include <zlib.h>
#include <Sequence/samrecord.hpp>

namespace Sequence
{

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
    //! True if a record could not be read, or some sort of error occurred
    bool empty() const;

    //Member functions relating to stuff that users care about

    //! \return The read name
    std::string read_name() const;
    //! \return The sequence in base space
    std::string seq() const;
    //! \return Unpacked cigar string
    std::string cigar() const;
    /*!
      \return Quality score string.
      \note Same scale as what is stored in the BAM file
    */
    std::string qual() const;
    //! \return a Sequence::samflag
    samflag flag() const;
    //! \return The mapping quality score
    std::uint32_t mapq() const;
    //! \return The mapping position of the read.  O-offset. -1 = unmapped
    std::int32_t pos() const;
    //! \return The mapping position of the read's mate.  0-offset. -1 = mate unmapped
    std::int32_t next_pos() const;
    //! \return The ID number of the reference sequence where this read maps.  0-offset, -1 = unmapped
    std::int32_t refid() const;
    //! \return The ID number of the reference sequence where this read's mate maps.  0-offset, -1 = mate unmapped
    std::int32_t next_refid() const;
    //! \return Template length
    std::int32_t tlen() const;
    //! \return The length of the read
    std::int32_t l_seq() const;
    //! \return All auxillary data as a single string
    std::string aux() const;
    //Iterator member functions
    //! Beginning of encoded sequence data. Each integer contains 2 bases, with first base stored in the "high nibble"
    const std::uint8_t * seq_cbeg() const;
    //! One past end of encoded sequence data.
    const std::uint8_t * seq_cend() const;
    //! Beginning of quality data
    const char * qual_cbeg() const;
    //! One past end of quality data
    const char * qual_cend() const;

    /*! Search auxillary data for a specific tag
      \return The first position of the match if it exists, nullptr if it does not
    */
    const char * hasTag(const char * tag);
    /*! Search auxillary data for a specific tag, beginning at position start
      \return The first position of the match if it exists, nullptr if it does not
    */
    const char * hasTag(const char * start, const char * tag);
  };

}

#endif
