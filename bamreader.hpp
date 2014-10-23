#ifndef __SEQUENCE__BAMREADER_HPP__
#define __SEQUENCE__BAMREADER_HPP__

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <bamrecord.hpp>

namespace Sequence 
{

  //fwd declaration
  class bamreaderImpl;

  class bamreader
  {
  public:
    using I32 = std::int32_t;
  private:
    std::unique_ptr<bamreaderImpl> __impl;
  public:
    //! Initialize a bamreader object with a file name
    bamreader( const char * bamfilename = nullptr);
    ~bamreader();

    //! \retun The next alignment in the file
    bamrecord next_record() const;
    //! \return True if bam file is at EOF, false otherwise
    bool eof() const;
    //! \return True if an error was encountered while reading the bam file, false otherwise
    bool error() const;
 
    //! Iterator type (const only!)
    using refdata_citr = std::vector< std::pair<std::string,I32> >::const_iterator;
    //! size_type for the container of reference data
    using size_type = std::vector< std::pair<std::string,I32> >::size_type;
    //! Typedef for reference data (sequence name, length)
    using refdataObj = std::pair<std::string,I32>;
    //! \return Const iterator pointing to info for first reference sequence
    refdata_citr ref_cbegin() const;
    //! \return Const iterator to end of reference data
    refdata_citr ref_cend() const;
    //! \return the i-th reference name/length pair
    refdataObj operator[](const size_type & i);
    //! \return The complete header
    std::string header() const;
    //! \return The number of sequences in the reference
    std::uint32_t n_ref() const;
  };

}

#endif
