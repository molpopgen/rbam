#ifndef __SEQUENCE__BAMUTIL_HPP__
#define __SEQUENCE__BAMUTIL_HPP__

namespace bamutil {
  static const char int2seq[16] = {'=','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'};
  static const char bamCig[9] = {'M','I','D','N','S','H','P','=','X'};
  //! What stored in __rdata1 of a bamrecordImpl
  enum class I32s : size_t {REFID,POS,L_SEQ,NEXT_REFID,NEXT_POS,TLEN};
  static const size_t lI32s = 6;
  //! What stored in __rdata2 of a bamrecordImpl
  enum class U32s : size_t {BIN_MQ_NL,FLAG_NC,BIN,MAPQ,L_READ_NAME,FLAG,NC};
  static const size_t lU32s = 7;
}

#endif
