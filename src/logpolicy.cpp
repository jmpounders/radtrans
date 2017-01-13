#include "logpolicy.h"

namespace Logging
{

  FileLogPolicy::~FileLogPolicy()
  {
    if (outStream)
      closeOStream();
  }

  void
  FileLogPolicy::openOStream(const std::string& name)
  {
    outStream->open( name.c_str() );
    if ( !outStream->is_open() )
      throw(std::runtime_error("Logger: unable to open an output stream."));
  }

  void
  FileLogPolicy::closeOStream()
  {
    if ( outStream )
      outStream->close();
  }

  void
  FileLogPolicy::write(const std::string& msg)
  {
    (*outStream) << msg << std::endl;
    outStream->flush();
  }


}
