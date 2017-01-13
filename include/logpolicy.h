#ifndef LOGPOLICY_H
#define LOGPOLICY_H

#include <string>
#include <fstream>
#include <sstream>
#include <memory>
#include <stdexcept>

namespace Logging
{

  /// Generic interface for log policy
  class LogPolicyInterface
  {
  public:
    virtual void openOStream(const std::string& name) = 0;
    virtual void closeOStream() = 0;
    virtual void write(const std::string& msg) = 0;
  };


  /// File logging policy
  class FileLogPolicy : public LogPolicyInterface
  {
  public:
  FileLogPolicy()
    : outStream(new std::ofstream) {};
    ~FileLogPolicy();

    void openOStream(const std::string& name);
    void closeOStream();
    void write(const std::string& msg);
  
  private:
    //std::unique_ptr< std::ofstream > outStream;
    std::ofstream* outStream;
  };


  enum SeverityType { info, debug, error, warning };

}

#endif
