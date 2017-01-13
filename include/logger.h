#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <fstream>
#include <sstream>
//#include <mutex>

#include "logpolicy.h"

namespace Logging
{

  /// Logger class
  template< typename LogPolicy >
    class Logger
    {
    public:
      Logger(const std::string& name);
      ~Logger();

      template< SeverityType severity, typename...Args >
        void print( Args...args);
  
    private:
      std::stringstream _logStream;
      LogPolicy* _policy;
      //std::mutex _writeMutex;
  
      void _printImpl();

      template< typename First, typename...Rest >
        void _printImpl(First param1, Rest...param);
    };


  template< typename LogPolicy >
    Logger<LogPolicy>::Logger(const std::string& name)
    {
      _policy = new LogPolicy;
      if ( !_policy )
        throw std::runtime_error("Unable to create logger instance");
      _policy->openOStream(name);
      _policy->write("Radiation Transport Log");
    }

  template< typename LogPolicy >
    Logger<LogPolicy>::~Logger()
    {
      if ( _policy) {
        _policy->closeOStream();
        delete _policy;
      }
    }


  template< typename LogPolicy >
    template< SeverityType severity, typename...Args >
    inline void
    Logger<LogPolicy>::print( Args...args )
  {
    //_writeMutex.lock();
    switch( severity )
      {
      case info :
        _logStream << "INFO    : ";
        break;
      case debug :
        _logStream << "DEBUG   : ";
        break;
      case warning :
        _logStream << "WARNING : ";
        break;
      case error :
        _logStream << "ERROR   : ";
        break;
      }
    _printImpl(args...);
    //_writeMutex.unlock();
  }


  template< typename LogPolicy >
    inline void Logger<LogPolicy>::_printImpl()
    {
      _policy->write( _logStream.str() );
      _logStream.str("");
    }


  template< typename logPolicy >
    template< typename First, typename...Rest >
    inline void Logger< logPolicy >::_printImpl(First param1, Rest...param)
  {
    _logStream << param1;
    _printImpl(param...);
  }

}

#endif
