#ifndef LOG_H
#define LOG_H

#include "logger.h"

extern Logging::Logger<Logging::FileLogPolicy> logInstance;

#define LOG      logInstance.print<Logging::info>
#define LOG_DBG  logInstance.print<Logging::debug>
#define LOG_ERR  logInstance.print<Logging::error>
#define LOG_WARN logInstance.print<Logging::warning>

#endif
