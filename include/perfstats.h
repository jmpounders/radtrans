#ifndef PERFSTATS_H
#define PERFSTATS_H

#include <string>
#include <map>

#include "timing.h"

/// Collect performance statistics
/**
 *  Collect wall and cpu time measurements for the life of this
 *  object.  
 **/
class PerfStats
{
 public:
  PerfStats(std::string const obj)
  {
    objectName = obj;
    _wallTimeStart = Timing::get_wall_time();
    _cpuTimeStart = Timing::get_cpu_time();
  };

  ~PerfStats()
  {
    wallTime[objectName] += Timing::get_wall_time() - _wallTimeStart;
    cpuTime[objectName] += Timing::get_cpu_time() - _cpuTimeStart;
    ++count[objectName];
  };

  static void display();

  static std::map<std::string, double> wallTime;
  static std::map<std::string, double> cpuTime;
  static std::map<std::string, int> count;
  
 private:
  std::string objectName;
  
  double _wallTimeStart;
  double _cpuTimeStart;

  static void _print(const std::string& obj, const double& time, const std::string& units);
};

#endif
