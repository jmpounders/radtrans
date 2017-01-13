#include "perfstats.h"

#include <iostream>
#include <iomanip>

std::map<std::string, double> PerfStats::wallTime;
std::map<std::string, double> PerfStats::cpuTime;
std::map<std::string, int> PerfStats::count;


void
PerfStats::display()
{
  std::cout << std::endl;
  std::cout << "Performance Statistics" << std::endl;
  std::cout << "================================" << std::endl;
  std::cout << "Wall times:" << std::endl ;
  for (std::map<std::string, double>::iterator it = wallTime.begin(); it != wallTime.end(); ++it)
    _print(it->first, it->second, "seconds");
  std::cout << std::endl;
  std::cout << "CPU times:" << std::endl;
  for (std::map<std::string, double>::iterator it = cpuTime.begin(); it != cpuTime.end(); ++it)
    _print(it->first, it->second, "seconds");
  std::cout << std::endl;
  std::cout << "Call counts:" << std::endl;
  for (std::map<std::string, int>::iterator it = count.begin(); it != count.end(); ++it)
    _print(it->first, it->second, "calls");
  std::cout << "================================" << std::endl;
  std::cout << std::endl;
}

void
PerfStats::_print(const std::string& obj, const double& time, const std::string& units)
{
  if (obj.size() < 40)
    std::cout << "  " << obj << std::string(40-obj.size(), '.') << ": " << time << " " << units << std::endl;
  else
    std::cout << "  " << obj  << ": " << time << " " << units << std::endl;
}
