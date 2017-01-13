#include "timing.h"

#include <iostream>

namespace Timing
{
  Timers timer;

  
  Timers::~Timers()
  {
    for (std::map<std::string, Stopwatch*>::iterator it=_timerList.begin(); it!=_timerList.end(); ++it) {
      delete it->second;
    }
  }

  
  void
  Timers::display()
  {
    std::cout << std::endl << "Timing Results" << std::endl;
    Stopwatch* sw;
    for (std::map<std::string, Stopwatch*>::iterator it=_timerList.begin(); it!=_timerList.end(); ++it) {
      sw = it->second;
      std::cout << "  " << it->first << " = " << sw->getWallTime() << std::endl;
    }
  }

  
}
