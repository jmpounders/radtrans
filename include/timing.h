#ifndef TIMING_H
#define TIMING_H

#include <sys/time.h>
#include <map>
#include <string>

namespace Timing
{

  static double get_wall_time() {
    struct timeval time;
    if (gettimeofday( &time, NULL )) return 0;
    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
  }

  static double get_cpu_time() {
    return (double)clock() / CLOCKS_PER_SEC;
  }

  /// Captures elapsed wall and CPU times
  class Stopwatch
  {
   public:
    Stopwatch() : clockRunning(false), totalTime_wall(0), totalTime_cpu(0) {};
    double punch()
    {
      double currentTime_wall = get_wall_time();
      double currentTime_cpu  = get_cpu_time();
      if (clockRunning) {
        totalTime_wall += currentTime_wall - lastPunch_wall;
        totalTime_cpu  += currentTime_cpu  - lastPunch_cpu;
      }
      lastPunch_wall = currentTime_wall;
      lastPunch_cpu  = currentTime_cpu;
      clockRunning = !clockRunning;
      return totalTime_wall;
    };

    double getWallTime() { return totalTime_wall; };
    double getCPUTime()  { return totalTime_cpu;  };

   private:
    bool clockRunning;

    double totalTime_wall;
    double lastPunch_wall;

    double totalTime_cpu;
    double lastPunch_cpu;
  };

  /// Manges multiple Stopwatch classes
  class Timers
  {
   public:
    ~Timers();
    void addTimer(std::string timerName) { _timerList[timerName] = new Stopwatch; };
    double punch(std::string timerName) { return _timerList[timerName]->punch(); };
    void display();

   private:
    std::map<std::string, Stopwatch*> _timerList;
  };

  /// Global timing class
  /**
   *  This class will be accessible to the entire program.
   **/
  extern Timers timer;

}
 
#endif // TIMING_H
