#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <chrono>
#include <string>

#include "Log/Log.h"

class Timer {
private:
  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> duration;
  std::string message;
  double ms;
public:
  Timer();
  Timer(const std::string &);
  void Start();
  void Stop();
  ~Timer();
};

#endif
