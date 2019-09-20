#include "Timer.h"

Timer::Timer()
{
  ms = 0.0;
};

Timer::Timer(const std::string &mes)
  : message(mes)
{
  ms = 0.0;
};


Timer::~Timer()
{
  std::string report = message + std::to_string(ms) + "ms ";
  std::cout << report <<"\n";
};


void Timer::Start()
{
  start = std::chrono::high_resolution_clock::now();
}

void Timer::Stop()
{
  end = std::chrono::high_resolution_clock::now();
  duration = end - start;
  ms = ms + duration.count() * 1000.0f;
}
