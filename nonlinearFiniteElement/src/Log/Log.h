#ifndef LOG_H
#define LOG_H

#include <iostream>
#include <string>
#include <vector>


#if _DEBUG == 1
#define ARRAY(x,y) Log::Logger().Info(x,y)
#define INFO(x) Log::Logger().Info(x)
#define WARNING(x) Log::Logger().Warning(x)
#elif _DEBUG == 2
#define ARRAY(x,y) Log::Logger().Info(x,y)
#define INFO(x)
#define WARNING(x)
#elif _RELEASE == 1
#define ARRAY(x,y)
#define INFO(x)
#define WARNING(x)
#endif

class Log {
public:
  static const unsigned int LevelError = 0;
  static const unsigned int LevelWarning = 1;
  static const unsigned int LevelInfo = 2;
private:
  static unsigned int m_LogLevel;

public:
  static Log& Logger()
  {
    static Log instanceOfLog;
    return instanceOfLog;
  };
  Log();
  void setLevel(unsigned int);
  void Error(const std::string&);
  void Warning(const std::string&);
  void Info(const std::string&);
  void Info(const double);
  void Info(const unsigned int *, unsigned int);
  void Info(const double *, unsigned int);
};
#endif
