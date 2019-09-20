#include "Log.h"

Log::Log()
{
};

unsigned int Log::m_LogLevel = LevelInfo;

void Log::setLevel(unsigned int level)
{
  m_LogLevel = level;
};

void Log::Error(const std::string& message)
{
  if (m_LogLevel >= LevelError)
    std::cout << "\033[1;31m[ERROR]: \033[0m" << message << std::endl;
};

void Log::Warning(const std::string& message)
{
  if (m_LogLevel >= LevelWarning)
    std::cout << "\033[1;33m[WARNING]: \033[0m" << message << std::endl;
};

void Log::Info(const std::string& message)
{
  if (m_LogLevel >= LevelInfo)
    std::cout << "\033[1;32m[Info]: \033[0m" << message << std::endl;
};

void Log::Info(const double d)
{
  if (m_LogLevel >= LevelInfo)
    std::cout << "\033[1;32m[Info]: \033[0m" << std::to_string(d) << std::endl;
};


void Log::Info(const unsigned int *array, unsigned int printSize)
{
  if (m_LogLevel >= LevelInfo)                                                  
    {                                                                           
      std::cout << "\033[1;34m[Array]: \033[0m" << std::endl;                   
      for (unsigned int i = 0; i < printSize; i++)                                      
	std::cout << "a[" << i << "]= " << array[i] << std::endl;        
    }                                                                           
      
}

void Log::Info(const double *array, unsigned int printSize)
{
  if (m_LogLevel >= LevelInfo)
    {
      std::cout << "\033[1;34m[Array]: \033[0m" << std::endl;
      for (unsigned int i = 0; i < printSize; i++)
        std::cout << "a[" << i << "]= " << array[i] << std::endl;
    }

}
/*
void Log::Info(const double *array, unsigned int printSize)
{
  if (m_LogLevel >= LevelInfo)
    {
      std::cout << "\033[1;34m[Array]: \033[0m" << std::endl;
      for (unsigned int i = 0; i < printSize; i++)
        std::cout << "a[" << i << "]= " << array[i] << std::endl;
    }

}
*/
