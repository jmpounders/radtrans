#ifndef OUTPUTGENERATOR_H
#define OUTPUTGENERATOR_H

#include <vector>
#include <string>

#include "boost/config/warning_disable.hpp"
#include "boost/spirit/include/karma.hpp"

class OutputGenerator
{
 public:
  OutputGenerator();
  ~OutputGenerator();

  bool generateNamedArray(std::back_insert_iterator<std::string>& sink,
                           std::string& name,
                           std::vector<double>& values);

};

#endif
