#include "outputgenerator.h"

#include "boost/spirit/include/karma_string.hpp"

  /*
  // Test OutputParser
  double d1 = 1.234;
  vector<double> v1 = {1.1,2.2,3.3};
  vector<double> v2 = {10,9,8,7,6,5,4,3,2,1};

  OutputGenerator outputer;
  string generated;
  back_insert_iterator<string> sink(generated);
  std::string name("v1");
  bool generateResult = outputer.generateNamedArray(sink,
						    name,
						    *data.getData<vector<double> >(name));
  if (generateResult)
    {
      cout << "Generate success " << endl;
      cout << generated << endl;
    }
  else
    {
      cout << "Generate failed" << endl;
    }
  data.printDataSet(outputer);
  // End OutputParser test
  */

namespace karma = boost::spirit::karma;
namespace ascii = boost::spirit::ascii;

OutputGenerator::OutputGenerator()
{
}

OutputGenerator::~OutputGenerator()
{
  // Need to delete data sets
}


/*
 *  Generate an input string of the form
 *  "name 1,2,3,..."
 */
bool
OutputGenerator::generateNamedArray(std::back_insert_iterator<std::string>& sink,
                                    std::string& name,
                                    std::vector<double>& values)
{
  return karma::generate(sink,
                         ascii::string << karma::lit(" ") << (karma::double_ % ','),
                         name,values);
}
