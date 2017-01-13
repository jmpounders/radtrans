#ifndef INPUTPARSER_H
#define INPUTPARSER_H

#include <vector>
#include <string>
#include <fstream>

#include "boost.h"
#include "dataset.h"


// Input parser grammar
/*
    The following simple grammar defines one rule for reading named arrays, i.e.
    the data produced from this rule (the synthesized attribute) will be of type
    namedArrayData_t.
    Note that these types of grammars are composable.
 */

/// Grammar for reading named double arrays
template <typename Iterator>
struct NamedArrayGrammar : qi::grammar<Iterator, namedArrayData_t(), ascii::space_type >
{
  NamedArrayGrammar() : NamedArrayGrammar::base_type(namedArray)
  {
    text  %= qi::lexeme[+(qi::char_ - " ")];
    doubleArray %= qi::double_ % ',';
    namedArray %= text >> doubleArray;
  }

  qi::rule<Iterator, std::string(), ascii::space_type> text;
  qi::rule<Iterator, std::vector<double>(), ascii::space_type> doubleArray;
  qi::rule<Iterator, namedArrayData_t(), ascii::space_type> namedArray;
};

/// Grammar for reading named string parameters
template <typename Iterator>
struct NamedStringGrammar : qi::grammar<Iterator, namedStringData_t(), ascii::space_type >
{
  NamedStringGrammar() : NamedStringGrammar::base_type(namedString)
  {
    text  %= qi::lexeme[+(qi::char_ - ' ' - '<' - qi::eol)];
    namedString %= text >> text;
  }

  qi::rule<Iterator, std::string(), ascii::space_type> text;
  qi::rule<Iterator, namedStringData_t(), ascii::space_type> namedString;
};

/// Grammar for reading a data set consisting of zero or more data nodes
/**
 *  Each data node is either a named array, named string, or another data set.
 **/
template <typename Iterator>
struct DataSetGrammar : qi::grammar<Iterator, DataSetStruct(), qi::locals<std::string>, ascii::space_type >
{
  DataSetGrammar() : DataSetGrammar::base_type(DataSet_r)
    {
      using namespace qi::labels;

      DataNode_r %= DataSet_r | nag | nsg;
      
      startSet %= '<' >> !qi::lit('/')
                >> qi::lexeme[+(qi::char_ - '>')]
                >> '>';

      endSet = "</" >> ascii::string(_r1) >> '>';

      DataSet_r %= startSet[_a = _1]
                >> *DataNode_r
                >> endSet(_a);
    }

  qi::rule<Iterator, DataSetStruct(), qi::locals<std::string>, ascii::space_type> DataSet_r;
  qi::rule<Iterator, DataNode(),        ascii::space_type> DataNode_r;
  qi::rule<Iterator, std::string(),     ascii::space_type> startSet;
  qi::rule<Iterator, void(std::string), ascii::space_type> endSet;
  NamedArrayGrammar<Iterator> nag;
  NamedStringGrammar<Iterator> nsg;
};


/// Input parser object
/**
 *  This class reads in all of the data sets from the input file and stores the result
 *  as a DataSetStruct.
 **/
class InputParser
{
 public:
  InputParser();
  ~InputParser();

  bool parseInputFile(std::string fileName);

  int countSets(std::vector<std::string> path)
    { return setRead.countSets(inputData, path, 0); };
  
  std::vector<double> getVector(std::vector<std::string> path,
                                std::string varName,
                                int setNumber = 1)
    { return setRead.getVector(inputData, path, varName, 0, setNumber); };

  std::string getString(std::vector<std::string> path,
                        std::string varName,
                        int setNumber = 1)
    { return setRead.getString(inputData, path, varName, 0, setNumber); };

  void setVector(std::vector<std::string> path,
                 std::string varName,
                 std::vector<double>& varData,
                 int setNumber = 1)
  { setRead.setVector(inputData, path, varName, varData, 0, setNumber); };

  void echoInput();

 private:
  bool _readInputFile(std::string fileName);
  bool _parseInput();
  bool _check();
  bool _fillIn();

  std::string fileContents;
  DataSetGrammar<std::string::const_iterator> inputGrammar;
  DataSetStruct inputData;
  GetFromSet setRead;

};

#endif
