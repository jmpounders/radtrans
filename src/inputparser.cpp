#include "inputparser.h"
#include "log.h"

InputParser::InputParser()
{
}

InputParser::~InputParser()
{
  // Need to delete data sets
}

bool
InputParser::parseInputFile(std::string fileName)
{
  _readInputFile(fileName);
  _parseInput();
  return _check();
}

void
InputParser::echoInput()
{
 // Print the result
  DataSetPrinter printer;
  printer(inputData);
}

bool
InputParser::_parseInput()
{
  // Parse file into data set structure using grammar
  std::string::const_iterator iter = fileContents.begin();
  std::string::const_iterator end = fileContents.end();
  bool parseResult = qi::phrase_parse(iter, end, inputGrammar, ascii::space, inputData);
  return parseResult;
}


bool
InputParser::_check()
{
  std::vector<std::string> path(1,"Transport");
  std::vector<double> v;
  std::string s;

  bool abort = false;
  int dataSetCounter = 0;
  int solnManagers = 0;

  path[0] = "Mesh";
  dataSetCounter = countSets(path);
  if (dataSetCounter == 0) {
    LOG_ERR("Must have at least one Mesh data set.");
    abort = true;
  }
 
  path[0] = "Material";
  dataSetCounter = countSets(path);
  if (dataSetCounter == 0) {
    LOG_ERR("Must have at least one Material data set.");
    abort = true;
  }

  path[0] = "Transport";
  dataSetCounter = countSets(path);
  if (dataSetCounter != 1) {
    LOG_ERR("Must have one and only one Transport data set.");
    abort = true;
  }

  path[0] = "Solver";
  dataSetCounter = countSets(path);
  if (dataSetCounter != 1) {
    LOG_ERR("Must have one and only one Solver data set.");
    abort = true;
  }
  
  s = getString(path, "type");
  if (s == "empty") {
    LOG_ERR("Must specify a solver type; use the type parameter.");
    abort = true;
  }

  
  path[0] = "FixedSource";
  dataSetCounter = countSets(path);
  solnManagers += dataSetCounter;
  if (dataSetCounter > 1) {
    LOG_ERR("Found more than one FixedSource solvers.  This is not allowed.");
    abort = true;
  }

  path[0] = "Eigenvalue";
  dataSetCounter = countSets(path);
  solnManagers += dataSetCounter;
  if (dataSetCounter > 1) {
    LOG_ERR("Found more than one Eigenvalue solvers.  This is not allowed.");
    abort = true;
  }

  path[0] = "Transient";
  dataSetCounter = countSets(path);
  solnManagers += dataSetCounter;
  if (dataSetCounter > 1) {
    LOG_ERR("Found more than one Transient solvers.  This is not allowed.");
    abort = true;
  }

  if (solnManagers != 1) {
    LOG_ERR("Must have one and only solution manager; numer found = ", solnManagers);
    abort = true;
  }

  return !abort;
}

bool
InputParser::_fillIn()
{
  std::vector<std::string> path(1,"Transport");
  std::vector<double> v;
  std::string s;

  // Set quadOrder = 0 if it is not present.  This may be the case with a diffusion solver.
  v = getVector(path, "quadOrder");
  if (v.size() == 0) {
    v.push_back(0);
    setVector(path, "quadOrder", v);
  }

  return true;
}

bool
InputParser::_readInputFile(std::string fileName)
{
  // Open input file and don't skip white space
  std::ifstream inputFile;
  inputFile.open(fileName, std::ios_base::in);
  inputFile.unsetf(std::ios::skipws);

  // Read entire file into memory
  std::copy(std::istream_iterator<char>(inputFile),
            std::istream_iterator<char>(),
            std::back_inserter(fileContents));
  inputFile.close();

  return true;
}


