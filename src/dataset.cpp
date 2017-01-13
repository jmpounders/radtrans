#include "dataset.h"
#include <iostream>

int
GetFromSet::countSets(DataSetStruct const& data,
                      std::vector<std::string> const& setName,
                      int level) const
{
  int counter = 0;
  GetNodeNameVisitor getName;
  std::string nodeName;

  if ( setName.size() == level+1 ) {
    // Count the number of data sets with this name
    BOOST_FOREACH(DataNode const& node, data.children)
      {
        nodeName = boost::apply_visitor(getName, node);
        if ( nodeName == setName[level] ) {
          counter++;
        }
      }
  }
  else {
    // Recurse
    BOOST_FOREACH(DataNode const& node, data.children)
      {
        nodeName = boost::apply_visitor(getName, node);

        if ( nodeName == setName[level] ) {
          return countSets(boost::get<DataSetStruct>( node ), setName, level+1);
        }
      }
  }
  return counter;
}

std::vector<double>
GetFromSet::getVector(DataSetStruct const& data,
                std::vector<std::string> const& setName,
                std::string const& varName,
                int level,
                int setNumber) const
{
  GetNodeNameVisitor getName;
  GetNodeValueVisitor getVal;
  std::string nodeName;
  int setCount = 0;

  BOOST_FOREACH(DataNode const& node, data.children)
    {
      nodeName = boost::apply_visitor(getName, node);

      if (level < setName.size()) {
        if ( nodeName == setName[level] ) {
          setCount++;
          if ( setNumber == setCount) {
            return getVector(boost::get<DataSetStruct>( node ), setName, varName, level+1);
          }
        }
      }
      else if ( nodeName == varName ) {
        return boost::apply_visitor(getVal, node);;
      }
    }
  return getVal.empty;
}

std::string
GetFromSet::getString(DataSetStruct const& data,                      
                std::vector<std::string> const& setName,
                std::string const& varName,
                int level,
                int setNumber) const
{
  GetNodeNameVisitor getName;
  GetNodeStringVisitor getStr;
  std::string nodeName;
  int setCount = 0;

  BOOST_FOREACH(DataNode const& node, data.children)
    {
      nodeName = boost::apply_visitor(getName, node);
      
      if (level < setName.size()) {
        if ( nodeName == setName[level] ) {
          setCount++;
          if ( setNumber == setCount )
            return getString(boost::get<DataSetStruct>( node ), setName, varName, level+1);
        }
      }
      else if ( nodeName == varName ) {
        return boost::apply_visitor(getStr, node);;
      }
    }
  return getStr.empty;
}


void
GetFromSet::setVector(DataSetStruct& data,
                      std::vector<std::string> const& setName,
                      std::string const& varName,
                      std::vector<double>& varData,
                      int level,
                      int setNumber)
{
  GetNodeNameVisitor getName;
  GetNodeValueVisitor getVal;
  std::string nodeName;
  int setCount = 0;

  std::cerr << level << " " << setName.size() << std::endl;

  if (level == setName.size()) {
    std::cerr << "adding data to " << data.setName << std::endl;
    namedArrayData_t namedArray(varName, varData);
    data.children.push_back( DataNode(namedArray) );
  }
  else {
    BOOST_FOREACH(DataNode& node, data.children)  {
      nodeName = boost::apply_visitor(getName, node);

      if ( nodeName == setName[level] ) {
        setCount++;
        if ( setNumber == setCount)
          setVector(boost::get<DataSetStruct>( node ), setName, varName, varData, level+1);
      }
    }
  }
}



// Printer
void
tab(int indent)
{
  for (int i=0; i<indent; ++i) {
    std::cout << ' ';
  }
}

// Following are two overloads of a node printer
void
DataNodePrinter::operator()(DataSetStruct const& data) const
{
  DataSetPrinter(indent+tabsize)(data);
}

void
DataNodePrinter::operator()(namedStringData_t const& data) const
{
  std::string stringName;
  std::string stringValue;

  stringName = fusion::at_c<0>(data);
  stringValue = fusion::at_c<1>(data);

  tab(indent+tabsize);
  std::cout << stringName << " = " << stringValue;
  std::cout << std::endl;
}

void
DataNodePrinter::operator()(namedArrayData_t const& data) const
{
  std::string arrayName;
  std::vector<double> arrayValues;

  arrayName = fusion::at_c<0>(data);
  arrayValues = fusion::at_c<1>(data);

  tab(indent+tabsize);
  std::cout << arrayName << " = ";
  for (int i=0; i<arrayValues.size(); ++i) {
    std::cout << arrayValues[i] << " ";
  }
  std::cout << std::endl;
}


// Following is the data set printer
void
DataSetPrinter::operator()(DataSetStruct const& data) const
{
  tab(indent);
  std::cout << "Data Set: " << data.setName << std::endl;
  tab(indent);
  std::cout << '{' << std::endl;
  BOOST_FOREACH(DataNode const& node, data.children)
  {
    boost::apply_visitor(DataNodePrinter(indent), node);
  }

  tab(indent);
  std::cout << '}' << std::endl;
}
