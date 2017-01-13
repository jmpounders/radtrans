#ifndef DATASET_H
#define DATASET_H

#include <vector>
#include <string>
#include <fstream>

#include "boost.h"

////////////////////////////////////////////////////////////////////////
// Data storage structures
////////////////////////////////////////////////////////////////////////

struct DataSetStruct;

/*
    The following FUSION vector types store heterogenous collections of data,
    in this case, a string and a vector or a string and a string.
*/
typedef fusion::vector<std::string, std::vector<double> > namedArrayData_t;
typedef fusion::vector<std::string, std::string> namedStringData_t;

/// Data node type
/*
 *  This BOOST variant structure represents a node in a data set (tree).  Each
 *  node is either a named array, a named string or another data set.
 */
typedef boost::variant<boost::recursive_wrapper<DataSetStruct>,
                       namedStringData_t,
                       namedArrayData_t> DataNode;

/// Data set structure
/**
 *  A data set consists of a name and a set of nodes.  This struct is
 *  adapted into a BOOST structure so that BOOST methods can be used
 *  to set and retrieve data.
 */
struct DataSetStruct
{
  std::string setName;
  std::vector<DataNode> children;
};
BOOST_FUSION_ADAPT_STRUCT(DataSetStruct,
			  (std::string, setName)
			  (std::vector<DataNode>, children));

////////////////////////////////////////////////////////////////////////
// Data access structures
////////////////////////////////////////////////////////////////////////

/// Structure to get values from a data set.
/**
 *  This struct provides functions to get either a vector of doubles
 *  or a string from a DataSetStruct. The implementation uses a set
 *  a set of pre-defined visitors that are applied to the data
 *  hierarchy using boost.
 */
struct GetFromSet
{
  GetFromSet() {};
  int countSets(DataSetStruct const& data,
                std::vector<std::string> const& setName,
                int level) const;
  std::vector<double> getVector(DataSetStruct const& data,
                                std::vector<std::string> const& setName,
                                std::string const& varName,
                                int level,
                                int setNumber = 1) const;
  std::string getString(DataSetStruct const& data,
                        std::vector<std::string> const& setName,
                        std::string const& varName,
                        int level,
                        int setNumber = 1) const;
  void setVector(DataSetStruct& data,
                 std::vector<std::string> const& setName,
                 std::string const& varName,
                 std::vector<double>& varData,
                 int level,
                 int setNumber = 1);
};

/// Visitor to get the name of a node within a data set.
struct GetNodeNameVisitor : boost::static_visitor<std::string>
{
  std::string const& operator()(DataSetStruct const& data) const { return data.setName; };
  std::string const& operator()(namedStringData_t const& data) const { return fusion::at_c<0>(data); };
  std::string const& operator()(namedArrayData_t const& data) const { return fusion::at_c<0>(data); };
};

/// Visitor to get a vector from a node within a data set
struct GetNodeValueVisitor : boost::static_visitor<std::vector<double> >
{
  std::vector<double> const& operator()(DataSetStruct const& data) const { return empty; };
  std::vector<double> const& operator()(namedStringData_t const& data) const { return empty; };
  std::vector<double> const& operator()(namedArrayData_t const& data) const { return fusion::at_c<1>(data); };
  std::vector<double> empty;
};

/// Visitor to get a string from a node within a data set
struct GetNodeStringVisitor : boost::static_visitor<std::string>
{
 GetNodeStringVisitor() : doubleVector("double vector"), dataStruct("data set"), empty("empty") {};
  std::string const& operator()(DataSetStruct const& data) const { return dataStruct; };
  std::string const& operator()(namedStringData_t const& data) const { return fusion::at_c<1>(data);; };
  std::string const& operator()(namedArrayData_t const& data) const { return doubleVector; };
  std::string doubleVector;
  std::string dataStruct;
  std::string empty;
};


////////////////////////////////////////////////////////////////////////
// Data printers
////////////////////////////////////////////////////////////////////////

int const tabsize = 2;

/// Structure for printing the contents of a data set
/**
 *  This structure provides a wrapper for recusively looping through
 *  the nodes of a data set and calling visitor functions to print
 *  the contents.
 */
struct DataSetPrinter
{
  DataSetPrinter(int indent=0) : indent(indent) {}
  void operator()(DataSetStruct const& data) const;
  int indent;
};

/// Visitor that will print the contents of a data node
struct DataNodePrinter : boost::static_visitor<>
{
  DataNodePrinter(int indent = 0) : indent(indent) {}
  void operator()(DataSetStruct const& data) const;
  void operator()(namedStringData_t const& data) const;
  void operator()(namedArrayData_t const& data) const;

  int indent;
};



#endif
