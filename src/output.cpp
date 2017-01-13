#include "output.h"

#include<sstream>
#include<vector>

Output::Output(std::string fileName, int format) :
  outputFormat(format), outputFileName(fileName), solutionPoint(0)
{
  if ( format & ASCII )
    outputFile.open(outputFileName);

  if ( format & HDF5 ) {
    hdf.open(outputFileName);
    _createStructure();
  }
}

Output::~Output()
{
  if (outputFormat & ASCII)
    outputFile.close();

  if (outputFormat & HDF5)
    hdf.close(outputFileName);

  for (std::map<std::string, HDF5StreamableArray*>::iterator it=_listStreamables.begin();
       it!=_listStreamables.end();
       ++it) {
    delete it->second;
  }

}


void
Output::registerStreamingOutput(std::string varName)
{
  if (outputFormat & HDF5) {
    std::vector<std::string> path;
    _listStreamables[varName] = new HDF5StreamableArray(outputFileName, path, varName);
  }
}


void
Output::streamOutput(std::string varName, double datum)
{
  if (outputFormat & HDF5) {
    *_listStreamables[varName] << datum;
  }
}

std::string
Output::getName()
{
  return outputFileName;
}


void
Output::_createStructure()
{
  hdf.createDataGroup(outputFileName, "/Mesh");
}


