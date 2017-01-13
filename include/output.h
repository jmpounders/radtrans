#ifndef OUTPUT_H
#define OUTPUT_H

#include<iostream>
#include<fstream>
#include<string>

#include "hdf5interface.h"

/// Output object, currently supporting HDF5 and ASCII output
class Output
{
 public:
  enum OutputFormat { ASCII=1, HDF5=2, MESH=4 };
  int outputFormat;

  Output(std::string fileName, int format = ASCII);
  ~Output();

  template<class T>
  void writeData(std::string varName, T* data, int size)  {
    if (outputFormat & HDF5)
      hdf.writeData(outputFileName, varName, data, size);
  };

  void registerStreamingOutput(std::string varName);
  void streamOutput(std::string varName, double datum);

  std::string getName();

  
 private:
  void _createStructure();
  
  std::string outputFileName;
  std::ofstream outputFile;
  HDF5Interface hdf;

  int solutionPoint;

  std::map<std::string, HDF5StreamableArray*> _listStreamables;

};

#endif
