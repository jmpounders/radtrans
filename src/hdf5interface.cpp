#include "hdf5interface.h"

#include <iostream>

// Instantiation of static file list
std::map<std::string, H5::H5File*> HDF5Interface::_fileList;

void
HDF5Interface::open(const std::string& fileName, char mode)
{
  // Create new file
  H5::H5File* file = _getFile(fileName);
  if ( !file )
    {
      // Note these access modes are not checked later during reads/writes
      if (mode == 'W') {
        file = new H5::H5File( fileName+".h5", H5F_ACC_TRUNC );
      }
      else if (mode == 'R') {
        file = new H5::H5File( fileName+".h5", H5F_ACC_RDONLY );
      }
      _fileList[fileName] = file;
    }
  else
    {
      std::cerr << "File has already been opened\n";
    }
}

void
HDF5Interface::close(const std::string& fileName)
{
  std::map<std::string, H5::H5File*>::iterator _fileIt = _fileList.find( fileName );
  if ( _fileIt == _fileList.end() )
    {
      std::cerr << "File does not exist in memory\n";
    }
  else
    {
      delete _fileIt->second;
      _fileList.erase( _fileIt );
    }
}

void
HDF5Interface::createDataGroup(const std::string& fileName, const std::string groupName)
{
  H5::H5File* file = _getFile(fileName);
  if ( !file )
    {
      std::cerr << "File has not been opened\n";
      return;
    }
  H5::Group* newGroup = new H5::Group( file->createGroup( groupName  ) );
}

/**
 *  Write a 1D vector of doubles to file
 */
void
HDF5Interface::writeData(const std::string& fileName,
                         const std::string& varName,
                         const std::vector<double>& data)
{
  // Store the vector as contiguous memory
  double* _data = new double [data.size()];
  for (int i=0; i<data.size(); i++) {
    _data[i] = data[i];
  }

  writeData(fileName, varName, _data, data.size());
  
  delete [] _data;

  // For the following the data set must already exist
  //H5::DataSet dataset = file->openDataSet(dataSetName);
}


HDFData*
HDF5Interface::readData(const std::string& fileName,
                        const std::string& varName)
{
  // Make sure the file has been opened
  H5::H5File* file = _getFile(fileName);
  if ( !file ) {
    std::cerr << "File does not exist in memory\n";
    // Should throw error and return here...
  }

  H5::DataSet dataset = file->openDataSet(varName);
  H5::DataSpace dataspace = dataset.getSpace();

  HDFData* dataStruct = new HDFData;

  // Get the rank and dimensions of the data
  dataStruct->rank = dataspace.getSimpleExtentNdims();
  int ndims = dataspace.getSimpleExtentDims(dataStruct->dims, NULL);

  // Create the memory space for the data
  H5::DataSpace memspace(dataStruct->rank, dataStruct->dims);
  int totalLength = 1.0;
  for (int i=0; i<dataStruct->rank; i++)
    totalLength = totalLength*dataStruct->dims[i];

  // Determine the data type and read the data
  H5T_class_t dataType = dataset.getTypeClass();
  if (dataType == H5T_INTEGER) {
    dataStruct->data = new int [totalLength];
    dataset.read(dataStruct->data, H5::PredType::NATIVE_INT, memspace, dataspace);
  }
  else if (dataType == H5T_FLOAT) {
    dataStruct->data = new double [totalLength];
    dataset.read(dataStruct->data, H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
  }
  else {
    std::cout << "data type is unknown\n" << dataType << std::endl;
  }

  return dataStruct;
}


H5::H5File*
HDF5Interface::_getFile(const std::string& fileName)
{
  std::map<std::string, H5::H5File*>::iterator _fileIt;
  _fileIt = _fileList.find( fileName );
  if ( _fileIt == _fileList.end() )
    {
      return NULL;
    }
  else
    {
      return _fileIt->second;
    }

}


HDF5StreamableArray::HDF5StreamableArray(std::string& fileName,
                                         std::vector<std::string> path,
                                         std::string varName,
                                         hsize_t chunkSize) :
  _path(path), _varName(varName)
{
  _dims[0] = 0;
  _maxDims[0] = H5S_UNLIMITED;
  HDF5Interface hdf;
  _file = hdf._getFile(fileName);
  
  _chunkDims[0] = chunkSize;

  _data = new double [_chunkDims[0]];
  _dataIndex = 0;
  _numChunksWritten = 0;
  _size[0] = 0;
  _offset[0] = 0;

  _dataSpace = new H5::DataSpace(1, _dims, _maxDims);

  H5::DSetCreatPropList prop;
  prop.setChunk(1, _chunkDims);

  _dataSet = new H5::DataSet(_file->createDataSet(varName,
                                                  H5::PredType::NATIVE_DOUBLE,
                                                  *_dataSpace,
                                                  prop));
}


HDF5StreamableArray::~HDF5StreamableArray()
{
  if (_dataIndex > 0) {
    // write data
    _size[0] += _dataIndex;
    _dataSet->extend(_size);

    H5::DataSpace *fileSpace = new H5::DataSpace(_dataSet->getSpace());
    _offset[0] = _numChunksWritten*_chunkDims[0];
    
    _chunkDims[0] = _dataIndex;
    fileSpace->selectHyperslab(H5S_SELECT_SET, _chunkDims, _offset);

    H5::DataSpace *memSpace = new H5::DataSpace(1, _chunkDims, NULL);

    _dataSet->write(_data, H5::PredType::NATIVE_DOUBLE, *memSpace, *fileSpace);

    _numChunksWritten++;
    _dataIndex = 0;
  }
}

void
HDF5StreamableArray::operator<<(double datum)
{
  _data[_dataIndex] = datum;
  _dataIndex++;
  if (_dataIndex == _chunkDims[0]) {
    // write data
    _size[0] += _chunkDims[0];
    _dataSet->extend(_size);

    H5::DataSpace *fileSpace = new H5::DataSpace(_dataSet->getSpace());
    _offset[0] = _numChunksWritten*_chunkDims[0];
    fileSpace->selectHyperslab(H5S_SELECT_SET, _chunkDims, _offset);

    H5::DataSpace *memSpace = new H5::DataSpace(1, _chunkDims, NULL);

    _dataSet->write(_data, H5::PredType::NATIVE_DOUBLE, *memSpace, *fileSpace);

    _numChunksWritten++;
    _dataIndex = 0;
  }
}
