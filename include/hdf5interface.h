#ifndef HDF5INTERFACE_H
#define HDF5INTERFACE_H

#include <string>
#include <map>
#include <vector>

#include "H5Cpp.h"

class HDF5StreamableArray;

///HDF data storage container
/**
 *  This struct holds HDF data of any time.  It is interpreted (typed)
 *  by the HDFDataStruct object
 **/
struct HDFData
{
  ~HDFData() {};
  int rank;
  hsize_t dims[3];
  void* data;
};

/// HDF data object
/**
 *  This struct stores the data from a single HDF (floating point) data read and the
 *  associated descriptive parameters.
 **/
template<typename T>
struct HDFDataStruct
{
  ~HDFDataStruct() {if (data) delete data;};
  T* getData() { return (T*)data->data; };
  T operator[](int i) { return *((T*)data->data + i) ; };
  T operator()(int i) { return *((T*)data->data + i) ; };
  T operator()(int i, int j) { return *((T*)data->data + data->dims[1]*i + j) ; };
  void clear()
  {
    if (data->data) delete [] (T*)data->data;
    data->rank = 0;
    data->dims[0] = 0;
    data->dims[1] = 0;
    data->dims[2] = 0;
  };

  HDFData* data;
};

/// Interface class for HDF5 file I/O
/**
 *  This class provides basic functionality for reading and writing HDF5 files.
 **/
class HDF5Interface
{
  friend class HDF5StreamableArray;
 public:
  void open(const std::string& fileName, char mode = 'W');
  void close(const std::string& filename);

  void createDataGroup(const std::string& fileName,
                       const std::string groupName);

  template<class T>
  void writeData(const std::string& fileName,
                 const std::string& varName,
                 T* data,
                 long size)
    {
      // Make sure the file has been opened
      H5::H5File* file = _getFile(fileName);
      if ( file ) {
        // Create the dataspace
        hsize_t dims[] = { hsize_t(size) };
        H5::DataSpace dataspace(1, dims);
        _createAndWrite(file, dataspace, varName, data);
      }
    }
  
  void writeData(const std::string& fileName,
                 const std::string& varName,
                 const std::vector<double>& data);

  HDFData* readData(const std::string& fileName,
                    const std::string& varName);

 private:
  void _createAndWrite(H5::H5File* file, H5::DataSpace& dataspace, const std::string& varName, double* data)
  { H5::DataSet dataset = file->createDataSet(varName, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(data, H5::PredType::NATIVE_DOUBLE); };
  
  void _createAndWrite(H5::H5File* file, H5::DataSpace& dataspace, const std::string& varName, int* data)
  { H5::DataSet dataset = file->createDataSet(varName, H5::PredType::NATIVE_INT, dataspace);
    dataset.write(data, H5::PredType::NATIVE_INT); };

  void _createAndWrite(H5::H5File* file, H5::DataSpace& dataspace, const std::string& varName, long* data)
  { H5::DataSet dataset = file->createDataSet(varName, H5::PredType::NATIVE_ULONG, dataspace);
    dataset.write(data, H5::PredType::NATIVE_INT); };

  H5::H5File* _getFile(const std::string& fileName);

  static std::map<std::string, H5::H5File*> _fileList;         //!< Map of file names to file pointers
                                                               // This is static to prevent the same file
                                                               // from being opened more than once.
};

/// HDF5 streamable arrays object
/**
 *  This class provides a mechanism for streaming data continuously to an HDF5 file.
 *  This is useful when the size of the data is not known beforehand, but is generated
 *  semi-continuously as the solution proceeds.
 **/
class HDF5StreamableArray
{
 public:
  HDF5StreamableArray(std::string& fileName,
                      std::vector<std::string> path,
                      std::string varName,
                      hsize_t chunkSize = 100);
  ~HDF5StreamableArray();
  void operator<<(double datum);
  
 private:
  hsize_t _dims[1];
  hsize_t _maxDims[1];
  hsize_t _chunkDims[1];

  H5::H5File* _file;
  H5::DataSpace* _dataSpace;
  H5::DataSet* _dataSet;
  
  std::vector<std::string> _path;
  std::string _varName;

  double* _data;
  int _dataIndex;
  int _numChunksWritten;
  hsize_t _size[1];
  hsize_t _offset[1];
};

#endif
