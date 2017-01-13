#ifndef ELEMENT_H
#define ELEMENT_H

#include<vector>

#include "node.h"
#include "material.h"

/// Virtual element class for MeshInterface
class Element
{
 public:
  Element(int nNodes) : _numNodes(nNodes)  {};
  virtual ~Element() {};

  virtual void addNode( Node* newNode ) = 0;
  virtual const Node* getNode( int i ) const = 0;
  virtual double getVolume() = 0;

  void setMat( Material* m ) { _mat = m; };
  Material* getMat() const { return _mat; };

 protected:
  int _numNodes;
  Material* _mat;
};

/// Compact element class (deprecated)
/**
 *  These elments are "compact" in the sense that they do not reference
 *  any memory external to the class.  They are thus more "portable".
 */
class CompactElement : public Element
{
 public:
  CompactElement(int nNodes) : Element(nNodes) { _nodeList.reserve(size_t(nNodes)); };
  ~CompactElement() {};

  void clear() { _nodeList.clear(); _edgeIDList.clear(); _neighborIDList.clear(); };
  void setElementID( long elemID ) { _elementID = elemID; };
  long getElementID() const { return _elementID; };
  void addNode( Node newNode ) { _nodeList.push_back(newNode); };
  void addEdgeID( int edgeID ) { _edgeIDList.push_back(edgeID); };
  long getEdgeID( int edge ) const { return _edgeIDList[edge]; };
  void addNeighbor ( int elemID ) { _neighborIDList.push_back(elemID); };
  long getNeighbor ( int edge ) const { return _neighborIDList[edge]; };

  // Concrete implementations of virtual base class
  void addNode( Node* newNode ) { _nodeList.push_back(*newNode); };
  const Node* getNode( int i ) const { return &_nodeList[i]; };
  double getVolume() { return 0.0; };

 private:
  long _elementID;
  std::vector<Node> _nodeList;
  std::vector<long> _edgeIDList;
  std::vector<long> _neighborIDList;

};


/// Ultra light element class
struct UltraLightElement
{
  long elementID;
  double x[3];
  double y[3];
  double z[3];
  long edgeID[3];
  long neighborID[3];
  long vertexID[3];
};


/// Line (1D) Element class
/**
 *  Currently elements are restricted to one dimension, defined by two points.
 *  These elements store nodes as pointers, so they must be used in a setting
 *  where there is stable external storage of nodes (e.g. the Mesh class).
 */
class LineElement : public Element
{
 public:
  LineElement( Node *n1, Node *n2 );

  // Concrete implementation of virtual base class
  void addNode( Node* newNode ) { _nodeList.push_back(newNode); };
  const Node* getNode( int i ) const { return _nodeList[i]; };
  double getVolume() { return _nodeList[1]->get_x() - _nodeList[0]->get_x(); };

 private:
  std::vector<Node *> _nodeList;
};

#endif
