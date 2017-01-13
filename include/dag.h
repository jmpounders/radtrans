#ifndef DAG_H
#define DAG_H

#include "vertex.h"

#include <map>
#include <deque>

class DAG
{
 public:
  DAG();
  ~DAG();

  void addVertex(long id) { vertexList[id] = new Vertex(id); };
  void addEdge(long parentID, long childID);
  void removeEdge(long parentID, long childID);
  void getRoots();
  void orderVertices();
  
 private:
  std::map<long,Vertex*> vertexList;
  std::deque<long> roots;
  std::deque<long> orderedVertices;
};

#endif
