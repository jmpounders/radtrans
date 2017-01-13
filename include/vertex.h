#ifndef VERTEX_H
#define VERTEX_H

#include <list>

class Vertex
{
 public:
  Vertex(long idIn) { id = idIn; };
  ~Vertex() {};

  bool operator<(const Vertex& v) const { return id<v.id; };
  bool operator==(const Vertex& v) { return id==v.id; };

  void addDownwind(long dwVert) { downwindVertices.push_back(dwVert); };
  void addUpwind(long uwVert) { upwindVertices.push_back(uwVert); };

  void removeDownwind(long dwVert) { downwindVertices.remove(dwVert); };
  void removeUpwind(long uwVert) { upwindVertices.remove(uwVert); };

  std::list<long> downwindVertices;
  std::list<long> upwindVertices;

  long id;
  
};

#endif
