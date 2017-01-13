#include "dag.h"

void
DAG::addEdge(long parentID, long childID) {
  if (vertexList.find(parentID) == vertexList.end()) addVertex(parentID);
  if (vertexList.find(childID) == vertexList.end()) addVertex(childID);
  vertexList[parentID]->addDownwind(childID);
  vertexList[childID]->addUpwind(parentID);
}

void
DAG::removeEdge(long parentID, long childID) {
  vertexList[parentID]->removeDownwind(childID);
  vertexList[childID]->removeUpwind(parentID);
}

void
DAG::getRoots() {
  if (roots.size() > 0) roots.clear();
  for (std::map<long,Vertex*>::iterator it=vertexList.begin();
       it != vertexList.end();
       it++) 
    if (it->second->upwindVertices.size() == 0) roots.push_back(it->first);
}

void
DAG::orderVertices() {
  if (orderedVertices.size() > 0) orderedVertices.clear();
  getRoots();
  while (roots.size() > 0) {
    long n = roots.back();
    orderedVertices.push_back(n);
    roots.pop_back();
    for (std::list<long>::iterator it=vertexList[n]->downwindVertices.begin();
         it != vertexList[n]->downwindVertices.end();
         it++) {
      removeEdge(n, *it);
      if (vertexList[*it]->upwindVertices.size() == 0) roots.push_back(*it);
    }
  }
}
