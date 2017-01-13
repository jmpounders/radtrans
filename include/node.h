#ifndef NODE_H
#define NODE_H

/// Vertex object
class Node
{
 public:
  Node(double x, double y=0.0, double z=0.0, unsigned int id=0) : _x(x), _y(y), _z(z), _id(id) {};
  double get_x() const { return _x; };
  double get_y() const { return _y; };
  double get_z() const { return _z; };
  unsigned int get_id() const { return _id; };

  void getCoords(double* coord) { coord[0] = _x; coord[1] = _y; coord[2] = _z; };

 private:
  unsigned int _id;
  double _x;
  double _y;
  double _z;
};

#endif
