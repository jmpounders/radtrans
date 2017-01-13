#include "element.h"

LineElement::LineElement(Node *n1, Node *n2) : 
  Element(2)
{
  _nodeList.reserve(2);
  _nodeList.push_back(n1);
  _nodeList.push_back(n2);
}
