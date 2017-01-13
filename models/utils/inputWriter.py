#!/bin/usr/python

from __future__ import print_function
import sys

class InputSet :
    def __init__(self, name) :
        self.name = name;
        self.subList = []
        self.data = {}

    def addData(self, varName, varValue) :
        self.data[varName] = varValue
                
    def add(self, subInput) :
        self.subList.append(subInput)
        
    def write(self, fileName=sys.stdout, indent=0) :
        spaces = " "*indent
        print(spaces + "<" + self.name + ">", file=fileName)
        for item in self.data :
            print(spaces + "  " + item + " ", end="", file=fileName)
            if (type(self.data[item]) is list) :
                print(', '.join(map(str, self.data[item])), file=fileName)
            else :
                print(self.data[item], file=fileName)
        
        for set in self.subList :
            set.write(fileName, indent+2)
            
        print(spaces + "</" + self.name + ">", file=fileName)
        print(file=fileName)
