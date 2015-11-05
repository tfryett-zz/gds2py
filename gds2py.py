#!/usr/bin/env python

# This file is a part of gdscad (2D Geometry suite)               
#                                                                 
# Copyright 2015 Taylor Fryett                                    
#                                                                 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np
import ctypes
lib = ctypes.cdll.LoadLibrary('./build/libgdscad.so')

class CoordPnt:
    def __init__(self, x, y):
        self.obj = lib.CoordPnt_new(ctypes.c_double(x),
                                    ctypes.c_double(y))
        self._x = x
        self._y = y

    def __enter__(self):
        return self

    def __exit__(self):
        lib.CoordPnt_delete(self.obj)

    def __str__(self):
        return "(" + str(self._x) + ", " + str(self._y) + ")"

class LineSeg:
    def __init__(self, startCoord, endCoord):
        self.obj = lib.LineSeg_new(startCoord.obj, endCoord.obj)
        self._startCoord = startCoord
        self._endCoord = endCoord

    def __enter__(self):
        return self

    def __exit__(self):
        lib.LineSeg_delete(self.obj)

    def __str__(self):
        return str(self._startCoord) + "--->" + str(self._endCoord)

class Polygon:
    def __init__(self, coordList, layer = 1, datatype = 0):
        x = np.array([coord._x for coord in coordList], dtype = np.float)
        y = np.array([coord._y for coord in coordList], dtype = np.float)
        self.obj = lib.Polygon_new(ctypes.c_void_p(x.ctypes.data),
                                   ctypes.c_void_p(y.ctypes.data),
                                   ctypes.c_int(len(x)), ctypes.c_int(layer),
                                   ctypes.c_int(datatype))
        self._coordList = coordList
        self._layer = layer
        self._datatype = datatype

    def __enter__(self):
        return self

    def __exit__(self):
        lib.Polygon_delete(self.obj)

class Oval:
    def __init__(self, center, minorRadius, majorRadius,
                 numPoints = 64, layer = 1, datatype = 0):
        self.obj = lib.Oval_new
        self._center = center
        self._minorRadius = minorRadius
        self._majorRadius = majorRadius
        self._numPoints = numPoints
        self._layer = layer
        self._datatype = datatype

    def __enter__(self):
        return self

    def __exit__(self):
        lib.Oval_delete(self.obj)

class Circle:
    def __init__(self, center, radius, numPoints = 64, layer = 1,
                 datatype = 0):
        self.obj = lib.Circle_new(center, radius, numPoints, layer,
                                  datatype)
        self._center = center
        self._radius = radius
        self._numPoints = numPoints
        self._layer = layer
        self._datatype = datatype

    def __enter__(self):
        return self

    def __exit__(self):
        lib.Circle_delete(self.obj)

class Square:
    def __init__(self, center, sidelength, layer = 1, datatype = 0):
        self.obj = lib.Square_new(center.obj, sidelength, layer, datatype)
        self._center = center
        self._sideLength = sideLength
        self._layer = layer
        self._datatype = datatype

    def __enter__(self):
        return self

    def __exit__(self):
        lib.Square_delete(self.obj)

class Rectangle:
    def __init__(self, center, sideXLength, sideYLength, layer = 1,
                 datatype = 0):
        self.obj = lib.Rectangle_new(center.obj, sideXLength,
                                     sideYLength, layer, datatype)
        self._center = center
        self._sideXLength = sideXLength
        self._sideYLength = sideYLength
        self._layer = layer
        self._datatype = datatype

    def __enter__(self):
        return self

    def __exit__(self):
        lib.Rectangle_delete(self.obj)

class Path:
    def __init__(self, coordList, width, layer, datatype):
        x = np.array([coord._x for coord in coordList], dtype = np.float)
        y = np.array([coord._y for coord in coordList], dtype = np.float)
        self.obj = lib.Path_new(ctypes.c_void_p(x.ctypes.data),
                                   ctypes.c_void_p(y.ctypes.data),
                                   ctypes.c_int(len(x)),
                                   ctypes.c_double(width),
                                   ctypes.c_int(layer),
                                   ctypes.c_int(datatype))
        self._coordList = coordList
        self._width = width
        self._layer = layer
        self._datatype = datatype

    def __enter__(self):
        return self

    def __exit__(self):
        lib.Path_delete(self.obj)

class Cell:
    def __init__(self, cellname):
        self.obj = lib.Cell_new(cellname)
        self._cellname = cellname
        self._polygonList = []
        self._refCellList = []
        self._cellArrayList= []
        self._pathList = []

    def __enter__(self):
        return self

    def __exit__(self):
        lib.Cell_delete(self.obj)

    def addPolygon(self, polygon):
        self._polygonList.append(polygon)
        lib.Cell_addPolygon(self.obj, polygon.obj)
        
    def addCellReference(self, cellRef):
        self._refCellList.append(cellRef)
        lib.Cell_addCellReference(self.obj, cellRef.obj)

    def addCellArray(self, cellArray):
        self._cellArrayList.append(cellArray)
        lib.Cell_addCellArray(self.obj, cellArray.obj)

    def addPath(self, path):
        self._pathList.append(path)
        lib.Cell_addPath(self.obj, path.obj)

class CellReference:
    def __init__(self, referencedCell, centerPos, xspacing, yspacing,
                 magnification, rotation):
        self.obj = lib.CellReference_new(referencedCell.obj, centerPos.obj,
                                         xspacing, yspacing, 
                                         magnification, rotation)
        self._referencedCell = referencedCell
        self._centerPos = centerPos
        self._xspacing = xspacing
        self._yspacing = yspacing
        self._magnification = magnification
        self._rotation = rotation

        def __enter__(self):
            return self

        def __exit__(self):
            lib.CellReference_delete(self.obj)

class CellArray:
    def __init__(self, referencedCell, startingPos, numCols, numRows,
                 magnification, rotation):
        self.obj = lib.CellArray_new(referencedCell.obj, startingPos.obj,
                                     numCols, numRows, magnification,
                                     rotation)
        self._referencedCell = referencedCell
        self._startingPos = startingPos
        self._magnification = magnification
        self._rotation = rotation
        self._numCols = numCols
        self._numRows = numRows

    def __enter__(self):
        return self

    def __exit__(self):
        lib.CellArray_delete(self.obj)

class Layout:
    def __init__(self):
        lib.Layout_new.argtypes = None
        lib.Layout_new.restype = ctypes.c_void_p
        self.obj = lib.Layout_new()
        
    def __enter__(self):
        return self

    def __exit__(self):
        lib.Layout_delete(self.obj)

    def addCell(self, cell):
        lib.Layout_addCell(self.obj, cell.obj)

    def write(self, filename):
        lib.Layout_write(self.obj, filename)
        

#### For testing purposes

if __name__ == "__main__":
    fc = CoordPnt(1,1)
    sc = CoordPnt(2,1)
    tc = CoordPnt(2,2)
    foc = CoordPnt(1,2)
    
    polylist = [fc,sc,tc,foc]
    poly = Polygon(polylist)

    ovl = Oval(fc, 5., 3.)

    cell = Cell("mycell")
    cell.addPolygon(poly)

    lay = Layout()
    lay.addCell(cell)

    lay.write("test.gds")
