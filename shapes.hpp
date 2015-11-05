// This file is a part of silhouette (2D Geometry suite)
//
// Copyright 2015 Taylor Fryett
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef SHAPES_HPP
#define SHAPES_HPP

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <stdexcept>

//##################################################################//
//####################### CoordPnt #################################//

/// \brief The Coordinate Point Class - Is not a drawable object - for utility
/// purposes
/// 
/// This is used to make drawable objectsa and perform abstract operations 
/// such as rotation about a point and forming lines for creating reflections
/// accross said line.
class CoordPnt {
private:
  double x; //!< The x Coordinate point
  double y; //!< The y Coordinate point

protected:

public:
  /// \brief default contructor - set Coordinate to the origin
  CoordPnt(void);

  /// \brief constructor with user specified x and y values
  ///
  /// @usrX the x Coordinate point (double) of the new object
  /// @usrY the y Coordinate point (double) of the new object
  CoordPnt(double usrX, double usrY);

  /// \brief copy initializer - creates a new identical object
  CoordPnt(const CoordPnt &p);

  /// \brief Allows two Coord objects to be summed together
  CoordPnt operator + (const CoordPnt &p);

  /// \brief Allows the coord1 += coord2 operation 
  // (as in coord1 = coord1 + coord2)
  CoordPnt operator += (const CoordPnt &p);

  /// \brief Allows two CoordPnt objects to be subtracted from one another
  CoordPnt operator - (const CoordPnt &p);

  /// \brief Allows the coord1 -= coord2 operation 
  /// (as in coord1 = coord1 - coord2)
  CoordPnt operator -= (const CoordPnt &p);

  /// \brief Allows an assingment of one coord object to another
  CoordPnt &operator = (const CoordPnt &p);

  /// \brief returns the value of the x field in the CoordPnt object
  double getX(void) const;

  /// \brief returns the value of the y field in the CoordPnt object
  double getY(void) const;

  /// \brief Resets the x value of the Coordinate to the user specified value
  ///
  /// @usrX the x (double) value to set the x Coordinate value to
  void setX(double usrX);

  /// \brief Resets the y value of the Coordinate to the user specified value
  ///
  /// @usrY the y (double) value to set the y Coordinate value to
  void setY(double usrY);

  /// \brief creates a string "(x, y)" that creates a printable 
  std::string toString(void) const;
};

/// \brief Operator to "scale" a point (move further from origin).
CoordPnt operator*(int multiplier, const CoordPnt& coord);

/// \brief Operator to "scale" a point (move further from origin).
CoordPnt operator*(const CoordPnt& coord, int multiplier);

/// \brief Operator to add two points.
CoordPnt operator+(const CoordPnt& coord1, const CoordPnt& coord2);

/// \brief Operator to subtract two points.
CoordPnt operator-(const CoordPnt& coord1, const CoordPnt& coord2);

//##################################################################//
//####################### LineSeg ##################################//

/// \brief The Line Class - it is not a drawable quantitiy - It is not a 
/// drawable object - for utility purposes
///
/// This type of object is used for boolean processes and for utility methods
/// such as reflecting an object across a line object.
class LineSeg {
private:
  std::vector<CoordPnt> limits; //!< defines the limits of line segment

protected:
public:

  /// \brief Creates a line defined by two end points defined by two CoordPnt 
  /// objects.
  ///
  /// @startSeg The beginning of the line segment
  /// @endSeg The end of the line segment
  LineSeg(CoordPnt startSeg, CoordPnt endSeg);
    
  /// \brief Reset the start of the line segment.
  ///
  /// @startSeg where the line segment will have its new start
  void setStartSeg(const CoordPnt &startSeg);

  /// \brief Restart the end of the line segment.
  ///
  /// @setEndSeg the new end of the line segment
  void setEndSeg(const CoordPnt &endSeg);

  /// \brief Retrieve a copy of the CoordPnt object that is at the start of the line
  /// segment.
  CoordPnt getStartPnt(void) const;

  /// Retrieve a copy of the CoordPnt object that is at the start of the line
  /// segment.
  CoordPnt getEndPnt(void) const;

  /// Return the length of the LineSeg object
  double length(void);

  /// Returns the rotation of the LineSeg object with respect to the positive
  /// x-axis (radians, -pi <= rad <= pi)
  double angleOffset(void);

  /// Returns the rotation of the LineSeg object with respect to the positive
  /// x-axis (degrees -180 <= degree <= 180 )
  double angleOffsetDegree(void);

  /// Returns the coordinate point that lies at the center of the line segment.
  CoordPnt getCenter(void) const;

}; // class LineSeg 

LineSeg operator+(const LineSeg& lin1, const LineSeg& lin2);

LineSeg operator-(const LineSeg& lin1, const LineSeg& lin2);

//----------------------------------------------------------------//
// Line utility functions
// contains functions which do not make sense to have as methods/
// class memebers or are better represented as functions than
// methods. They are all included in this header as they are 
// intimately related to the LineSeg class, and will be expected
// to be available when the LinSeg class is available.

/// \brief Determines if two lines intersect.
///
/// @lin1 The first line to consider.
/// @lin2 The second line to consider.
///
/// Returns true if lin1 crosses lin2. Note that we purposefully
/// implement this function so that two lines that are colinear
/// and overlap return false.
bool lineSegIntersect(LineSeg lin1, LineSeg lin2);

/// \brief Finds if the point pnt is on the line segment seg.
///
/// @seg The line segment in question.
/// @pnt The coordinate point in question.
///
/// Returns true if pnt is on seg, false otherwise.
bool pointOnLineSeg(LineSeg seg, CoordPnt pnt);

/// \brief Finds the orientation of the coordinate point to the line
/// segment.
///
/// @seg The line segment.
/// @pnt The coordinate point.
///
/// The function returns following values
/// 0 --> lin1.getStartPnt() pnt and lin1.getEndPnt() are colinear
/// 1 --> Clockwise
/// 2 --> Counterclockwise
int orientation(LineSeg seg, CoordPnt pnt);

//##################################################################//
//####################### Polygon ##################################//

class Polygon {
private:

protected:
  CoordPnt center; //!< Coordinate of the center of the polyon.
  int Layer; //!< The layer that the polygon belongs to. Usually corresponds to a fabrication step.
  int DataType; //!< Another number that details information, such as metal to be used.
  std::vector<CoordPnt> boundingBox; //!< The bounding box of the structure. Good for collision detection.
    
  //!< a vector of verticies that can be used to define the nature of varius
  // objects
  std::vector<CoordPnt> vertices;

  /// \brief Recalculates and relocates the center of the Polygon object.
  void findResetCenter(void);

  /// \brief Calcultes and returns the area of the polygon.
  double getArea(void) const;

  /// \brief Checks the vertices for their validity and conformity
  /// to the GDSII standards for BOUNDARY records, and then sets
  /// the polygons vertices to be usrVertices if all checks out.
  void setVertices(std::vector<CoordPnt> usrVertices);

  /// \brief Determines the bounding box from the "vertices" field.
  ///
  /// This method itterates over the "vertices" field to find the maximum and
  /// minimum displacements in the x and y directions. The resulting bounding
  /// box will have its sides parallel to the x and y axis.
  std::vector<CoordPnt> findBoundingBox(void) const;

  /// \brief Determines if the set of vertices passed in define a 
  ///        Polygon that contains an internal void.
  ///
  /// @usrVertices The vector of vertices that defines a polygon.
  ///
  /// As per the GDSII standard we must avoid any polygon (i.e 
  /// BOUNDARY record) that contains internal voids (same condition
  /// as having edges that cross one another. If such a condition
  /// exists this method returns true, otherwise it returns false.
  bool containsInternalVoid(std::vector<CoordPnt> usrVertices);

  /// \brief Thd default constructor for the Polygon class.
  ///
  /// Class is made protected so that child classes may inherit it 
  /// but it is not directly callable since it doess not make any
  /// shapes.
  Polygon(void);

public:

  /// \brief The full constructor for the polygon class.
  ///
  /// @usrVertices The list of vertices that defines the polygon.
  /// @usrLayer The layer at which the polygon will reside.
  /// @usrDataType the data type of which the polygon will belong.
  Polygon(std::vector<CoordPnt> usrVertices, int usrLayer, 
	  int usrDataType);

  /// \brief The full constructor for the polygon class.
  ///
  /// @usrVertices The list of vertices that defines the polygon.
  /// @usrLayer The layer at which the polygon will reside.
  Polygon(std::vector<CoordPnt> usrVertices, int usrLayer);

  /// \brief Generates a polygon defined by the vector of vertices passed in
  /// as the argument.
  ///
  /// @usrVertices The list of vertices that defines the polygon.
  ///
  /// Constructs the polygon from the vector of points passed into the 
  /// constructor. The order of the points matters greatly as the polygon 
  /// will be constructed by the points in the exact order they are passed in.
  Polygon(std::vector<CoordPnt> usrVertices);

  /// \brief Allows the user to rotate any Polygon objects.
  ///
  /// @rotatePnt The CoordPnt object that defines the point about which the 
  /// Polygon object should be rotated.
  /// @rotationAngle The angle (in radians) to rotate the Polygon object by.
  void rotate(CoordPnt rotatePnt, double rotationAngle);

  /// \brief Returns the int corresponding the the layer it is located on.
  int getLayer(void) const;
    
  /// \brief Allows the user to set the layer to the specified int.
  void setLayer(int newLayer);

  /// \brief Returns the int corresponding the the datatype of the polygon.
  int getDataType(void) const;
    
  /// \brief Allows the user to set the datatype to the specified int.
  void setDataType(int newDataType);

  /// \brief Returns a reference to the coordinate point list that comprises
  /// a polygon.
  std::vector<CoordPnt> &getVertices(void) const;

  /// \brief Determines if the coordinate point passed in lies
  /// within the polygon.
  ///
  /// @pnt The coordinate point in question.
  ///
  /// Returns true only if @pnt is inside the polygon ("this").
  /// In this definition inside does not include the boundaries,
  /// therefore if @pnt lies on this's boundary this method will
  /// return false. To note this algorithm uses the "casting ray"
  /// approach.
  bool pointInsidePolygon(CoordPnt pnt);

  friend void Polygon_addVertex(Polygon* poly, CoordPnt* newVertex);

  friend Polygon* Polygon_new(void);

}; // class polygon


/*
  Polygon operator+(Polygon& poly1, Polygon& poly2);

  Polygon operator-(Polygon& poly1, Polygon& poly2);
*/


//##################################################################//
//########################## Path ##################################//

class Path {
private:
  std::vector<CoordPnt> coordPath; //!< The points through which the path traverses.
  double pathWidth; //!< The width of the path.
  int pathType; //!< Determines the shape of end points. Can be 0, 1, or 2.
  int layer;
  int dataType;

protected:
    
public:
  //! Full constructor for the path class.
  Path(std::vector<CoordPnt> usrCoordPath, double usrPathWidth, 
       int usrPathType, int usrLayer, int usrDataType);

  Path(std::vector<CoordPnt> usrCoordPath, double usrPathWidth, 
       int usrPathType, int usrLayer);

  Path(std::vector<CoordPnt> usrCoordPath, double usrPathWidth, 
       int usrPathType);

  Path(std::vector<CoordPnt> usrCoordPath, double usrPathWidth);

  Path(std::vector<CoordPnt> usrCoordPath);

  int getPathType(void) const;

  void setPathType(int newPathType);

  double getPathWidth(void) const;

  void setPathWidth(double newPathWidth);

  std::vector<CoordPnt> getCoordPath(void) const;

  void setCoordPath(std::vector<CoordPnt> newCoordPath);

  void appendToCoordPath(CoordPnt nextCoordPnt);

  void setDataType(int newDataType);

  int getDataType(void) const;

  void setLayer(int newLayer);

  int getLayer(void) const;

};

//##################################################################//
//########################## Oval ##################################//

/// \brief A class to create ovals of any specification.
///
/// This is a class that generates ovals, typically defined by their center,
/// minor axis, major axis, and their orientation. Along with these quantities
/// the number of points that lies on the circumfrence also needs to be 
/// defined.
///
/// All objects of this class satisfy the general equation for an ellipse:
/// x^2/a^2 + y^2/b^2 = 1
/// where a and b are the lengths of the minor or major axis.
class Oval : public Polygon {
private:

protected:
  double eccentricity; //!< The ellipse eccentricity (0 for a circle, 1 for a line segment)
  double minorLength; //<! The minor axis length
  double majorLength; //!< The major axis length
  //! \brief Returns the eccentricity of the referenced oval object.
  double findEccentricity(void) const;
public:

  /// \brief Defines an arbitrary oval.
  /// 
  /// @majorAxis The LineSeg object that defines the placement and length of the
  /// major axis.
  /// @minorAxisLength The length of the minor axis.
  /// @numCoordPnts The number of CoordPnt objects that will represent the 
  /// discretization of the oval objct.
  ///
  /// Allows the user to construct an arbitrary oval by specifying the 
  /// major and minor axis. Since specifying the minor axis by a lineseg would 
  /// at best be redudant, instead we use a length to avoid such a problem.
  Oval(LineSeg majorAxis, double minorAxisLength, int numCoordPnts = 64,
       int layer = 1, int datatype = 0);

  /// \brief Defines an arbitrary oval.
  /// 
  /// @majorAxisLength The length of the major axis (centered on usrCenter).
  /// @minorAxisLength The length of the minor axis (centered on usrCenter).
  /// @numCoordPnts The number of CoordPnt objects that will represent the 
  /// discretization of the oval objct.
  /// 
  /// Allows the user to construct an arbitrary oval by specifying the 
  /// major and minor axis lengths. This should be the constructor to use the 
  /// majority of the time, unless you have some specific need to use another 
  /// type of constructor.
  /// raised.
  Oval(CoordPnt usrCenter, double majorAxisLength, double minorAxisLength,
       int numCoordPnts = 64, int layer = 1, int datatype = 0);


};

//##################################################################//
//######################### Circle #################################//

class Circle : public Oval {
private:

protected:

public:
    
  /// \brief The recommended constructor for this class. 
  ///
  /// @usrCenter The coordinate of the center of the circle.
  /// @usrRadius The radius of the constructed circle.
  Circle(CoordPnt usrCenter, double usrRadius);
    
  /// \brief The recommended constructor for this class. 
  ///
  /// @usrCenter The coordinate of the center of the circle.
  /// @usrRadius The radius of the constructed circle.
  /// @numVertices The number of vertices in the polygon to represent the circle.
  Circle(CoordPnt usrCenter, double usrRadius, int numCoordPnts = 64,
	 int layer = 1, int datatype = 0);


}; // End declaration of class "Circle"


//##################################################################//
//######################## Rectangle ###############################//

/// \brief Class to create Rectangles.
class Rectangle : public Polygon {
private:

protected:

public:
  /// \brief Construct the rectangle specified by width, height, and the 
  /// center of the rectangle.
  /// 
  /// @usrCenter The coordinate of the center of the new Rectangle object.
  /// @width The width of the Rectangle object
  /// @height The height of the Rectangle object.
  Rectangle(CoordPnt usrCenter, double width, double height, int layer,
	    int datatype);
};


//##################################################################//
//########################## Square ################################//

/// \brief Class for making squares.
///
/// This is a specialized version of the Rectangle class.
class Square : public Rectangle {
private:

protected:

public:
  /// \brief Constructs a square from the Rectangle constructor.
  Square(CoordPnt usrCenter, double sideLength, int layer, int datatype);

};


extern "C" {

  CoordPnt* CoordPnt_new(double usrX, double usrY);

  void CoordPnt_delete(CoordPnt* coord);

  LineSeg* LineSeg_new(CoordPnt* startSeg, CoordPnt* endSeg);

  void LineSeg_delete(LineSeg* lin);

  Polygon* Polygon_new(double* xcoords, double* ycoords, int numcoords,
		       int layer, int datatype);

  void Polygon_delete(Polygon* poly);

  Square* Square_new(CoordPnt* center, double sidelength, int layer,
		     int datatype);

  void Square_delete(Square* sq);
  
  Rectangle* Rectangle_new(CoordPnt* center, double width,
			   double height, int layer, int datatype);

  void Rectangle_delete(Rectangle* rect);
  
  Circle* Circle_new(CoordPnt* center, double radius, int numVertices,
		     int layer, int datatypeo);

  void Circle_delete(Circle* circ);

  Oval* Oval_new(CoordPnt* center, double majorAxisLength, double minorAxisLength,
		 int numCoordPnts, int layer, int datatype);

  void Oval_delete(Oval* ov);

  Path* Path_new(double* xcoords, double* ycoords, int numcoords,
		 double width, int layer, int datatype);

  void Path_delete(Path* path);
}

#endif // SHAPES_HPP

