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

#include "shapes.hpp"

//##################################################################//
//####################### CoordPnt #################################//

/// The default constructor. Creates a coordinate point at the origin
///    
/// Example use:
//  coord origin = CoordPnt();
CoordPnt::CoordPnt(void) {
  x = 0;
  y = 0;
}
  
CoordPnt::CoordPnt(double usrX, double usrY) {
  x = usrX;
  y = usrY;
}

// constructor that copies the passed in coord point object
CoordPnt::CoordPnt(const CoordPnt &p) {
  x = p.x;
  y = p.y;
}

CoordPnt CoordPnt::operator + (const CoordPnt &p) {
  return CoordPnt(x + p.getX(), y + p.getY());
}

CoordPnt CoordPnt::operator += (const CoordPnt &p) {
  x += p.getX();
  y += p.getY();
  return *this;
}

CoordPnt CoordPnt::operator - (const CoordPnt &p) {
  return CoordPnt(x - p.getX(), y - p.getY());
}

CoordPnt CoordPnt::operator -= (const CoordPnt &p) {
  x -= p.getX();
  y -= p.getY();
  return *this;
}

CoordPnt &CoordPnt::operator = (const CoordPnt &p) {
  x = p.x;
  y = p.y;
  return *this;
}


// returns a copy of the field "x"
double CoordPnt::getX(void) const {
  return this->x;
}

// allows user to reset the field "x"
void CoordPnt::setX(double usrX) {
  this->x = usrX;
}

// returns a copy of the field "y"
double CoordPnt::getY(void) const {
  return this->y;
}

// allows user to reset the field "y"
void CoordPnt::setY(double usrY) {
  this->y = usrY;
}

// allows the user to turn a CoordPnt object into a print friendly string
std::string CoordPnt::toString(void) const {
  std::ostringstream output;
  output << "(" << this->getX() << ", " << this->getY() << ")";
  return output.str();
}

CoordPnt operator+(const CoordPnt& coord1, const CoordPnt& coord2) {
  return CoordPnt(coord1.getX() + coord2.getX(), 
		  coord1.getY() + coord2.getY());
}
  
CoordPnt operator-(const CoordPnt& coord1, const CoordPnt& coord2) {
  return CoordPnt(coord1.getX() - coord2.getX(), 
		  coord1.getY() - coord2.getY());
}

CoordPnt operator*(int multiplier, const CoordPnt& coord) {
  return CoordPnt(coord.getX()*multiplier, coord.getY()*multiplier);
}

CoordPnt operator*(const CoordPnt& coord, int multiplier) {
  return CoordPnt(coord.getX()*multiplier, coord.getY()*multiplier);
}

//##################################################################//
//####################### LineSeg ##################################//

// normal constructor 
LineSeg::LineSeg(CoordPnt startSeg, CoordPnt endSeg) {
  limits.push_back(startSeg);
  limits.push_back(endSeg);
}

// sets the starting part of the line segment
void LineSeg::setStartSeg(const CoordPnt &startSeg) {
  this->limits[0] = startSeg;
}

// sets the ending part of the line segment
void LineSeg::setEndSeg(const CoordPnt &endSeg) {
  this->limits[1] = endSeg;
}

// returns the start of the line segment
CoordPnt LineSeg::getStartPnt() const {
  CoordPnt returnPnt = limits[0];
  return returnPnt;
}

// returns the end of the line segment
CoordPnt LineSeg::getEndPnt() const {
  CoordPnt returnPnt = limits[1];
  return returnPnt;
}

// Returns the length of the line segment
double LineSeg::length() {
  CoordPnt startPnt = this->getStartPnt();
  CoordPnt endPnt = this->getEndPnt();
  double xSeparation = std::abs(startPnt.getX() - endPnt.getX());
  double ySeparation = std::abs(startPnt.getY() - endPnt.getY());
  // use the standard distance equation:
  //     d = sqrt(x^2 + y^2)
  return std::sqrt(std::pow(xSeparation, 2) + std::pow(ySeparation, 2));
}

// Returns the rotation of the LineSeg object with respect to the positive
// x-axis (radians)
double LineSeg::angleOffset() {
  CoordPnt startPnt = this->getStartPnt();
  CoordPnt endPnt = this->getEndPnt();
  double xSeparation = std::abs(startPnt.getX() - endPnt.getX());
  double ySeparation = std::abs(startPnt.getY() - endPnt.getY());
  // atan is the arctangent function
  return std::atan2(ySeparation, xSeparation);
}

// Returns the rotation of the LineSeg object with respect to the positive
// x-axis in degrees
double LineSeg::angleOffsetDegree() {
  double rotInRad = angleOffset();
  // acos is the arccosine function. This provides the most precise value of
  // the constant pi allowed by the user's compiler.
  double pi = std::acos(-1.); 
  // Degrees = 180/pi * Radians
  return rotInRad*180./pi;
}

// Returns a CoordPnt that is at the average point of the two limit points of
// the LineSeg object.
CoordPnt LineSeg::getCenter() const {
  return CoordPnt((limits[0].getX() + limits[1].getX())/2., 
		  (limits[0].getY() + limits[1].getY())/2.);
}

LineSeg operator+(const LineSeg& lin1, const LineSeg& lin2) {
  return LineSeg(lin1.getStartPnt() + lin2.getStartPnt(),
		 lin1.getEndPnt() + lin2.getEndPnt());
}

LineSeg operator-(const LineSeg& lin1, const LineSeg& lin2) {
  return LineSeg(lin1.getStartPnt() - lin2.getStartPnt(),
		 lin1.getEndPnt() - lin2.getEndPnt());
}

bool pointOnLineSeg(LineSeg seg, CoordPnt pnt) {
  if ( pnt.getX() <= std::max(seg.getStartPnt().getX() , seg.getEndPnt().getX()) &&
       pnt.getX() >= std::min(seg.getStartPnt().getX() , seg.getEndPnt().getX()) &&
       pnt.getY() <= std::max(seg.getStartPnt().getY() , seg.getEndPnt().getY()) &&
       pnt.getY() >= std::min(seg.getStartPnt().getY() , seg.getEndPnt().getY())
       )
    return true;
  else
    return false;
}
 
int orientation(LineSeg seg, CoordPnt pnt) {
  // See 10th slides from following link for derivation of the formula
  // http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
  CoordPnt firstToPnt = pnt - seg.getStartPnt();
  CoordPnt lastToPnt = seg.getEndPnt() - pnt;
  int val = firstToPnt.getY()*lastToPnt.getX() -
    firstToPnt.getX()*lastToPnt.getY();
  //    int val = (q.y - p.y) * (r.x - q.x) -
  //         (q.x - p.x) * (r.y - q.y);
 
  if (val == 0) // colinear 
    return 0;  
  else
    return (val > 0)? 1: 2; // clock or counterclock wise
}

bool lineSegIntersect(LineSeg lin1, LineSeg lin2) {
  int orient1 = orientation(lin1, lin2.getStartPnt());
  int orient2 = orientation(lin1, lin2.getEndPnt());
  int orient3 = orientation(lin2, lin1.getStartPnt());
  int orient4 = orientation(lin2, lin1.getEndPnt());

  // General case:
  // The line segments are at different angles and they intersect.
  if (orient1 != orient2 && orient3 != orient4)
    return true;
  else 
    return false; 

  /*
    Provided for reference if there is ever the need to determine
    if two parallel lines intersect.

    // Special Cases
    // lin1 and lin2.getStartpnt() are coliniar and 
    // lin2.getStartpnt() lies on lin1
    if (orient1 == 0 && pointOnLineSeg(lin1, lin2.getStartPnt()))
    return true;
 
    // lin1 and lin2.getEndpnt() are coliniar and 
    // lin2.getEndpnt() lies on lin1
    if (orient2 == 0 && pointOnLineSeg(lin1, lin2.getEndPnt())) 
    return true;
 
    // lin2 and lin1.getStartpnt() are coliniar and 
    // lin1.getStartpnt() lies on lin2
    if (orient3 == 0 && pointOnLineSeg(lin2, lin1.getStartPnt())) 
    return true;
 
    // lin2 and lin1.getEndpnt() are coliniar and 
    // lin1.getEndpnt() lies on lin2
    if (orient4 == 0 && pointOnLineSeg(lin2, lin1.getEndPnt())) 
    return true;
  */
}

//##################################################################//
//####################### Polygon ##################################//

Polygon::Polygon(std::vector<CoordPnt> usrVertices, int usrLayer, 
		 int usrDataType) {
  this->setVertices(usrVertices);
  this->setLayer(usrLayer);
  this->setDataType(usrDataType);
}

Polygon::Polygon(std::vector<CoordPnt> usrVertices, int usrLayer) :
  Polygon(usrVertices, usrLayer, 0) {}

Polygon::Polygon(std::vector<CoordPnt> usrVertices) :
  Polygon(usrVertices, 1) {}

Polygon::Polygon() {
  this->setLayer(1);
  this->setDataType(0);
}

std::vector<CoordPnt> Polygon::findBoundingBox() const {
  // Initialize the variables to store the maximum and minimum values of 
  // x and y that we find
  int maxY = vertices[0].getY();
  int minY = vertices[0].getY();
  int maxX = vertices[0].getX();
  int minX = vertices[0].getX();
  // Test to find out what the minimum and maximum extent of the x and y
  // are and store them in the corresponding variables
  for (std::vector<CoordPnt>::const_iterator it = vertices.begin(); 
       it != vertices.end(); ++it) {
    if (maxY < it->getY())
      maxY = it->getY();
    else if (minY > it->getY())
      minY = it->getY();
    if (maxX < it->getX())
      maxX = it->getX();
    else if (minX > it->getX())
      minX = it->getX();
  }
  // create the std::vector to store the appropriate points in
  std::vector<CoordPnt> newBoundingBox;
  newBoundingBox.push_back(CoordPnt(minX, maxY)); // upper left corner
  newBoundingBox.push_back(CoordPnt(maxX, maxY)); // upper right corner
  newBoundingBox.push_back(CoordPnt(maxX, minY)); // lower right corner
  newBoundingBox.push_back(CoordPnt(minX, minY)); // lower left corner
  // return completed bounding box
  return newBoundingBox;
}

// Rotates the Polygon objects about rotatePnt by rotationAngle radians.
void Polygon::rotate(CoordPnt rotatePnt, double rotationAngle) {
  // Create the rotation matrix
  double rotationMatrix[2][2];
  rotationMatrix[0][0] = std::cos(rotationAngle);
  rotationMatrix[0][1] = -1.*std::sin(rotationAngle);
  rotationMatrix[1][0] = std::sin(rotationAngle);
  rotationMatrix[1][1] = std::cos(rotationAngle);
  // create a position vector (column) for rotatePnt
  double rotateVec[2];
  rotateVec[0] = rotatePnt.getX();
  rotateVec[1] = rotatePnt.getY();
  // create a position column vector for the vertices
  double vertexVec[2];
  // Now rotate each vertex point about the rotatePnt
  for (std::vector<CoordPnt>::iterator it = vertices.begin();
       it != vertices.end(); ++it) {
    // set up the position vector for this vertex
    vertexVec[0] = it->getX();
    vertexVec[1] = it->getY();
    // make the rotation point the origin
    vertexVec[0] -= rotateVec[0];
    vertexVec[1] -= rotateVec[1];
    // rotate the vertex
    vertexVec[0] = rotationMatrix[0][0]*vertexVec[0] + 
      rotationMatrix[0][1]*vertexVec[1];
    vertexVec[1] = rotationMatrix[1][0]*vertexVec[0] + 
      rotationMatrix[1][1]*vertexVec[1];
    // now shift the matrix back
    vertexVec[0] += rotateVec[0];
    vertexVec[1] += rotateVec[1];
    // now reassign the vertices to their new locations
    it->setX(vertexVec[0]);
    it->setY(vertexVec[1]);
  }
  // update the center field
  this->findResetCenter();
}
   
// Resets the center field
void Polygon::findResetCenter() {
  int numVertices = 0;
  double xSum = 0.;
  double ySum = 0.;
  for (std::vector<CoordPnt>::iterator it = vertices.begin(); 
       it != vertices.end(); ++it) {
    xSum += it->getX();
    ySum += it->getY();
    numVertices++;
  }
  // Reset the center field to be the average values from the vertices
  this->center = CoordPnt(xSum/numVertices, ySum/numVertices);
}

int Polygon::getLayer() const {
  return this->Layer;
}

void Polygon::setLayer(int newLayer) {
  if (newLayer >= 0 && newLayer < 64)
    this->Layer = newLayer;
  else {
    std::stringstream errorMsg;
    errorMsg << "The layer assigned out of range. User attempted to"
	     << " assign a value of " << newLayer << ", but " 
	     << "layers must be equal to or between 0 and 63.\n";
    std::invalid_argument(errorMsg.str());
  }
}

int Polygon::getDataType() const {
  return this->DataType;
}

void Polygon::setDataType(int newDataType) {
  if (newDataType >= 0 && newDataType < 64)
    this->Layer = newDataType;
  else {
    std::stringstream errorMsg;
    errorMsg << "The data type assigned out of range. User attempted"
	     << " to assign a value of " << newDataType << ", but"
	     << " data types must be equal to or between 0 and 63.\n";
    std::invalid_argument(errorMsg.str());
  }
}

std::vector<CoordPnt> &Polygon::getVertices() const {
  return const_cast<std::vector<CoordPnt> &> (this->vertices);
}

bool Polygon::containsInternalVoid(std::vector<CoordPnt> usrVertices) {
  // test to make sure what was passed in couuld actually be a polygon
  // if we don't the next test will always throw an error and it will
  // be harder to determine what went wrong than if we catch it here
  // like this.
  if (usrVertices.size() < 3) {
    std::stringstream errorMsg;
    errorMsg << "The number of vertices must be equal to or larger"
	     << " than three. User passed in a vector of "
	     << usrVertices.size() << ".\n";
    throw std::invalid_argument(errorMsg.str());
  }

  // even a single intersection implies an internal void. neither
  // are allowed by gdsII standards so quit at the first one found.
  uint upperLimit = usrVertices.size();
  for (uint i = 0; i <= upperLimit; i++) {
    LineSeg lin1(usrVertices[i], usrVertices[i + 1]);
    for (uint j = i + 2; j <= upperLimit; j++) {
      LineSeg lin2(usrVertices[j], usrVertices[j + 1]);
      if (lineSegIntersect(lin1, lin2))
	return true;
    }
  }

  return false;
}

void Polygon::setVertices(std::vector<CoordPnt> usrVertices) {
  if (usrVertices.size() < 3 || usrVertices.size() > 199) {
    std::stringstream errorMsg;
    errorMsg << "Each polygon may only have 3-199 vertices. User"
	     << " tried to initiate a polygon with " 
	     << usrVertices.size() << " vertices.\n";
    throw std::invalid_argument(errorMsg.str());
  }

  if (!this->containsInternalVoid(usrVertices))
    this->vertices = usrVertices;
  // update the bounding box
  this->boundingBox = this->findBoundingBox();
}

bool Polygon::pointInsidePolygon(CoordPnt pnt) {
  // for the ray casting techneique we need to use a point that is
  // garunteed to be outside of the polygon (we want to reuse the
  // algorithm for determining the intersection of two line segments)
  // Use the point that has the same y coordinate as @pnt but has
  // an x coordinate that is one unit towards +inf than the largest
  // x coordinate in the vertex list.
  uint numVertices = this->boundingBox.size();
  int maxXCoord = std::abs(this->boundingBox[0].getX());
  int maxYCoord = std::abs(this->boundingBox[0].getY());
  for (uint i = 1; i < numVertices; i++) {
    if (std::abs(this->boundingBox[i].getX()) > maxXCoord) 
      maxXCoord = std::abs(this->boundingBox[i].getX());
    if (std::abs(this->boundingBox[i].getY()) > maxYCoord)
      maxYCoord = std::abs(this->boundingBox[i].getY());
  }
  CoordPnt rayCastCoord = CoordPnt(maxXCoord*1.1, pnt.getY());
  LineSeg rayCast = LineSeg(pnt, rayCastCoord);

  // now detect if the ray (rayCast) intersects any of the sides of
  // the polygon. If rayCast intersects an odd number of sides then
  // @pnt is inside the polygon, otherwise it is outside.
  numVertices = this->vertices.size();
  int intersectionCount = 0;
  for (uint i = 0; i < numVertices - 1; i++) {
    LineSeg polySide = LineSeg(this->vertices[i], this->vertices[i + 1]);
    if (lineSegIntersect(rayCast, polySide)) // linesegintersect in line.cxx & line.hxx
      intersectionCount++;
    // If it hits a vertex determining what happes can be tricky, so
    // let us just first solve for a ray that does not intersect with
    // any vertices and then solve. This is more computationally
    // heavy so save it for if the first attempt did not work.
    if (pointOnLineSeg(rayCast, this->vertices[i])) { // from line.cxx line.hxx
      // first cast a ray from @pnt to each of the vertices and then
      // determine the place where the maximum angular separation is
      // so we can aim for half way inbetween those vertices which
      // will give us the best chance of success.
      double rayCast1Angle = LineSeg(pnt, this->vertices[0]).angleOffset();
      double rayCast2Angle = LineSeg(pnt, this->vertices[1]).angleOffset();
      double maxAngularSep = std::abs(rayCast2Angle - rayCast1Angle);
      double bestAngle = (rayCast2Angle + rayCast1Angle)/2.0;
      for (uint j = 1; j < numVertices - 1; j++) {
	rayCast1Angle = LineSeg(pnt, this->vertices[j]).angleOffset();
	rayCast2Angle = LineSeg(pnt, this->vertices[j + 1]).angleOffset();
	double curAngleSep = std::abs(rayCast2Angle - rayCast1Angle);
	if (curAngleSep > maxAngularSep)
	  bestAngle = (rayCast1Angle + rayCast2Angle)/2.0;
      }

      double boundingCircRadius = std::sqrt(pow(maxXCoord, 2) +
					    pow(maxYCoord, 2));
      CoordPnt rayPnt = CoordPnt(boundingCircRadius*std::cos(bestAngle),
				 boundingCircRadius*std::sin(bestAngle));
      rayCast = LineSeg(pnt, rayPnt);

      // retry the same algorithm with the new rayCast line
      intersectionCount = 0; // restart the count
      for (uint j = 0; j < numVertices - 1; j++) {
	LineSeg polySide = LineSeg(this->vertices[j], this->vertices[j + 1]);
	if (lineSegIntersect(rayCast, polySide)) // linesegintersect in line.cxx & line.hxx
	  intersectionCount++;
      }
      break; // get out of the loop since we have our answer now.
    }
  }

  return (intersectionCount % 2 == 1);
}

/*
// Uses the Weiler-Atherton Algorithm
Polygon operator-(Polygon& poly1, Polygon& poly2) {

// We will use this struct to keep track of the intersection
// points. It is important to mark the points as either leaving
// or entering the clipping region (@poly2). We will use this
// bool value to mark when to add points from @poly1 or @poly2
// to the clipped polygon we will return. To that end we will 
// keep track of the indices of the points that come right before
// the intersection for both polygons to speed the algorithm up.
struct intersectPnt {
Polygon Pnt;
bool enteringPoly2;
int poly1FirstIndex;
int poly2FirstIndex;
};
    
// where we will store all of the points that are the
// intersection between sides of the two different polygons
std::vector<intersectPnt> intersectionPoints;

int numVerticesPoly1 = poly1.getVertices().size();
int numVerticesPoly2 = poly2.getVertices().size();
for (int p1VertexNum = 0; p1VertexNum < numVerticesPoly1; p1VertexNum++) {
int secondVertexNumPoly1; // use to make sure we wrap around
if (p1VertexNum == numVerticesPoly1 - 1)
secondVertexNumPoly1 = 0;
else
secondVertexNumPoly1 = p1VertexNum + 1;
LineSeg p1Lin = LineSeg(poly1.getVertices()[p1VertexNum],
poly1.getVertices()[secondVertexNumPoly1]);
for (int p2VertexNum = 0; p2VertexNum < numVerticesPoly2; p2VertexNum++) {
int secondVertexNumPoly2;
if (p2VertexNum == numVerticesPoly2 - 1)
secondVertexNumPoly2 = 0;
else
secondVertexNumPoly2 = p1VertexNum + 1;
LineSeg p2Lin = LineSeg(poly1.getVertices()[p1VertexNum],
poly1.getVertices()[secondVertexNumPoly2]);

// test if p2Lin and p1Lin intersect
if (lineSegIntersect(p1Lin, p2Lin)) {
// take advantage of the fact that the function
// lineSegIntersect does not allow both LineSegs to have
// the same slope
double errorMargin = 1e-6;
// i.e. p1Lin has slope of infinity. We know p2Lin must
// have a defined slope
if (std::abs(p1Lin.getStartPnt().getX() -
p1Lin.getEndPnt().getX()) < errorMargin) {
	    
}
}
	
} // ends loop over polygon 2 vertices
} // ends loop over polygon 1 vertices

}
*/

//##################################################################//
//########################## Path ##################################//

Path::Path(std::vector<CoordPnt> usrCoordPath, double usrPathWidth, 
	   int usrPathType, int usrLayer, int usrDataType) {
  this->setCoordPath(usrCoordPath);
  this->pathWidth = usrPathWidth;
  this->setPathType(usrPathType); // check for validity and then set value
  this->setLayer(usrLayer);
  this->setDataType(usrDataType);
}

Path::Path(std::vector<CoordPnt> usrCoordPath, double usrPathWidth, 
	   int usrPathType, int usrLayer) :
  Path::Path(usrCoordPath, usrPathWidth, usrPathType, usrLayer, 0) {}

Path::Path(std::vector<CoordPnt> usrCoordPath, double usrPathWidth, 
	   int usrPathType) :
  Path::Path(usrCoordPath, usrPathWidth, usrPathType, 0) {}

Path::Path(std::vector<CoordPnt> usrCoordPath, double usrPathWidth) :
  Path::Path(usrCoordPath, usrPathWidth, 0) {}

Path::Path(std::vector<CoordPnt> usrCoordPath) :
  Path::Path(usrCoordPath, 0) { 
}

void Path::setPathType(int newPathType) {
  if (newPathType >= 0 && newPathType <= 2)
    this->pathType = newPathType;
  else  {
    std::ostringstream warning;
    warning << "Only path types 0, 1, and 2 are defined. User " 
	    << "tried to use a path type of " << newPathType 
	    << "." << std::endl;
    throw std::invalid_argument(warning.str());
  }
}

int Path::getPathType(void) const {
  return this->pathType;
}

double Path::getPathWidth(void) const {
  return this->pathWidth;
}

void Path::setPathWidth(double newPathWidth) {
  this->pathWidth = newPathWidth;
}

std::vector<CoordPnt> Path::getCoordPath(void) const {
  return this->coordPath;
}

void Path::setCoordPath(std::vector<CoordPnt> newCoordPath) {
  if (newCoordPath.size() > 1 && newCoordPath.size() < 200) 
    this->coordPath = newCoordPath;
  else {
    std::stringstream errorMsg;
    errorMsg << "The number of coordinates in a Path must"
	     << " be in the range 2-199. User tried to set" 
	     << " a path with " << newCoordPath.size() << ".\n";
    throw std::invalid_argument(errorMsg.str());
  }
}

void Path::appendToCoordPath(CoordPnt nextCoordPnt) {
  if (this->coordPath.size() < 199) // can have up to 199 coordiantes
    this->coordPath.push_back(nextCoordPnt);
  else {
    std::stringstream errorMsg;
    errorMsg << "Path has its maximum number of coordinates already.\n";
    throw std::invalid_argument(errorMsg.str());
  }
}

void Path::setDataType(int newDataType) {
  if (newDataType >= 0 && newDataType < 64)
    this->dataType = newDataType;
  else {
    std::stringstream errorMsg;
    errorMsg << "Invalid data type. Data types must be "
	     << "in the range 0-63, but user specified "
	     << "data type of " << newDataType << ".\n";
    throw std::invalid_argument(errorMsg.str());
  }
}

int Path::getDataType(void) const {
  return this->dataType;
}

void Path::setLayer(int newLayer) {
  if (newLayer >= 0 && newLayer < 64)
    this->dataType = newLayer;
  else {
    std::stringstream errorMsg;
    errorMsg << "Invalid layer. Layers must be "
	     << "in the range 0-63, but user specified "
	     << "layer as " << newLayer << ".\n";
    throw std::invalid_argument(errorMsg.str());
  }
}

int Path::getLayer(void) const {
  return this->layer;
}

//##################################################################//
//########################## Oval ##################################//

Oval::Oval(CoordPnt usrCenter, double majorAxisLength, double minorAxisLength,
	   int numCoordPnts, int layer, int datatype) {
  // Throw an exception if numCoordPnts is not reasonably high (at least 4)
  if (numCoordPnts < 4) 
    throw std::invalid_argument("numCoordPnts is too small to be effective.");
  if (majorAxisLength < minorAxisLength)
    throw std::invalid_argument("The major radius must be larger than the minor radius");
  // The angular separation of the vertices
  const double PI = std::acos(-1); // compiler limited representation of pi
  double angularSep = 2*PI/numCoordPnts;
  double angularPos = 0.;
  for (int i = 0; i < numCoordPnts; i++) {
    double localRadius = majorAxisLength*minorAxisLength/
      std::sqrt(std::pow(minorAxisLength*std::cos(angularPos), 2) +
		std::pow(majorAxisLength*std::sin(angularPos), 2));
    double vertexX = localRadius*std::cos(angularPos);
    double vertexY = localRadius*std::sin(angularPos);
    vertices.push_back(CoordPnt(vertexX + usrCenter.getX(), 
				vertexY + usrCenter.getY()));
    angularPos += angularSep;
  }
    
  // set all of the oval field values for this object
  // area of an ellipse is a*b*pi (acos(-1) = most precise version of pi)
  this->minorLength = minorAxisLength;
  this->majorLength = majorAxisLength;
  this->eccentricity = findEccentricity();
  this->center = usrCenter;
  this->setDataType(datatype);
  this->setLayer(layer);
}

Oval::Oval(LineSeg majorAxis, double minorAxisLength, int numCoordPnts,
	   int layer, int datatype) :
  Oval(majorAxis.getCenter(), minorAxisLength, numCoordPnts, layer, datatype) {
  this->rotate(majorAxis.getCenter(), majorAxis.angleOffset());
}

// The eccentricity is defined as:
// e = sqrt((a^2 - b^2)/a^2)
// where e is the eccentricity, a is the major axis length, b is the minor
// axis length.
double Oval::findEccentricity(void) const {
  return std::sqrt((pow(majorLength, 2) - pow(minorLength, 2))
		   /pow(majorLength, 2));
}

//##################################################################//
//######################### Circle #################################//

Circle::Circle(CoordPnt usrCenter, double usrRadius) :
  Oval(usrCenter, usrRadius, usrRadius) {}

Circle::Circle(CoordPnt usrCenter, double usrRadius, int numCoordPnts,
	 int layer, int datatype) :
  Oval(usrCenter, usrRadius, usrRadius, numCoordPnts, layer, datatype) {}    


//##################################################################//
//######################## Rectangle ###############################//

// The usual constructor for this class
Rectangle::Rectangle(CoordPnt usrCenter, double width, double height,
		     int layer, int datatype) {
  this->center = usrCenter;
  double maxX = usrCenter.getX() + width/2;
  double minX = usrCenter.getX() - width/2;
  double maxY = usrCenter.getY() + height/2;
  double minY = usrCenter.getY() - height/2;
  this->vertices.push_back(CoordPnt(minX, maxY)); // upper left corner
  this->vertices.push_back(CoordPnt(maxX, maxY)); // upper right corner
  this->vertices.push_back(CoordPnt(maxX, minY)); // lower right corner
  this->vertices.push_back(CoordPnt(minX, minY)); // lower left corner
  this->boundingBox = vertices; // for rectangles these are the same
  this->setLayer(layer);
  this->setDataType(datatype);
}

//##################################################################//
//########################## Square ################################//

// Use the Rectangle constructor but assume make sure both sides have the 
// same length
Square::Square(CoordPnt usrCenter, double sideLength, int layer, int datatype) : 
  Rectangle(usrCenter, sideLength, sideLength, layer, datatype) {};


//##################################################################//
//######################### Extern C ###############################//

// actually implement the python wrapping functions.
extern "C" {

  CoordPnt* CoordPnt_new(double usrX, double usrY) { return new CoordPnt(usrX, usrY); }

  void CoordPnt_delete(CoordPnt* coord) { delete coord; }

  LineSeg* LineSeg_new(CoordPnt* startSeg, CoordPnt* endSeg) {
    return new LineSeg(*startSeg, *endSeg);
  }

  void LineSeg_delete(LineSeg* lin) { delete lin; }
  
  //  Polygon* Polygon_new(void) { return new Polygon(); }

  Polygon* Polygon_new(double* xcoords, double* ycoords, int numcoords,
		       int layer, int datatype) {
    std::vector<CoordPnt> vertices;
    for (int i = 0; i < numcoords; i++)
      vertices.push_back(CoordPnt(xcoords[i], ycoords[i]));
    return new Polygon(vertices, layer, datatype);
  }

  void Polygon_delete(Polygon* poly) { delete poly; }
  
  Square* Square_new(CoordPnt* center, double sidelength, int layer,
		     int datatype) {
    return new Square(*center, sidelength, layer, datatype);
  }

  void Square_delete(Square* sq) { delete sq; }
  
  Rectangle* Rectangle_new(CoordPnt* center, double width,
			   double height, int layer, int datatype) {
    return new Rectangle(*center, width, height, layer, datatype);
  }

  void Rectangle_delete(Rectangle* rect) { delete rect; }

  Circle* Circle_new(CoordPnt* center, double radius, int numVertices,
		     int layer, int datatype) {
    return new Circle(*center, radius, numVertices, layer, datatype);
  }

  void Circle_delete(Circle* circ) { delete circ; }

  Oval* Oval_new(CoordPnt* center, double majorAxisLength, double minorAxisLength,
		 int numCoordPnts, int layer, int datatype) {
    return new Oval(*center, majorAxisLength, minorAxisLength,
		    numCoordPnts, layer, datatype);
  }

  void Oval_delete(Oval* ov) { delete ov; }

  Path* Path_new(double* xcoords, double* ycoords, int numcoords,
		 double width, int layer, int datatype) {
    std::vector<CoordPnt> vertices;
    for (int i = 0; i < numcoords; i++)
      vertices.push_back(CoordPnt(xcoords[i], ycoords[i]));
    return new Path(vertices, width, layer, datatype);
  }

  void Path_delete(Path* path) { delete path; }
  
}
