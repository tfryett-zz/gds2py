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

#include "containers.hpp"

//##################################################################//
//##################### CellReference ##############################//

CellReference::CellReference(Cell& referenceCell, int centerX, 
			     int centerY) : refCell(referenceCell) {
  this->center = CoordPnt(centerX, centerY);
  this->magnification = 1.0;
  this->rotation = 0.0;
}

CellReference::CellReference(Cell& referenceCell) :
  refCell(referenceCell) {
  this->center = CoordPnt(0, 0);
  this->magnification = 1.0;
  this->rotation = 0.0;
}

CellReference::CellReference(Cell& referenceCell, CoordPnt centerPnt) :
  refCell(referenceCell) {
  this->center = centerPnt;
  this->magnification = 1.0;
  this->rotation = 0.0;
}

CellReference::CellReference(Cell& referenceCell, CoordPnt centerPnt,
			     double magnification, double rotation) :
  refCell(referenceCell) {
  this->center = centerPnt;
  this->setMagneification(magnification);
  this->setRotation(rotation);
}

std::vector<CoordPnt> CellReference::findVertices(void) {
  std::vector<CoordPnt> boundingVertex;
  std::vector<Polygon> polygonVec = this->refCell.getPolygonList();
  double totMax = std::numeric_limits<double>::max();
  double totMin = std::numeric_limits<double>::min();
  double minX = totMax;
  double minY = totMax;
  double maxX = totMin;
  double maxY = totMin;
  // find the minimum and maximum values of x and y coordinates
  for(std::vector<Polygon>::iterator polyIt = polygonVec.begin(); 
      polyIt != polygonVec.end(); ++polyIt) {
    std::vector<CoordPnt> vertexVec = polyIt->getVertices();
    for (std::vector<CoordPnt>::iterator vertIt = vertexVec.begin();
	 vertIt != vertexVec.end(); ++vertIt) {
      if (vertIt->getX() > maxX)
	maxX = vertIt->getX();
      if (vertIt->getX() < minX)
	minX = vertIt->getX();
      if (vertIt->getY() > maxY)
	maxY = vertIt->getY();
      if (vertIt->getY() < minY)
	minY = vertIt->getY();
    }
  }
  // check to make sure they actually changed into something
  if (minX == totMax && minY == totMax && 
      maxX == totMin && maxY == totMin)
    throw std::logic_error("Invalid polygons detected.");
  else if ((minX == maxX) && (minY == maxY))
    throw std::logic_error("Invalid polygons detected.");

  boundingVertex.push_back(CoordPnt(minX, minY));
  boundingVertex.push_back(CoordPnt(minX, maxY));
  boundingVertex.push_back(CoordPnt(maxX, minY));
  boundingVertex.push_back(CoordPnt(maxX, maxY));
  return boundingVertex;
}

double CellReference::getMagnification(void) const {
  return this->magnification;
}

void CellReference::setMagneification(double newMagnification) {
  this->magnification = newMagnification;
}

double CellReference::getRotation(void) const {
  return this->rotation;
}

void CellReference::setRotation(double newRotation) {
  this->rotation = newRotation;
}

CoordPnt CellReference::getCenter(void) const {
  return this->center;
}

void CellReference::setCenter(CoordPnt newCenter) {
  this->center = newCenter;
}

std::string CellReference::getCellname() const {
  return this->refCell.getCellname();
}

//##################################################################//
//######################### CellArray ##############################//

CellArray::CellArray(Cell& cell, CoordPnt usrStartingPos, 
		     int usrNumCol, int usrNumRow, double usrXSpacing, 
		     double usrYSpacing, double magnification,
		     double rotation) : refCell(cell) {
  this->startingPos = usrStartingPos;
  this->setNumCol(usrNumCol);
  this->setNumRow(usrNumRow);
  this->xSpacing = usrXSpacing;
  this->ySpacing = usrYSpacing;
  this->setMagnification(magnification);
  this->setRotation(rotation);
}

void CellArray::setStartingPos(CoordPnt newStartingPos) {
  this->startingPos = newStartingPos;
}

CoordPnt CellArray::getStartingPos() const {
  return this->startingPos;
}

void CellArray::setXSpacing(double newXSpacing) {
  this->xSpacing = newXSpacing;
}
    
double CellArray::getXSpacing(void) const {
  return this->xSpacing;
}

void CellArray::setYSpacing(double newYSpacing) {
  this->ySpacing = newYSpacing;
}
    
double CellArray::getYSpacing(void) const {
  return this->ySpacing;
}

void CellArray::setNumCol(int newNumCol) {
  // the magic number 32,767 is from the GDSII standard. No reason
  // was given, but we will adhere to the rule.
  if (newNumCol > 0 && newNumCol <= 32767)
    this->numCol = newNumCol;
  else {
    std::stringstream errorMsg;
    errorMsg << "The number of columns may not exceed 32,767, and"
	     << " must be nonzero. User entered " << newNumCol 
	     << ".\n";
    std::invalid_argument(errorMsg.str());
  }
}

int CellArray::getNumCol() const {
  return this->numCol;
}

void CellArray::setNumRow(int newNumRow) {
  // the magic number 32,767 is from the GDSII standard. No reason
  // was given, but we will adhere to the rule.
  if (newNumRow > 0 && newNumRow <= 32767)
    this->numRow = newNumRow;
  else {
    std::stringstream errorMsg;
    errorMsg << "The number of rows may not exceed 32,767, and"
	     << " must be nonzero. User entered " << newNumRow 
	     << ".\n";
    std::invalid_argument(errorMsg.str());
  }  
}

int CellArray::getNumRow() const {
  return this->numRow;
}

double CellArray::getRotation() const {
  return this->rotation;
}
    
void CellArray::setRotation(double newRotation) {
  this->rotation = newRotation;
}

double CellArray::getMagnification() const {
  return this->magnification;
}

void CellArray::setMagnification(double newMagnification) {
  if (newMagnification > 0)
    this->numRow = newMagnification;
  else {
    std::stringstream errorMsg;
    errorMsg << "The number of rows may not exceed 32,767, and"
	     << " must be nonzero. User entered " << newMagnification 
	     << ".\n";
    std::invalid_argument(errorMsg.str());
  }  
}

std::string CellArray::getCellname() const {
  return this->refCell.getCellname();
}


//##################################################################//
//########################## Cell ##################################//

// The full blooded Cell object constructor.
Cell::Cell(std::string usrCellname) {
  this->setCellname(usrCellname);
  time_t now = time(0);
  tm *ltm = localtime(&now);
  this->timeCreated.year = 1900 + ltm->tm_year;
  this->timeCreated.month =  1 + ltm->tm_mon;
  this->timeCreated.day =  ltm->tm_mday;
  this->timeCreated.hour =  1 + ltm->tm_hour;
  this->timeCreated.minute = 1 + ltm->tm_min;
  this->timeCreated.second = 1 + ltm->tm_sec;
}

// Allows the user to reset this cell's name
void Cell::setCellname(std::string usrCellname) {
  if (usrCellname.size() > 32)
    throw std::invalid_argument("Cell name must be under 32 characters.");
  if (usrCellname.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_?$") !=
      std::string::npos)
    throw std::invalid_argument("Cell name contains invalid characters.");
  this->cellname = usrCellname;
}

// Simple member function that returns the specified cell's name
std::string Cell::getCellname() const {
  return this->cellname;
}

std::vector<Polygon>& Cell::getPolygonList() const {
  return const_cast<std::vector<Polygon> &> (this->polyList);
}
  
// Returns the data about when this cell was created 
timeData Cell::getTimeData() const {   
  return this->timeCreated;
}

void Cell::addPolygon(Polygon usrPolygon) {
  this->polyList.push_back(usrPolygon);
}

void Cell::addPath(Path usrPath) {
  this->pathList.push_back(usrPath);
}

void Cell::addCellReference(CellReference usrCellReference) {
  this->cellReferenceList.push_back(usrCellReference);
}

void Cell::addCellArray(CellArray usrCellArray) {
  this->cellArrayList.push_back(usrCellArray);
}

std::vector<Path>& Cell::getPathList() const {
  return const_cast<std::vector<Path> &> (this->pathList);
}

std::vector<CellReference>& Cell::getCellReferenceList() const {
  return const_cast<std::vector<CellReference> &> (this->cellReferenceList);
}

std::vector<CellArray>& Cell::getCellArrayList(void) const {
  return const_cast<std::vector<CellArray> &> (this->cellArrayList);
}

//##################################################################//
//######################### Layout #################################//

Layout::Layout() {}

void Layout::addCell(Cell& usrCell) {
  this->cellVec.push_back(&usrCell);
}

std::vector<Cell*> Layout::getCells() const {
  return this->cellVec;
}

void Layout::write(std::string usrFilename) {
  utils::GDS_File myFile(usrFilename);
  myFile.Write(this->cellVec);
}


//##################################################################//
//######################## Extern C ################################//


extern "C" {
  
  Layout* Layout_new(void) { return new Layout(); }
  
  void Layout_delete(Layout* lay) { delete lay; }
  
  void Layout_addCell(Layout* layout, Cell* cell) { layout->addCell(*cell); }

  void Layout_write(Layout* layout, const char* filename) {
    layout->write(std::string(filename, strlen(filename)));
  }

  Cell* Cell_new(const char* cellname) {
    return new Cell(std::string(cellname, strlen(cellname)));
  }

  void Cell_delete(Cell* cell) { delete cell; }

  void Cell_addPolygon(Cell* cell, Polygon* usrPolygon) { cell->addPolygon(*usrPolygon); }

  void Cell_addPath(Cell* cell, Path* usrPath) { cell->addPath(*usrPath); }

  void Cell_addCellReference(Cell* cell, CellReference* usrCellReference) {
    cell->addCellReference(*usrCellReference);
  }
  
  void Cell_addCellArray(Cell* cell, CellArray* usrCellArray) {
    cell->addCellArray(*usrCellArray);
  }

  CellArray* CellArray_new(Cell* cell, CoordPnt* usrStartingPos, int usrNumCol, 
			   int usrNumRow, double xSpacing, double ySpacing,
			   double magnification, double rotation) {
    return new CellArray(*cell, *usrStartingPos, usrNumCol, 
			 usrNumRow, xSpacing, ySpacing, magnification,
			 rotation);
  }

  void CellArray_delete(CellArray* cellArray) { delete cellArray; }


  CellReference* CellReference_new(Cell* referenceCell, CoordPnt* center,
				   double magnification, double rotation) {
    return new CellReference(*referenceCell, *center, magnification, rotation);
  }

  void CellReference_delete(CellReference* referenceCell) { delete referenceCell; }

}
