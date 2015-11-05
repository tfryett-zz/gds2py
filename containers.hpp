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

#ifndef CONTAINERS_HPP
#define CONTAINERS_HPP

#include <limits>
#include <vector>
#include <string>
#include <ctime>
#include <functional>
#include "shapes.hpp"
#include "gdsfile.hpp"

/// \brief A struct to contain all of the date and time information needed by GDSII standards.
struct timeData {
  int year; 
  int month; 
  int day; 
  int hour; 
  int minute; 
  int second; 
};


//##################################################################//
//##################### CellReference ##############################//

class Cell;

class CellReference {
private:
  const Cell& refCell;
  std::vector<CoordPnt> findVertices(void);
  CoordPnt center;

protected:
  double magnification; //!< The amount the referenced cell will be magnified compared to the origianl.
  double rotation; //!< The ammount the reference cell will be rotated
  std::vector<CoordPnt> vertices;

public:
  CellReference(Cell& referenceCell);
  CellReference(Cell& referenceCell, CoordPnt centerPnt);
  CellReference(Cell& referenceCell, int centerX, int centerY);
  CellReference(Cell& referenceCell, CoordPnt centerPnt,
		double magnification, double rotation);

  double getMagnification(void) const;

  void setMagneification(double newMagnification);

  double getRotation(void) const;

  void setRotation(double newRotation);

  CoordPnt getCenter(void) const;

  void setCenter(CoordPnt newCenter);

  std::string getCellname(void) const;

};


//##################################################################//
//######################### CellArray ##############################//

/// class CellArray
///
/// \brief The GDSII represenation of an array of referenced cells.
///
/// This class is used primarily to save memory and to reduce human
/// error. The concept is very similar to constant references in
/// C and C++ in that it is a low memory "stand in" that does not
/// allow for direct modification of the referenced cell, but does
/// become updated when the instance of the referenced cell is
/// changed.
class CellArray {
private:
  Cell& refCell; //!< The actual reference that makes it a cell reference.
  int numCol; //!< The number of times @refCell will be written in the x direction.
  int numRow; //!< The number of times @refCell will be written in the y direction.
  double xSpacing; //!< The spacing between each refCell in the x direction.
  double ySpacing; //!< The spacing between each refCell in the y direction.
  CoordPnt startingPos; //!< The coordinate of the lower left @refCell of the array of @refCell
  double magnification; //!< The magnification factor to apply to each @refCell
  double rotation; //!< The Rotation to apply to each of @refCell.

protected:

public:
  /// \brief The sole constructor for the CellArray class
  ///
  /// @cell The reference to the cell we want to reference.
  /// @usrStartingPos The coordinate of the lower left part of the cell array.
  /// @usrNumCol The number of columns in the cell array.
  /// @usrNumRow The number of rows in the cell array.
  /// @xSpacing The spacing between each column.
  /// @ySpacing The spacing between each row.
  /// @magnification The magnification of the referenced cell in the array.
  /// @rotation The rotation of the referenced cell in the array.
  CellArray(Cell& cell, CoordPnt usrStartingPos, int usrNumCol, 
	    int usrNumRow, double xSpacing, double ySpacing,
	    double magnification = 1, double rotation = 0);

  /// \brief Sets the position of the lower left corner to @newStartingPos.
  ///
  /// @newStartingPos the new position of the lower left corner of the array.
  void setStartingPos(CoordPnt newStartingPos);

  /// \brief Returns the position of the lower left corner of the array.
  CoordPnt getStartingPos(void) const;

  /// \brief Sets the spacing between array columns.
  ///
  /// @newXSpacing the new spacing between array columns.
  ///
  /// The x spacing is allowed to be negative, which will cause the
  /// array to propogate along the negative x direction.
  void setXSpacing(double newXSpacing);
    
  /// \brief Returns the spacing inbetween columns.
  double getXSpacing(void) const;

  /// \brief Sets the spacing between array rows.
  ///
  /// @newXSpacing the new spacing between array rows.
  ///
  /// The x spacing is allowed to be negative, which will cause the
  /// array to propogate along the negative y direction.
  void setYSpacing(double newYSpacing);
    
  /// \brief Returns the spacing inbetween rows.
  double getYSpacing(void) const;

  /// \brief Sets a new spacing inbetween columns.
  ///
  /// @newNumCol The number of columns to be set.
  ///
  /// First checks if @newNumCol is a valid number of columns 
  /// (that is @newNumCol > 0) and then sets the member @numCol 
  /// or if it is not a valid number of columns throws an
  /// std::invalid_argument exception.
  void setNumCol(int newNumCol);

  /// \brief Returns the number of columns in the array.
  int getNumCol(void) const;

  /// \brief Sets a new spacing inbetween rows.
  ///
  /// @newNumRow The number of rows to be set.
  ///
  /// First checks if @newNumRow is a valid number of rows 
  /// (that is @newNumRow > 0) and then sets the member @numRow 
  /// or if it is not a valid number of rows throws an
  /// std::invalid_argument exception
  void setNumRow(int newNumRow);

  /// \brief Returns the number of rows in the array.
  int getNumRow(void) const;

  /// \brief Returns the name of the referenced cell.
  std::string getCellname(void) const;

  /// \brief Returns the value of rotation of the array members (in radians).
  double getRotation(void) const;
    
  /// \brief Sets a value for the rotation of the array members (in radians.)
  ///
  /// @newRotation The rotation of the array members (in radians).
  void setRotation(double newRotation);

  /// \brief Returns the value of magnification of each array element.
  ///
  /// A magnification of 1.0 means no magnification.
  double getMagnification(void) const;

  /// \brief Sets the value of magnification of each array element to be @newMagnification.
  ///
  /// @newMagnification The magnification to be set.
  ///
  /// A magnification of 1.0 means no magnification. Only positive
  /// magnification is allowed, and is only limted by the capacity
  /// of an eight byte double.
  void setMagnification(double newMagnification);

};


//##################################################################//
//########################## Cell ##################################//

/// Cells that are added to a Layout can be considered to be a
/// component of a larger fabrication flow. An example would be
/// the design for a transistor that is to be repeatedly used.
/// Having it be defined as a Cell would allow for a CellReference
/// to take the place of each transitor. This not only saves
/// space but also acts to counter mistakes in trying to define
/// many copies of the same thing.
class Cell {
private:

protected:
  std::string cellname; //!< The name this object.
  std::vector<Polygon> polyList; //<! The vector of objects the cell contains.
  std::vector<Path> pathList; //!< The vector of path objects in the cell.
  std::vector<CellReference> cellReferenceList; //!< The vector of CellReference objects that this cell contains.
  std::vector<CellArray> cellArrayList; //!< The vector of CellArray objects that this cell contains.
  timeData timeCreated; //!< The struct containing a year, month, ... second.

public:

  /// \brief Creates a Cell object with the specified cellname.
  Cell(std::string usrCellname);
    
  /// \brief Sets the cell name to be that of what the user specifies.
  ///
  /// @usrCellname The name of this cell.
  void setCellname(std::string usrCellname);

  /// \brief Returns the cell name for this object.
  std::string getCellname(void) const;

  /// \brief Returns a const reference to the Entity List
  std::vector<Polygon>& getPolygonList(void) const;

  /// \brief Returns a vector full of the info about when the cell was 
  /// created
  ///
  /// The elements of the returned vector are members of a struct
  /// of type timeData with fields as follows:
  /// Year
  /// Month
  /// Day
  /// Hour
  /// Minute
  /// Second
  timeData getTimeData(void) const;

  /// \brief Adds a polygon to the Cell.
  ///
  /// @usrPolygon The polygon to add to the Cell.
  void addPolygon(Polygon usrPolygon);

  /// \brief Adds a path to the Cell.
  ///
  /// @usrPath The path to add to the Cell.    
  void addPath(Path usrPath);

  /// \brief Returns the vector of Path objects in this Cell;
  std::vector<Path>& getPathList(void) const;

  /// \brief Adds a cell reference to the vector of cellReferences;
  ///
  /// @usrCellReference The CellReference object to add to this Cell.
  void addCellReference(CellReference usrCellReference);

  /// \brief Adds a cell reference array to the vector of cellArrays.
  ///
  /// @usrCellArray The CellArray object to add to this Cell.
  void addCellArray(CellArray usrCellArray);

  /// \brief Returns the vector of CellReferences contained in this Cell.
  std::vector<CellReference>& getCellReferenceList(void) const;

  /// \brief Returns the vector of CellArrays contained in this Cell.
  std::vector<CellArray>& getCellArrayList(void) const;
}; // class Cell


//##################################################################//
//######################### Layout #################################//

/// class Layout
///
/// This class is the overall container for sil. That is if you
/// wish to import from or write to a GDSII file you must do so
/// through this class. Think of this class as a library of Cell
/// objects which constitute the fabriaction flow you are wishing
/// to create.
class Layout {
private:

  std::vector<Cell*> cellVec; //!< \brief The collection of Cell pointers which constitute a Layout.

protected:

public:
  /// \brief The default constructor for this class.
  Layout(void);

  /// \brief The method used to add a Cell object to this Layout.
  ///
  /// @usrCell The Cell object you wish to add to this Layout.
  ///
  /// Cells that are added to a Layout can be considered to be a
  /// component of a larger fabrication flow. An example would be
  /// the design for a transistor that is to be repeatedly used.
  /// Having it be defined as a Cell would allow for a CellReference
  /// to take the place of each transitor. This not only saves
  /// space but also acts to counter mistakes in trying to define
  /// many copies of the same thing.
  void addCell(Cell& usrCell);

  /// \brief Writes all of the contained Cell objects to a file.
  ///
  /// @filename The name of the file to write to.
  void write(std::string filename);

  /// \brief Returns all of the Cell objects that are contained.
  std::vector<Cell*> getCells(void) const;

};

//##################################################################//
//######################## Extern C ################################//

extern "C" {
  
  Layout* Layout_new(void);
  
  void Layout_delete(Layout* lay);

  void Layout_addCell(Layout* layout, Cell* cell);

  void Layout_write(Layout* layout, const char* filename);

  Cell* Cell_new(const char* cellname);

  void Cell_delete(Cell* cell);

  void Cell_addPolygon(Cell* cell, Polygon* usrPolygon);

  void Cell_addPath(Cell* cell, Path* usrPath);

  void Cell_addCellReference(Cell* cell, CellReference* usrCellReference);

  void Cell_addCellArray(Cell* cell, CellArray* usrCellArray);
  
  CellArray* CellArray_new(Cell* cell, CoordPnt* usrStartingPos, int usrNumCol, 
			   int usrNumRow, double xSpacing, double ySpacing,
			   double magnification, double rotation);
 
  void CellArray_delete(CellArray* cellArray);

  CellReference* CellReference_new(Cell* referenceCell, CoordPnt* center,
				   double mangnification, double rotation);

  void CellReference_delete(CellReference* referenceCell);

}

#endif // CONTAINERS_HPP
