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

#include "gdsfile.hpp"
#include "containers.hpp"

/// Prevent users from accidentally using utility methods that they should
/// not normally be using by nesting it within an obvious nested namespace.
namespace utils {

  // Basic bare bones constructor for this class.
  GDS_File::GDS_File(std::string usrFilename) :
    outputFile(usrFilename.c_str(), std::ios::out | std::ios::ate | std::ios::binary) {
    this->filename = usrFilename;
    this->version = 0x0258; // version 600 aka 6.0
    this->libraryName = "MyLibrary";
    this->generations = 1; // keep only the last version of each structure
    this->format = 0; // by default be a archive formated file
    this->databaseUnits = 1e-3; // databaseUnits are one thousandth of a user units
    this->userUnits = 1e-9*this->databaseUnits; // one nanometer
    this->refLib1.resize(44, '\0');
    this->refLib2.resize(44, '\0');
    this->sysLittleEndian = this->sysIsLittleEndian();
  }

  // Writes the whole GDSII file from a cell instance.
  // NOTE: each record must have an even number of bytes
  void GDS_File::WriteCell(const Cell* cell) {

    // Labels the content as belonging to this cell
    this->WriteStructureHeaderRecords(cell);

    // Itteratively write the contents of each polygon element
    std::vector<Polygon>::const_iterator firstPoly = cell->getPolygonList().begin();
    std::vector<Polygon>::const_iterator lastPoly = cell->getPolygonList().end();
    for (std::vector<Polygon>::const_iterator polygon = firstPoly; 
	 polygon != lastPoly; ++polygon) {
      this->WriteElementHeaderRecords(BOUNDARY); // polygons are boundary typed
      this->WriteElementContentRecords(polygon);
      this->WriteElementTailRecords();
    }

    std::vector<Path>::const_iterator firstPath = cell->getPathList().begin();
    std::vector<Path>::const_iterator lastPath = cell->getPathList().end();
    for (std::vector<Path>::const_iterator path = firstPath;
	 path != lastPath; ++path) {
      this->WriteElementHeaderRecords(PATH);
      this->WriteElementContentRecords(path);
      this->WriteElementTailRecords();
    }

    std::vector<CellReference>::const_iterator firstRef = cell->getCellReferenceList().begin();
    std::vector<CellReference>::const_iterator lastRef = cell->getCellReferenceList().end();
    for (std::vector<CellReference>::const_iterator cellRef = firstRef;
	 cellRef != lastRef; ++cellRef) {
      this->WriteElementHeaderRecords(SREF);
      this->WriteElementContentRecords(cellRef);
      this->WriteFileTailRecords();
    }

    std::vector<CellArray>::const_iterator firstArray = cell->getCellArrayList().begin();
    std::vector<CellArray>::const_iterator lastArray = cell->getCellArrayList().end();
    for (std::vector<CellArray>::const_iterator cellArray = firstArray;
	 cellArray != lastArray; ++cellArray) {
      this->WriteElementHeaderRecords(AREF);
      this->WriteElementContentRecords(cellArray);
      this->WriteFileTailRecords();
    }

    // Write the ENDSTR that corresponds to this cell
    this->WriteStructureTailRecords();

  }

  void GDS_File::Write(const std::vector<Cell*> cellVec) {
      
    this->WriteFileHeaderRecords();
      
    for (uint i = 0; i < cellVec.size(); i++) 
      this->WriteCell(cellVec[i]);

    this->WriteFileTailRecords();
  }

  void GDS_File::WriteFileHeaderRecords() {
    // the size of the label, the in16_t that details the size of the record -
    // including itself - as well as the part off the label that details 
    // which record this is (i.e. HEADER, etc)
    const int16_t RECORD_LABEL_SIZE = 2*sizeof(int16_t); 
    std::string recordStr; // for recording the string data of a record
    int16_t recordSize; // total size of record
    std::vector<int16_t> recordInt; // contains the int data for the record

    //----------------------------------------------------------------------//
    // HEADER
    // Write the beginning of library record
    recordSize = sizeof(version) + RECORD_LABEL_SIZE; 
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(HEADER);
    this->writeInt16ToFile(this->version);

    //----------------------------------------------------------------------//
    // BGNLIB
    // output the date last modified
    time_t now = time(0);
    tm *ltm = localtime(&now);
    int16_t year =  1900 + ltm->tm_year;
    int16_t month = 1 + ltm->tm_mon;
    int16_t day =  ltm->tm_mday;
    int16_t hour =  1 + ltm->tm_hour;
    int16_t minute = 1 + ltm->tm_min;
    int16_t second = 1 + ltm->tm_sec;
    // mark that we are looking at BGNLIB
    // Factor comes from: 12 for the date/time
    recordSize = 12*sizeof(int16_t) + RECORD_LABEL_SIZE; 
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(BGNLIB);
    // mark the time last modified (now)
    this->writeInt16ToFile(year);
    this->writeInt16ToFile(month);
    this->writeInt16ToFile(day);
    this->writeInt16ToFile(hour);
    this->writeInt16ToFile(minute);
    this->writeInt16ToFile(second);

    // mark the time it was last accessed (now, since it was just created)
    this->writeInt16ToFile(year);
    this->writeInt16ToFile(month);
    this->writeInt16ToFile(day);
    this->writeInt16ToFile(hour);
    this->writeInt16ToFile(minute);
    this->writeInt16ToFile(second);
      
    //----------------------------------------------------------------------//
    // LIBNAME
    int16_t vecLength = this->libraryName.size();
    if (vecLength % 2 != 0) // every record must be even in number of bytes
      vecLength++;
    char libName[vecLength]; // must have terminating zero and be evenly sized
    for (int i = 0; i < vecLength; i++)
      libName[i] = this->libraryName[i];
    if (vecLength % 2 != 0)
      libName[vecLength - 1] = ' ';
    recordSize = sizeof(libName) + RECORD_LABEL_SIZE;
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(LIBNAME);
    this->writeCharToFile(libName, vecLength);

    //----------------------------------------------------------------------//
    // UNITS
    recordSize = 0x0014;
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(UNITS);
    this->writeFloat64ToFile(this->databaseUnits);
    this->writeFloat64ToFile(this->userUnits);

  } // GDS_File::Write

    // concludes each gdsii file
  void GDS_File::WriteFileTailRecords() {
    int16_t recordSize = sizeof(ENDLIB) + sizeof(int16_t); 
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(ENDLIB);
  } // WriteFileTailRecords

  void GDS_File::WriteStructureTailRecords() {
    int16_t recordSize = sizeof(ENDLIB) + sizeof(int16_t); 
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(ENDSTR);
  }

  void GDS_File::WriteStructureHeaderRecords(const Cell* cell) {
    const int16_t RECORD_LABEL_SIZE = 2*sizeof(int16_t); 
    int16_t recordSize;
    //----------------------------------------------------------------------//
    // BGNSTR
    time_t now = time(0);
    tm *ltm = localtime(&now);
    int16_t year =  1900 + ltm->tm_year;
    int16_t month =  1 + ltm->tm_mon;
    int16_t day =  ltm->tm_mday;
    int16_t hour =  1 + ltm->tm_hour;
    int16_t minute = 1 + ltm->tm_min;
    int16_t second = 1 + ltm->tm_sec;
    // 12 for 2X year, month etc, and 2 for the recordSize itself and BGNSTR
    recordSize = 12*sizeof(int16_t) + RECORD_LABEL_SIZE;
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(BGNSTR);
    // creation date
    this->writeInt16ToFile(year);
    this->writeInt16ToFile(month);
    this->writeInt16ToFile(day);
    this->writeInt16ToFile(hour);
    this->writeInt16ToFile(minute);
    this->writeInt16ToFile(second);

    this->writeInt16ToFile(year);
    this->writeInt16ToFile(month);
    this->writeInt16ToFile(day);
    this->writeInt16ToFile(hour);
    this->writeInt16ToFile(minute);
    this->writeInt16ToFile(second);
    //----------------------------------------------------------------------//
    // STRNAME
    int strSize = cell->getCellname().size();
    if (strSize % 2 != 0)
      strSize++;
    char cellnameChar[strSize];
    std::string cellname = cell->getCellname();
    int orgStrSize = cell->getCellname().size();
    for (int i = 0; i < orgStrSize; i++)
      cellnameChar[i] = cellname[i];
    if (orgStrSize % 2 != 0)
      cellnameChar[strSize - 1] = ' ';
    recordSize = strSize*sizeof(char) + sizeof(STRNAME) + sizeof(int16_t);
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(STRNAME);
    this->writeCharToFile(cellnameChar, strSize);
  }

  void GDS_File::WriteElementHeaderRecords(int16_t dataType) {
    // Each possible record is exactly four bytes long as none of the possible
    // records hold any data
    int16_t recordSize = 0x0004;
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(dataType);
  }

  void GDS_File::WriteElementContentRecords(std::vector<Polygon>::const_iterator polygon) {
    const int16_t RECORD_LABEL_SIZE = 2*sizeof(int16_t); 
    int16_t recordSize;
    // We already have wrote that we are in a "BOUNDARY" element
    recordSize = 0x0006;
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(LAYER);
    int16_t currentLayer = polygon->getLayer();
    this->writeInt16ToFile(currentLayer);
    // Each layer must be followed by a datatype
    recordSize = 0x0006;
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(DATATYPE);
    int16_t currentDataType = polygon->getDataType();
    this->writeInt16ToFile(currentDataType);
    // Now record each (x, y) coordinate pair
    std::vector<CoordPnt> myVertices = polygon->getVertices();
    recordSize = 2*(myVertices.size() + 1)*sizeof(int32_t) 
      + RECORD_LABEL_SIZE;
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(XY);
    for (std::vector<CoordPnt>::iterator xy = myVertices.begin(); 
	 xy != myVertices.end(); ++xy) {
      int32_t curX = xy->getX()/this->databaseUnits;
      this->writeInt32ToFile(curX);
      int32_t curY = xy->getY()/this->databaseUnits;
      this->writeInt32ToFile(curY);
    }
    // rerecord the first one as is GDS2 standard (marks the end of 
    // a polygon)
    int32_t firstX = myVertices[0].getX()/this->databaseUnits;
    this->writeInt32ToFile(firstX);
    int32_t firstY = myVertices[0].getY()/this->databaseUnits;
    this->writeInt32ToFile(firstY);
  }

  /// \brief Overridden for use with path elements
  void GDS_File::WriteElementContentRecords(std::vector<Path>::const_iterator path) {
    // Path
    const int16_t RECORD_LABEL_SIZE = 2*sizeof(int16_t); 
    int16_t recordSize;
    // -- Layer
    recordSize = 0x0006;
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(LAYER);
    int16_t currentLayer = path->getLayer();
    this->writeInt16ToFile(currentLayer);
    // -- Data Type
    // Each layer must be followed by a datatype
    recordSize = 0x0006;
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(DATATYPE);
    int16_t currentDataType = path->getDataType();
    this->writeInt16ToFile(currentDataType);
    // -- Path Type
    recordSize = 0x0006;
    int16_t pathtype = path->getPathType();
    // only need to record path type if it is not zero as zero is 
    // assumed if this record does not exist.
    if (pathtype != 0) {
      this->writeInt16ToFile(recordSize);
      this->writeInt16ToFile(PATHTYPE);
      this->writeInt16ToFile(pathtype);
    }
    // -- Width
    recordSize = 0x0008;
    int32_t width = path->getPathWidth();
    // only need to record path width if it is not zero as zero is 
    // assumed if this record does not exist.
    if (width != 0) {
      this->writeInt16ToFile(recordSize);
      this->writeInt16ToFile(WIDTH);
      this->writeInt32ToFile(width);
    }
    // -- XY
    // Now record each (x, y) coordinate pair
    std::vector<CoordPnt> myVertices = path->getCoordPath();
    recordSize = 2*(myVertices.size() + 1)*sizeof(int32_t) 
      + RECORD_LABEL_SIZE;
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(XY);
    for (std::vector<CoordPnt>::iterator xy = myVertices.begin(); 
	 xy != myVertices.end(); ++xy) {
      int32_t curX = xy->getX()/this->databaseUnits;
      this->writeInt32ToFile(curX);
      int32_t curY = xy->getY()/this->databaseUnits;
      this->writeInt32ToFile(curY);
    }
  }

  void GDS_File::WriteElementContentRecords(std::vector<CellReference>::const_iterator cellRef) {
    int16_t recordSize;
    const int16_t RECORD_LABEL_SIZE = 2*sizeof(int16_t); 

    // -- SREF
    recordSize = RECORD_LABEL_SIZE;
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(SREF);

    // -- SNAME
    int strSize = cellRef->getCellname().size();
    if (strSize % 2 != 0)
      strSize++;
    char cellnameChar[strSize];
    std::string cellname = cellRef->getCellname();
    int orgStrSize = cellRef->getCellname().size();
    for (int i = 0; i < orgStrSize; i++)
      cellnameChar[i] = cellname[i];
    if (orgStrSize % 2 != 0)
      cellnameChar[strSize - 1] = ' ';
    recordSize = strSize*sizeof(char) + sizeof(STRNAME) + sizeof(int16_t);
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(STRNAME);
    this->writeCharToFile(cellnameChar, strSize);      

    // -- strans
    char stransRecord[15];
    recordSize = RECORD_LABEL_SIZE + sizeof(stransRecord);
    this->writeInt16ToFile(recordSize);

    this->writeInt16ToFile(STRANS);

    // only certain parts of this record will not be zero, so make
    // sure that zero is the "default."
    for (uint i = 0; i < sizeof(stransRecord); i++)
      stransRecord[i] = 0;
    // setting bit 13 flags the reference has having non unity
    // manification
    if (cellRef->getMagnification() != 1.0)
      stransRecord[12] = 1;
    // Setting bit 14 flags the reference has having nonzero 
    // rotation
    if (cellRef->getRotation() != 0.0)
      stransRecord[13] = 1;
    this->writeCharToFile(stransRecord, sizeof(stransRecord));

    // If no MAG record is there then the magnification is assumed
    // to be one. Therefore, only write MAG if the magnification is
    // not one.
    if (cellRef->getMagnification() != 1.0) {
      recordSize = RECORD_LABEL_SIZE + sizeof(float64);
      this->writeInt16ToFile(recordSize);
      this->writeInt16ToFile(MAG);
      this->writeFloat64ToFile(cellRef->getMagnification());
    }

    // If no ANGLE record is present then the angle of rotation is 
    // assumed to be zero (note that there can still be a reflection
    // accross the x-axis which is in effect having a rotation of 
    // 180 degrees but I will keep all rotation information in ANGLE)
    if (cellRef->getRotation() != 0.0) {
      recordSize = RECORD_LABEL_SIZE + sizeof(float64);
      this->writeInt16ToFile(recordSize);
      this->writeInt16ToFile(ANGLE);
      this->writeFloat64ToFile(cellRef->getRotation());	
    }      

    // -- XY
    recordSize = RECORD_LABEL_SIZE + 2*sizeof(int32_t);
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(XY);
    int32_t curX = cellRef->getCenter().getX()/this->databaseUnits;
    this->writeInt32ToFile(curX);
    int32_t curY = cellRef->getCenter().getY()/this->databaseUnits;
    this->writeInt32ToFile(curY);
  }

  void GDS_File::WriteElementContentRecords(std::vector<CellArray>::const_iterator cellArray) {
    int16_t recordSize;
    const int16_t RECORD_LABEL_SIZE = 2*sizeof(int16_t); 

    // -- AREF
    recordSize = RECORD_LABEL_SIZE;
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(AREF);

    // -- SNAME
    int strSize = cellArray->getCellname().size();
    if (strSize % 2 != 0)
      strSize++;
    char cellnameChar[strSize];
    std::string cellname = cellArray->getCellname();
    int orgStrSize = cellArray->getCellname().size();
    for (int i = 0; i < orgStrSize; i++)
      cellnameChar[i] = cellname[i];
    if (orgStrSize % 2 != 0)
      cellnameChar[strSize - 1] = ' ';
    recordSize = strSize*sizeof(char) + sizeof(STRNAME) + sizeof(int16_t);
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(STRNAME);
    this->writeCharToFile(cellnameChar, strSize);      

    // -- strans
    char stransRecord[15];
    recordSize = RECORD_LABEL_SIZE + sizeof(stransRecord);
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(STRANS);
    // only certain parts of this record will not be zero, so make
    // sure that zero is the "default."
    for (uint i = 0; i < sizeof(stransRecord); i++)
      stransRecord[i] = 0;
    // setting bit 13 flags the reference has having non unity
    // manification
    if (cellArray->getMagnification() != 1.0)
      stransRecord[12] = 1;
    // Setting bit 14 flags the reference has having nonzero 
    // rotation
    if (cellArray->getRotation() != 0.0)
      stransRecord[13] = 1;
    this->writeCharToFile(stransRecord, sizeof(stransRecord));

    // If no MAG record is there then the magnification is assumed
    // to be one. Therefore, only write MAG if the magnification is
    // not one.
    if (cellArray->getMagnification() != 1.0) {
      recordSize = RECORD_LABEL_SIZE + sizeof(float64);
      this->writeInt16ToFile(recordSize);
      this->writeInt16ToFile(MAG);
      this->writeFloat64ToFile(cellArray->getMagnification());
    }

    // If no ANGLE record is present then the angle of rotation is 
    // assumed to be zero (note that there can still be a reflection
    // accross the x-axis which is in effect having a rotation of 
    // 180 degrees but I will keep all rotation information in ANGLE)
    if (cellArray->getRotation() != 0.0) {
      recordSize = RECORD_LABEL_SIZE + sizeof(float64);
      this->writeInt16ToFile(recordSize);
      this->writeInt16ToFile(ANGLE);
      this->writeFloat64ToFile(cellArray->getRotation());	
    }      

    // -- COLROW
    recordSize = RECORD_LABEL_SIZE + 2*sizeof(int16_t);
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(cellArray->getNumCol());
    this->writeInt16ToFile(cellArray->getNumRow());

    // -- XY
    // FIXME!!!!!! NEEDS MORE XY DATUMS
    recordSize = RECORD_LABEL_SIZE + 6*sizeof(int32_t);
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(XY);
    int32_t curX = cellArray->getStartingPos().getX()/this->databaseUnits;
    this->writeInt32ToFile(curX);
    int32_t curY = cellArray->getStartingPos().getY()/this->databaseUnits;
    this->writeInt32ToFile(curY);

    // record the furthest column
    int32_t extraDistance =  
      cellArray->getXSpacing()*cellArray->getNumCol()/this->databaseUnits;
    this->writeInt32ToFile(curX + extraDistance);
    this->writeInt32ToFile(curY);

    // record the furthest row
    this->writeInt32ToFile(curX);
    extraDistance = 
      cellArray->getYSpacing()*cellArray->getNumRow()/this->databaseUnits;
    this->writeInt32ToFile(curY + extraDistance);
  }

  void GDS_File::WriteElementTailRecords() {
    int16_t recordSize = 0x0004; // with no data it is exactly 4 bytes long
    this->writeInt16ToFile(recordSize);
    this->writeInt16ToFile(ENDEL);
  }

  void GDS_File::writeInt16ToFile(int16_t data) {
    // gdsII file format is written in big endian
    if (this->sysLittleEndian)
      data = int16Swap(data);
    this->outputFile.write(reinterpret_cast<char*>(&data), sizeof(int16_t));      
  }

  void GDS_File::writeInt32ToFile(int32_t data) {
    // gdsII file format is written in big endian
    if (this->sysLittleEndian)
      data = int32Swap(data);
    this->outputFile.write(reinterpret_cast<char*>(&data), sizeof(int32_t));      
  }

  void GDS_File::writeFloat64ToFile(float64 data) {
    int numBytes = 8; // 64 bit => 8 bytes with 8 bits each
    int numBitsPerByte = 8;
    int endOfByte = 7;
    char output[numBytes];
    // set all bits to 0
    for (int i = 0; i < numBytes; i++)
      output[i] = 0x00; 

    // use bitwise operations
    //
    // set bit x in number (set means make it be 1)
    // number |= 1 << x
    //
    // clear bit x in number
    // number &= ~(1 << x)
    //
    // toggle bit x in number
    // number ^= 1 << x
    //
    // check bit x
    // bit = (number >> x) & 1
  
    // leftmost bit specifies negative (1) or positive (0)
    if (data < 0) {
      output[0] |= 1 << endOfByte;
      data = -data;
    }

    int exponent;
    double mantissa;
    double error = 1e-9;
    if (data < 1e-77) // 16^{-64} smallest number possible
      output[0] = 0x00; // we only have touched first byte so only need to reset first byte
    else {
      exponent = (int) log(data)/log(16.0);
      mantissa = data/pow(16.0, exponent);
      // may be off by one power of 16
      if (mantissa < 1.0/16.0) {
	mantissa *= 16.0;
	exponent--;
      } else if (mantissa >= 1.0) {
	mantissa /= 16.0;
	exponent++;
      }

      if (!(mantissa < 1 && mantissa >= 1.0/16.0)) {
	std::stringstream errorMsg;
	errorMsg << "Mantissa must be inbetween 1 and 1/16. "
		 << "It is " << mantissa << "\n";
	throw std::logic_error(errorMsg.str());
      }
      if (!(exponent >= -64 && exponent < 63)) {
	std::stringstream errorMsg;
	errorMsg << "Exponent must be equal to or greater than "
		 << " 64, and less than 64. It is " << exponent
		 << "\n";
	throw std::logic_error(errorMsg.str());
      }

      // we need to have an exponential offset
      exponent += 64;

    }

    int numOfExponentBits = 7;
    for (int i = 1; i <= numOfExponentBits; i++) {
      if (exponent >= 128/pow(2, i)) {
	output[0] |= 1 << endOfByte - i;
	exponent -= 128/pow(2, i);
      }
    }


    int bitNum = 1;
    for (int byteNum = 1; byteNum <= numBytes; byteNum++) {
      for (int bitOffset = 0; bitOffset < numBitsPerByte; bitOffset++) {
	double bitValue = 1.0/pow(2, bitNum);
	bitNum++;
	if (mantissa >= bitValue) {
	  output[byteNum] |= 1 << endOfByte - bitOffset;
	  mantissa -= bitValue;
	}
      }
    }

    this->outputFile.write(reinterpret_cast<char*>(&output),
			   numBytes*sizeof(char));      
  }

  void GDS_File::writeCharToFile(char data[], int size) {
    // there is no difference between little endian and big endian data
    // for char data, as it is only one byte long thus there cannot be
    // byte swapping.
    this->outputFile.write(data, size);      
  }

  bool GDS_File::sysIsLittleEndian() {
    uint16_t num = 1;
    return (*(char *)&num == 1);
  }

  // switches the endianess of the input (num)
  int16_t int16Swap(int16_t num) {
    unsigned char byte1, byte2;
    byte1 = num & 255;
    byte2 = (num >> 8) & 255;
    return (byte1 << 8) + byte2;
  }

  // switches the endianess of the input (num)
  int32_t int32Swap(int32_t num) {
    unsigned char byte1, byte2, byte3, byte4;
    byte1 = num & 255;
    byte2 = (num >> 8) & 255;
    byte3 = (num >> 16) & 255;
    byte4 = (num >> 24) & 255;
    return ((int32_t) byte1 << 24) + ((int32_t) byte2 << 16) + 
      ((int32_t) byte3 << 8) + byte4;
  }

} // namespace utils

