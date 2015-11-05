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

#ifndef GDS_FILE_HXX
#define GDS_FILE_HXX

typedef double float64;

#include <fstream>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h> // cross-compiler integer datatypes
#include <stdexcept>
#include "shapes.hpp"

class Cell;
class CellReference;
class CellArray;

/// Prevent users from accidentally using utility methods that they should
/// not normally be using by nesting it within an obvious nested namespace.
namespace utils {

  struct Time {
    int16_t year;
    int16_t month;
    int16_t day;
    int16_t hour;
    int16_t second;
  };


  // Switches the input from little endian to big endian or vice versa.
  int16_t int16Swap(int16_t num);

  // Switches the input from little endian to big endian or vice versa.
  int32_t int32Swap(int32_t num);

  // Switches the input from little endian to big endian or vice versa.
  float64 float64Swap(float64 num);

  // These are the record markers needed for reading and writing to gds2 files
  // The first byte specifies the "record" and the second specifies the type
  // of data contained in the record.
  // E.g. for HEADER: 0x00 specifies header that this is a header record and 
  // 0x02 specifies that the record contains unsigned two byte integer data.
  const int16_t HEADER       = 0x0002; 
  const int16_t BGNLIB       = 0x0102;
  const int16_t LIBNAME      = 0x0206;
  const int16_t UNITS        = 0x0305;
  const int16_t ENDLIB       = 0x0400;
  const int16_t BGNSTR       = 0x0502;
  const int16_t STRNAME      = 0x0606;
  const int16_t ENDSTR       = 0x0700;
  const int16_t BOUNDARY     = 0x0800;
  const int16_t PATH         = 0x0900;
  const int16_t SREF         = 0x0A00;
  const int16_t AREF         = 0x0B00;
  const int16_t TEXT         = 0x0C00;
  const int16_t LAYER        = 0x0D02;
  const int16_t DATATYPE     = 0x0E02;
  const int16_t WIDTH        = 0x0F03;
  const int16_t XY           = 0x1003;
  const int16_t ENDEL        = 0x1100;
  const int16_t SNAME        = 0x1206;
  const int16_t COLROW       = 0x1302;
  const int16_t NODE         = 0x1500;
  const int16_t TEXTTYPE     = 0x1602; 
  const int16_t PRESENTATION = 0x1701;
  const int16_t STRING       = 0x1906;
  const int16_t STRANS       = 0x1A01;
  const int16_t MAG          = 0x1B05;
  const int16_t ANGLE        = 0x1C05;
  const int16_t REFLIBS      = 0x1F06;
  const int16_t FONTS        = 0x2006;
  const int16_t PATHTYPE     = 0x2102;
  const int16_t GENERATIONS  = 0x2202;
  const int16_t ATTRTABLE    = 0x2306;
  const int16_t ELFLAGS      = 0x2601;
  const int16_t NODETYPE     = 0x2A02;
  const int16_t PROPATTR     = 0x2B02;
  const int16_t PROPVALUE    = 0x2C06;
  const int16_t BOX          = 0x2D00;
  const int16_t BOXTYPE      = 0x2E02;
  const int16_t PLEX         = 0x2F02;
  const int16_t TAPENUM      = 0x3203;
  const int16_t TAPECODE     = 0x3302;
  const int16_t FORMAT       = 0x3602;
  const int16_t MASK         = 0x3706;
  const int16_t ENDMASKS     = 0x3800;

  class GDS_File {
  private:
    std::string filename; //!< The name of the output file.
    int16_t version; //!< The version of GDSII we will use
    std::string libraryName; //!< the Library name of the gdsII file
    std::string refLib1; //!< The first of two possible reference libraries. Does not have to be used. 
    std::string refLib2; //!< The second of two possible reference libraries. Does not have to be used. 
    std::vector<std::string> fonts; //!< List of fonts to use in gds structures
    int16_t generations; //!< The number of generation of a structure that are retained.
    int16_t format; //!< Specifies whether this is a Archive (0) or Filtered (1) formatted file.
    float64 databaseUnits; //!< Relative size of units stored in the database to the user's defined units.
    float64 userUnits; //!< Size of a unit in meters.
    std::ofstream outputFile; //!< Reference to the iostream to the output file.
    bool sysLittleEndian; //!< Is the global variable to 
    Time timeCreated;

    /// \brief Writes the data at the top of the GDSII file that specifies
    /// pertanent information about the rest of the file.
    ///
    /// Here are the exact contents of the File Header Record:
    ///
    /// Conent Name:    Hex Code:     Type of Data:
    /// HEADER	      0002          2-byte integer
    /// BGNLIB	      0102          12 2-byte integers
    /// LIBNAME	      0206          ASCII string
    /// REFLIBS         1F06          2 45-character ASCII strings
    /// FONTS	          2006          4 44-character ASCII strings
    /// ATTRTABLE	      2306          44-character ASCII string
    /// GENERATION      2202          2-byte integer
    /// FORMAT	  3602          2-byte integer
    /// UNITS	          0305          2 8-byte floats
    /// MASK	          3706          ASCII string
    /// ENDMASKS	  3800          No data
    void WriteFileHeaderRecords(void);

    /// \brief Writes the data at the end of the GDSII file that specifies
    /// that the file stream has ended.
    ///
    /// Here are the exact contents of the File Tail Record:
    ///
    /// Conent Name:    Hex Code:     Type of Data:
    /// ENDLIB	  0400   	    No data
    void WriteFileTailRecords(void);

    /// \brief Writes the Header that is required at the start of each
    /// Structure.
    ///
    /// Here are the exact contents of the Structure Header Record:
    ///
    /// Conent Name:    Hex Code:     Type of Data:
    /// BGNSTR	  0502          12 2-byte integers
    /// STRNAME	  0606          Up to 32-characters ASCII string
    void WriteStructureHeaderRecords(const Cell* cell);

    /// \brief Writes the conents that are required at the end of the each 
    /// Structure.
    ///
    /// Here are the exact contents of the Structure Tail Record:
    ///
    /// Conent Name:    Hex Code:     Type of Data:
    /// ENDSTR	  0700          No data
    void WriteStructureTailRecords(void);

    /// \brief Writes the conents that specify a new element in the GDSII file
    /// stream.
    ///
    /// @dataType is one of the listed "Content Name" variables which will 
    ///           correspond to the data in the element record.
    ///
    /// Here are the exact contents of the Element Header Record:
    ///
    /// Conent Name:    Hex Code:     Type of Data:
    /// BOUNDARY	  0800          No data
    /// PATH	          0900          No data
    /// SREF	          0A00          No data
    /// AREF	          0B00          No data
    /// TEXT	          0C00          No data
    /// NODE	          1500          No data
    /// BOX             2D00          No data
    void WriteElementHeaderRecords(int16_t dataType);

    /// \brief Writes the Element content which is the core data of a GDSII
    /// file - i.e. this is where the polygons are defined.
    ///
    /// Here are the possible contents of the Element Content Record:
    ///
    /// Conent Name:    Hex Code:     Type of Data:
    /// ELFLAGS	  2601          2-byte integer
    /// PLEX	          2F03          4-byte integer
    /// LAYER	          0D02          2-byte integers
    /// DATATYPE	  0E02          2-byte integer
    /// XY	          1003          Up to 200 4-byte integer pairs
    /// PATHTYPE	  2102          2-byte integer
    /// WIDTH	          0F03          4-byte integer
    /// SNAME	          1206          Up to 32-character ASCII string
    /// STRANS	  1A01          2-byte integer
    /// MAG	          1B05          8-byte float
    /// ANGLE	          1C05          8-byte float
    /// COLROW	  1302          2 2-byte integer
    /// TEXTTYPE	  1602          2-byte integer
    /// PRESENTATION    1701          2-byte integer
    /// ASCII STRING	  1906          Up to 512-character string
    /// NODETYPE	  2A02          2-byte integer
    /// BOXTYPE	  2E02          2-byte integer
    void WriteElementContentRecords(std::vector<Polygon>::const_iterator polygon);

    /// \brief Overridden for use with Path elements
    void WriteElementContentRecords(std::vector<Path>::const_iterator path);

    /// \brief Overridden for use with CellReference elements
    void WriteElementContentRecords(std::vector<CellReference>::const_iterator cellRef);

    /// \brief Overridden for use with CellArray elements
    void WriteElementContentRecords(std::vector<CellArray>::const_iterator cellArray);

    /// \brief Marks the end of each element record.
    ///
    /// Conent Name:    Hex Code:     Type of Data:
    /// ENDEL           08000         No data
    void WriteElementTailRecords(void);

    // Writes input data to the file specified by outputFile in a size sensitive
    // and endian sensitive manner (gdsII are always big endian).
    void writeInt16ToFile(int16_t data);

    // Writes input data to the file specified by outputFile in a size sensitive
    // and endian sensitive manner (gdsII are always big endian).
    void writeCharToFile(char data[], int size);

    // Writes input data to the file specified by outputFile in a size sensitive
    // and endian sensitive manner (gdsII are always big endian).
    void writeInt32ToFile(int32_t data);

    // Writes input data to the file specified by outputFile in a size sensitive
    // and endian sensitive manner (gdsII are always big endian).
    void writeFloat64ToFile(float64 data);

    // Tests if the system is little endian or big endian.
    // returns true if the system is a little endian system.
    bool sysIsLittleEndian(void);


    /// \brief Causes the GDS_File object to write the specified Cell object
    /// to the file specified by the private field filename.
    void WriteCell(const Cell* cell);

  protected:

  public:
    /// \brief Creates a GDS_File object which allows users to write their 
    /// data to a file of their choosing.
    /// 
    /// @usrFilename The file which the data will be written to.
    ///
    /// This constructor only creates an internal object - it does NOT 
    /// create or modify any files until it is told to do so by calling
    /// the object's member functions (e.g. Write()).
    GDS_File(std::string usrFilename);
      
    void Write(const std::vector<Cell*> cellVec);

  }; // class GDS_FILE
} // namespace utils

#endif // GDS_FILE_HXX
