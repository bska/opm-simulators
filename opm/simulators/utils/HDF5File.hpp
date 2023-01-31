/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef HDF5_FILE_HPP
#define HDF5_FILE_HPP

#include <hdf5.h>

#include <string>
#include <vector>

namespace Opm {

//! \brief Class handling simple output to HDF5.
class HDF5File {
public:
    //! \brief Enumeration of file opening modes.
    enum class OpenMode {
        APPEND,    //!< Append to an existing file (creates new if not)
        OVERWRITE, //!< Overwrite and write to an existing file
        READ       //!< Open existing file for reading
    };

    //! \brief Opens HDF5 file for I/O.
    //! \param fileName Name of file to open
    //! \param mode Open mode for file
    HDF5File(const std::string& fileName, OpenMode mode);

    //! \brief Destructor clears up any opened files.
    ~HDF5File();

    //! \brief Write a char buffer to a specified location in file.
    //! \param group Group ("directory") to write data to
    //! \param dset Data set ("file") to write data to
    //! \param buffer Data to write
    //! \details Throws exception on failure
    void write(const std::string& group,
               const std::string& dset,
               const std::vector<char>& buffer);

    //! \brief Read a char buffer from a specified location in file.
    //! \param group Group ("directory") to read data from
    //! \param dset Data set ("file") to read data from
    //! \param buffer Vector to store read data in
    //! \details Throws exception on failure
    void read(const std::string& group,
              const std::string& dset,
              std::vector<char>& buffer);

private:
    hid_t m_file = H5I_INVALID_HID; //!< File handle
};

}

#endif
