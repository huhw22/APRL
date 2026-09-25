#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>

namespace
{
  struct ParticleRecord
  {
    double position[3];
    double properVelocity[3];
    std::uint64_t sourceId;
    double macroWeight;
  };

  void requireHandle(hid_t handle, const std::string& message)
  {
    if (handle < 0) throw std::runtime_error(message);
  }

  void requireStatus(herr_t status, const std::string& message)
  {
    if (status < 0) throw std::runtime_error(message);
  }

  double lengthScale(const std::string& unit)
  {
    if (unit == "m" || unit == "meter" || unit == "metre") return 1.0;
    if (unit == "mm" || unit == "millimeter" || unit == "millimetre")
      return 1.0e-3;
    if (unit == "um" || unit == "micrometer" || unit == "micrometre")
      return 1.0e-6;
    if (unit == "nm" || unit == "nanometer" || unit == "nanometre")
      return 1.0e-9;
    throw std::invalid_argument("Unsupported --length-unit: " + unit);
  }

  bool parseRecord(const std::string& original, std::size_t lineNumber,
                   double scale, ParticleRecord& record)
  {
    std::string line = original;
    const std::string::size_type comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    if (line.find_first_not_of(" \t\r\n") == std::string::npos) return false;

    std::istringstream input(line);
    if (!(input >> record.position[0] >> record.position[1] >>
                  record.position[2] >> record.properVelocity[0] >>
                  record.properVelocity[1] >> record.properVelocity[2]))
      {
        std::ostringstream message;
        message << "Line " << lineNumber
                << " must contain x_plane y_plane zeta ux uy uz";
        throw std::runtime_error(message.str());
      }
    record.macroWeight = 1.0;
    if (!(input >> record.macroWeight))
      {
        if (!input.eof())
          {
            std::ostringstream message;
            message << "Line " << lineNumber
                    << " has an invalid optional macro_weight";
            throw std::runtime_error(message.str());
          }
        input.clear();
      }
    std::string trailing;
    if (input >> trailing)
      {
        std::ostringstream message;
        message << "Line " << lineNumber
                << " contains data after x_plane y_plane zeta ux uy uz [macro_weight]";
        throw std::runtime_error(message.str());
      }
    if (!(record.macroWeight > 0.0) ||
        !std::isfinite(record.macroWeight))
      throw std::runtime_error(
        "Line " + std::to_string(lineNumber) +
        " macro_weight must be positive and finite");
    for (unsigned int axis = 0; axis < 3; ++axis)
      {
        if (!std::isfinite(record.position[axis]) ||
            !std::isfinite(record.properVelocity[axis]))
          {
            std::ostringstream message;
            message << "Line " << lineNumber << " contains a non-finite value";
            throw std::runtime_error(message.str());
          }
        record.position[axis] *= scale;
      }
    return true;
  }

  hid_t memoryRecordType()
  {
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(ParticleRecord));
    requireHandle(type, "Cannot create HDF5 memory record type");
    hsize_t dimensions[1] = {3};
    hid_t vector = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dimensions);
    requireHandle(vector, "Cannot create HDF5 memory vector type");
    requireStatus(H5Tinsert(type, "position_m",
                            HOFFSET(ParticleRecord, position), vector),
                  "Cannot add position_m to HDF5 memory record");
    requireStatus(H5Tinsert(type, "proper_velocity",
                            HOFFSET(ParticleRecord, properVelocity), vector),
                  "Cannot add proper_velocity to HDF5 memory record");
    requireStatus(H5Tinsert(type, "source_id",
                            HOFFSET(ParticleRecord, sourceId),
                            H5T_NATIVE_UINT64),
                  "Cannot add source_id to HDF5 memory record");
    requireStatus(H5Tinsert(type, "macro_weight",
                            HOFFSET(ParticleRecord, macroWeight),
                            H5T_NATIVE_DOUBLE),
                  "Cannot add macro_weight to HDF5 memory record");
    H5Tclose(vector);
    return type;
  }

  hid_t fileRecordType()
  {
    const std::size_t f64 = 8;
    hid_t type = H5Tcreate(H5T_COMPOUND, 8 * f64);
    requireHandle(type, "Cannot create HDF5 file record type");
    hsize_t dimensions[1] = {3};
    hid_t vector = H5Tarray_create2(H5T_IEEE_F64LE, 1, dimensions);
    requireHandle(vector, "Cannot create HDF5 file vector type");
    requireStatus(H5Tinsert(type, "position_m", 0, vector),
                  "Cannot add position_m to HDF5 file record");
    requireStatus(H5Tinsert(type, "proper_velocity", 3 * f64, vector),
                  "Cannot add proper_velocity to HDF5 file record");
    requireStatus(H5Tinsert(type, "source_id", 6 * f64, H5T_STD_U64LE),
                  "Cannot add source_id to HDF5 file record");
    requireStatus(H5Tinsert(type, "macro_weight", 7 * f64, H5T_IEEE_F64LE),
                  "Cannot add macro_weight to HDF5 file record");
    H5Tclose(vector);
    return type;
  }

  void writeIntAttribute(hid_t object, const char* name, int value)
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create HDF5 attribute space");
    hid_t attribute = H5Acreate2(object, name, H5T_STD_I32LE, space,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute ") + name);
    requireStatus(H5Awrite(attribute, H5T_NATIVE_INT, &value),
                  std::string("Cannot write attribute ") + name);
    H5Aclose(attribute);
  }

  void writeStringAttribute(hid_t object, const char* name,
                            const std::string& value)
  {
    hid_t type = H5Tcopy(H5T_C_S1);
    requireHandle(type, "Cannot create HDF5 string type");
    requireStatus(H5Tset_size(type, value.size() + 1),
                  "Cannot size HDF5 string type");
    requireStatus(H5Tset_strpad(type, H5T_STR_NULLTERM),
                  "Cannot set HDF5 string padding");
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create HDF5 string attribute space");
    hid_t attribute = H5Acreate2(object, name, type, space,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute ") + name);
    requireStatus(H5Awrite(attribute, type, value.c_str()),
                  std::string("Cannot write attribute ") + name);
    H5Aclose(attribute);
    H5Tclose(type);
  }

  std::size_t countRecords(const std::string& filename, double scale)
  {
    std::ifstream input(filename.c_str());
    if (!input) throw std::runtime_error("Cannot open input file: " + filename);
    std::size_t count = 0;
    std::size_t lineNumber = 0;
    std::string line;
    ParticleRecord record = {};
    while (std::getline(input, line))
      {
        ++lineNumber;
        if (parseRecord(line, lineNumber, scale, record)) ++count;
      }
    if (count == 0)
      throw std::runtime_error("The text particle file contains no records");
    return count;
  }

  void convert(const std::string& inputFilename,
               const std::string& outputFilename, double scale)
  {
    const std::size_t total = countRecords(inputFilename, scale);
    if (total > static_cast<std::size_t>(
          std::numeric_limits<hsize_t>::max()))
      throw std::overflow_error("Particle count exceeds HDF5 size range");

    hid_t file = H5Fcreate(outputFilename.c_str(), H5F_ACC_TRUNC,
                           H5P_DEFAULT, H5P_DEFAULT);
    requireHandle(file, "Cannot create output file: " + outputFilename);
    hid_t group = H5Gcreate2(file, "/particles", H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    requireHandle(group, "Cannot create /particles group");
    writeIntAttribute(group, "format_version", 3);
    writeStringAttribute(group, "coordinate_frame",
                         "fixed-lab-plane-plus-longitudinal-offset");
    writeStringAttribute(group, "transverse_sampling",
                         "x-y-at-fixed-lab-input-plane");
    writeStringAttribute(group, "longitudinal_sampling",
                         "signed-offset-from-reconstructed-bunch-reference");
    writeStringAttribute(group, "position_unit", "m");
    writeStringAttribute(group, "proper_velocity_unit", "gamma*v/c");
    writeStringAttribute(group, "macro_weight_definition",
                         "positive relative weight normalized by simulator to beam.input.electrons");
    writeStringAttribute(group, "source_id_definition",
                         "one-based input data row");

    const hsize_t dimensions[1] = {static_cast<hsize_t>(total)};
    hid_t space = H5Screate_simple(1, dimensions, NULL);
    requireHandle(space, "Cannot create particle dataset space");
    hid_t fileType = fileRecordType();
    hid_t dataset = H5Dcreate2(group, "records", fileType, space,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    H5Tclose(fileType);
    requireHandle(dataset, "Cannot create /particles/records dataset");
    hid_t memoryType = memoryRecordType();

    const std::size_t chunkRecords = 65536;
    std::vector<ParticleRecord> buffer;
    buffer.reserve(std::min(chunkRecords, total));
    std::ifstream input(inputFilename.c_str());
    if (!input) throw std::runtime_error("Cannot reopen input file");
    std::size_t written = 0;
    std::size_t lineNumber = 0;
    std::string line;
    while (std::getline(input, line))
      {
        ++lineNumber;
        ParticleRecord record = {};
        if (!parseRecord(line, lineNumber, scale, record)) continue;
        record.sourceId = static_cast<std::uint64_t>(written + buffer.size() + 1);
        buffer.push_back(record);
        if (buffer.size() == chunkRecords)
          {
            hid_t fileSpace = H5Dget_space(dataset);
            const hsize_t start[1] = {static_cast<hsize_t>(written)};
            const hsize_t count[1] = {static_cast<hsize_t>(buffer.size())};
            requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
                                              start, NULL, count, NULL),
                          "Cannot select output particle hyperslab");
            hid_t memorySpace = H5Screate_simple(1, count, NULL);
            requireStatus(H5Dwrite(dataset, memoryType, memorySpace,
                                   fileSpace, H5P_DEFAULT, &buffer[0]),
                          "Cannot write output particle records");
            H5Sclose(memorySpace);
            H5Sclose(fileSpace);
            written += buffer.size();
            buffer.clear();
          }
      }
    if (!buffer.empty())
      {
        hid_t fileSpace = H5Dget_space(dataset);
        const hsize_t start[1] = {static_cast<hsize_t>(written)};
        const hsize_t count[1] = {static_cast<hsize_t>(buffer.size())};
        requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
                                          start, NULL, count, NULL),
                      "Cannot select final output particle hyperslab");
        hid_t memorySpace = H5Screate_simple(1, count, NULL);
        requireStatus(H5Dwrite(dataset, memoryType, memorySpace, fileSpace,
                               H5P_DEFAULT, &buffer[0]),
                      "Cannot write final output particle records");
        H5Sclose(memorySpace);
        H5Sclose(fileSpace);
        written += buffer.size();
      }
    H5Tclose(memoryType);
    H5Dclose(dataset);
    H5Gclose(group);
    requireStatus(H5Fclose(file), "Cannot close output HDF5 file");
    if (written != total)
      throw std::runtime_error("Particle input changed during conversion");
    std::cout << "Wrote " << written << " particle records to "
              << outputFilename << std::endl;
  }

  void usage(const char* executable)
  {
    std::cout << "Usage: " << executable
              << " --input particles.txt --output particles.h5 "
                 "[--length-unit m|mm|um|nm]\n"
              << "Each data line must contain: "
                 "x_plane y_plane zeta ux uy uz [macro_weight]\n";
  }
}

int main(int argc, char** argv)
{
  std::string input;
  std::string output;
  std::string unit("m");
  for (int index = 1; index < argc; ++index)
    {
      const std::string argument(argv[index]);
      if (argument == "--help" || argument == "-h")
        {
          usage(argv[0]);
          return 0;
        }
      if (index + 1 >= argc)
        {
          std::cerr << "Missing value after " << argument << std::endl;
          return 1;
        }
      if (argument == "--input") input = argv[++index];
      else if (argument == "--output") output = argv[++index];
      else if (argument == "--length-unit") unit = argv[++index];
      else
        {
          std::cerr << "Unknown option: " << argument << std::endl;
          return 1;
        }
    }
  if (input.empty() || output.empty())
    {
      usage(argv[0]);
      return 1;
    }
  try
    {
      convert(input, output, lengthScale(unit));
    }
  catch (const std::exception& error)
    {
      std::cerr << "Conversion failed: " << error.what() << std::endl;
      return 1;
    }
  return 0;
}
