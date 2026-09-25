#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>
#include <SDDS.h>

namespace
{
  struct Options
  {
    std::string input;
    std::string output;
    std::string weightColumn;
    std::string idColumn;
    int page;

    Options() : idColumn("particleID"), page(1) {}
  };

  struct ParticleRecord
  {
    double planePosition[2];
    double arrivalTimeOffset;
    double properVelocity[3];
    std::uint64_t sourceId;
    double macroWeight;
  };

  class CompensatedSum
  {
  public:
    CompensatedSum() : sum_(0.0L), correction_(0.0L) {}

    void add(long double value)
    {
      const long double adjusted = value - correction_;
      const long double next = sum_ + adjusted;
      correction_ = (next - sum_) - adjusted;
      sum_ = next;
    }

    long double value() const { return sum_; }

  private:
    long double sum_;
    long double correction_;
  };

  void requireHandle(hid_t handle, const std::string& message)
  {
    if (handle < 0) throw std::runtime_error(message);
  }

  void requireStatus(herr_t status, const std::string& message)
  {
    if (status < 0) throw std::runtime_error(message);
  }

  std::string sddsErrors()
  {
    int32_t count = 0;
    char** messages = SDDS_GetErrorMessages(&count, SDDS_ALL_GetErrorMessages);
    std::string result;
    for (int32_t index = 0; index < count; ++index)
      {
        if (!messages[index]) continue;
        if (!result.empty()) result += "; ";
        result += messages[index];
      }
    if (messages)
      {
        SDDS_FreeStringArray(messages, count);
        SDDS_Free(messages);
      }
    return result;
  }

  void failSdds(const std::string& message)
  {
    const std::string details = sddsErrors();
    throw std::runtime_error(details.empty() ? message :
                             message + ": " + details);
  }

  class SddsInput
  {
  public:
    explicit SddsInput(const std::string& filename) : active_(false)
    {
      dataset_ = SDDS_DATASET();
      std::vector<char> path(filename.begin(), filename.end());
      path.push_back('\0');
      if (!SDDS_InitializeInput(&dataset_, &path[0]))
        failSdds("Cannot initialize SDDS input " + filename);
      active_ = true;
    }

    ~SddsInput()
    {
      if (active_ && !SDDS_Terminate(&dataset_))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    }

    SDDS_DATASET* get() { return &dataset_; }

  private:
    SddsInput(const SddsInput&);
    SddsInput& operator=(const SddsInput&);

    SDDS_DATASET dataset_;
    bool active_;
  };

  class NumericColumn
  {
  public:
    NumericColumn(SDDS_DATASET* dataset, const std::string& name)
      : data_(NULL), type_(0), name_(name)
    {
      std::vector<char> mutableName(name.begin(), name.end());
      mutableName.push_back('\0');
      const int32_t index = SDDS_GetColumnIndex(dataset, &mutableName[0]);
      if (index < 0)
        throw std::runtime_error("SDDS column is missing: " + name);
      type_ = SDDS_GetColumnType(dataset, index);
      if (!SDDS_NUMERIC_TYPE(type_) && type_ != SDDS_CHARACTER)
        throw std::runtime_error("SDDS column is not numeric: " + name);
      data_ = SDDS_GetInternalColumn(dataset, &mutableName[0]);
      if (!data_) failSdds("Cannot access SDDS column " + name);
    }

    double value(std::int64_t row) const
    {
      switch (type_)
        {
        case SDDS_LONGDOUBLE:
          return static_cast<double>(
            static_cast<const long double*>(data_)[row]);
        case SDDS_DOUBLE:
          return static_cast<const double*>(data_)[row];
        case SDDS_FLOAT:
          return static_cast<const float*>(data_)[row];
        case SDDS_LONG64:
          return static_cast<double>(
            static_cast<const std::int64_t*>(data_)[row]);
        case SDDS_ULONG64:
          return static_cast<double>(
            static_cast<const std::uint64_t*>(data_)[row]);
        case SDDS_LONG:
          return static_cast<const std::int32_t*>(data_)[row];
        case SDDS_ULONG:
          return static_cast<const std::uint32_t*>(data_)[row];
        case SDDS_SHORT:
          return static_cast<const short*>(data_)[row];
        case SDDS_USHORT:
          return static_cast<const unsigned short*>(data_)[row];
        case SDDS_CHARACTER:
          return static_cast<const char*>(data_)[row];
        default:
          throw std::runtime_error(
            "Unsupported numeric SDDS type in column " + name_);
        }
    }

    std::uint64_t unsignedInteger(std::int64_t row) const
    {
      switch (type_)
        {
        case SDDS_LONG64:
          {
            const std::int64_t value =
              static_cast<const std::int64_t*>(data_)[row];
            if (value < 0) invalidId(row);
            return static_cast<std::uint64_t>(value);
          }
        case SDDS_ULONG64:
          return static_cast<const std::uint64_t*>(data_)[row];
        case SDDS_LONG:
          {
            const std::int32_t value =
              static_cast<const std::int32_t*>(data_)[row];
            if (value < 0) invalidId(row);
            return static_cast<std::uint64_t>(value);
          }
        case SDDS_ULONG:
          return static_cast<const std::uint32_t*>(data_)[row];
        case SDDS_SHORT:
          {
            const short value = static_cast<const short*>(data_)[row];
            if (value < 0) invalidId(row);
            return static_cast<std::uint64_t>(value);
          }
        case SDDS_USHORT:
          return static_cast<const unsigned short*>(data_)[row];
        case SDDS_CHARACTER:
          return static_cast<unsigned char>(
            static_cast<const char*>(data_)[row]);
        default:
          {
            const long double value = static_cast<long double>(this->value(row));
            const long double rounded = std::floor(value + 0.5L);
            if (!std::isfinite(static_cast<double>(value)) || value < 0.0L ||
                value != rounded ||
                value > static_cast<long double>(UINT64_MAX))
              invalidId(row);
            /* A floating SDDS ID cannot carry all uint64 values exactly. */
            if (value > 9007199254740992.0L)
              throw std::runtime_error(
                "Floating-point SDDS particle ID exceeds the exact 2^53 range in column " +
                name_);
            return static_cast<std::uint64_t>(rounded);
          }
        }
    }

  private:
    void invalidId(std::int64_t row) const
    {
      throw std::runtime_error(
        "SDDS source ID must be a nonnegative integer in column " + name_ +
        " at selected-page row " + std::to_string(row + 1));
    }

    const void* data_;
    int32_t type_;
    std::string name_;
  };

  bool hasColumn(SDDS_DATASET* dataset, const std::string& name)
  {
    if (name.empty()) return false;
    std::vector<char> mutableName(name.begin(), name.end());
    mutableName.push_back('\0');
    return SDDS_GetColumnIndex(dataset, &mutableName[0]) >= 0;
  }

  void requireColumn(SDDS_DATASET* dataset, const std::string& name,
                     const char* units)
  {
    std::vector<char> mutableName(name.begin(), name.end());
    mutableName.push_back('\0');
    std::vector<char> mutableUnits;
    char* unitPointer = NULL;
    if (units)
      {
        mutableUnits.assign(units, units + std::char_traits<char>::length(units));
        mutableUnits.push_back('\0');
        unitPointer = &mutableUnits[0];
      }
    if (SDDS_CheckColumn(dataset, &mutableName[0], unitPointer,
                         SDDS_ANY_NUMERIC_TYPE, NULL) != SDDS_CHECK_OK)
      {
        std::string message = "Required numeric SDDS column " + name;
        if (units) message += " must have units " + std::string(units);
        throw std::runtime_error(message);
      }
  }

  void requireElegantColumns(SDDS_DATASET* dataset)
  {
    requireColumn(dataset, "x", "m");
    requireColumn(dataset, "xp", NULL);
    requireColumn(dataset, "y", "m");
    requireColumn(dataset, "yp", NULL);
    requireColumn(dataset, "t", "s");
    try
      {
        requireColumn(dataset, "p", "m$be$nc");
      }
    catch (const std::exception&)
      {
        requireColumn(dataset, "p", NULL);
        std::cerr << "WARNING: SDDS column p has no standard Elegant "
                     "m$be$nc unit declaration; interpreting it as beta*gamma."
                  << std::endl;
      }
  }

  hid_t memoryRecordType()
  {
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(ParticleRecord));
    requireHandle(type, "Cannot create HDF5 memory record type");
    hsize_t positionDimensions[1] = {2};
    hsize_t velocityDimensions[1] = {3};
    hid_t position = H5Tarray_create2(
      H5T_NATIVE_DOUBLE, 1, positionDimensions);
    hid_t velocity = H5Tarray_create2(
      H5T_NATIVE_DOUBLE, 1, velocityDimensions);
    requireHandle(position, "Cannot create HDF5 plane-position type");
    requireHandle(velocity, "Cannot create HDF5 proper-velocity type");
    requireStatus(H5Tinsert(type, "plane_position_m",
                            HOFFSET(ParticleRecord, planePosition), position),
                  "Cannot add plane_position_m to HDF5 memory type");
    requireStatus(H5Tinsert(type, "arrival_time_offset_s",
                            HOFFSET(ParticleRecord, arrivalTimeOffset),
                            H5T_NATIVE_DOUBLE),
                  "Cannot add arrival_time_offset_s to HDF5 memory type");
    requireStatus(H5Tinsert(type, "proper_velocity",
                            HOFFSET(ParticleRecord, properVelocity), velocity),
                  "Cannot add proper_velocity to HDF5 memory type");
    requireStatus(H5Tinsert(type, "source_id",
                            HOFFSET(ParticleRecord, sourceId),
                            H5T_NATIVE_UINT64),
                  "Cannot add source_id to HDF5 memory type");
    requireStatus(H5Tinsert(type, "macro_weight",
                            HOFFSET(ParticleRecord, macroWeight),
                            H5T_NATIVE_DOUBLE),
                  "Cannot add macro_weight to HDF5 memory type");
    H5Tclose(position);
    H5Tclose(velocity);
    return type;
  }

  hid_t fileRecordType()
  {
    const std::size_t f64 = 8;
    hid_t type = H5Tcreate(H5T_COMPOUND, 8 * f64);
    requireHandle(type, "Cannot create HDF5 file record type");
    hsize_t positionDimensions[1] = {2};
    hsize_t velocityDimensions[1] = {3};
    hid_t position = H5Tarray_create2(
      H5T_IEEE_F64LE, 1, positionDimensions);
    hid_t velocity = H5Tarray_create2(
      H5T_IEEE_F64LE, 1, velocityDimensions);
    requireHandle(position, "Cannot create HDF5 file plane-position type");
    requireHandle(velocity, "Cannot create HDF5 file proper-velocity type");
    requireStatus(H5Tinsert(type, "plane_position_m", 0, position),
                  "Cannot add plane_position_m to HDF5 file type");
    requireStatus(H5Tinsert(type, "arrival_time_offset_s", 2 * f64,
                            H5T_IEEE_F64LE),
                  "Cannot add arrival_time_offset_s to HDF5 file type");
    requireStatus(H5Tinsert(type, "proper_velocity", 3 * f64, velocity),
                  "Cannot add proper_velocity to HDF5 file type");
    requireStatus(H5Tinsert(type, "source_id", 6 * f64, H5T_STD_U64LE),
                  "Cannot add source_id to HDF5 file type");
    requireStatus(H5Tinsert(type, "macro_weight", 7 * f64,
                            H5T_IEEE_F64LE),
                  "Cannot add macro_weight to HDF5 file type");
    H5Tclose(position);
    H5Tclose(velocity);
    return type;
  }

  void writeIntAttribute(hid_t object, const char* name, int value)
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create HDF5 integer attribute space");
    hid_t attribute = H5Acreate2(object, name, H5T_STD_I32LE, space,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute ") + name);
    requireStatus(H5Awrite(attribute, H5T_NATIVE_INT, &value),
                  std::string("Cannot write attribute ") + name);
    H5Aclose(attribute);
  }

  void writeDoubleAttribute(hid_t object, const char* name, double value)
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create HDF5 double attribute space");
    hid_t attribute = H5Acreate2(object, name, H5T_IEEE_F64LE, space,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute ") + name);
    requireStatus(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &value),
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

  void writeChunk(hid_t dataset, hid_t memoryType,
                  const std::vector<ParticleRecord>& records,
                  std::uint64_t offset)
  {
    if (records.empty()) return;
    hid_t fileSpace = H5Dget_space(dataset);
    requireHandle(fileSpace, "Cannot access output particle dataspace");
    const hsize_t start[1] = {static_cast<hsize_t>(offset)};
    const hsize_t count[1] = {static_cast<hsize_t>(records.size())};
    requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
                                      start, NULL, count, NULL),
                  "Cannot select output particle hyperslab");
    hid_t memorySpace = H5Screate_simple(1, count, NULL);
    requireHandle(memorySpace, "Cannot create output particle memory space");
    const herr_t status = H5Dwrite(dataset, memoryType, memorySpace,
                                   fileSpace, H5P_DEFAULT, &records[0]);
    H5Sclose(memorySpace);
    H5Sclose(fileSpace);
    requireStatus(status, "Cannot write output particle records");
  }

  void convert(const Options& options)
  {
    SddsInput input(options.input);
    SDDS_DATASET* dataset = input.get();
    const int headerVersion = dataset->layout.version;
    if (headerVersion < 1 || headerVersion > SDDS_VERSION)
      throw std::runtime_error("Official SDDS reader returned an invalid protocol version");
    requireElegantColumns(dataset);
    if (!options.weightColumn.empty())
      requireColumn(dataset, options.weightColumn, NULL);

    int currentPage = 0;
    int32_t readStatus = 0;
    while (currentPage < options.page &&
           (readStatus = SDDS_ReadPage(dataset)) > 0)
      ++currentPage;
    if (readStatus <= 0 || currentPage != options.page)
      failSdds("Requested SDDS page does not exist");

    const std::int64_t rows = SDDS_CountRowsOfInterest(dataset);
    if (rows <= 0)
      throw std::runtime_error("Selected SDDS page contains no particles");
    if (rows != dataset->n_rows)
      throw std::runtime_error(
        "Unexpected SDDS row filtering state in selected page");
    if (static_cast<std::uint64_t>(rows) >
        static_cast<std::uint64_t>(std::numeric_limits<hsize_t>::max()))
      throw std::overflow_error("Particle count exceeds HDF5 size range");

    NumericColumn x(dataset, "x");
    NumericColumn xp(dataset, "xp");
    NumericColumn y(dataset, "y");
    NumericColumn yp(dataset, "yp");
    NumericColumn time(dataset, "t");
    NumericColumn momentum(dataset, "p");
    const bool hasWeight = !options.weightColumn.empty();
    const bool hasId = hasColumn(dataset, options.idColumn);
    NumericColumn* weight = hasWeight ?
      new NumericColumn(dataset, options.weightColumn) : NULL;
    NumericColumn* id = hasId ?
      new NumericColumn(dataset, options.idColumn) : NULL;

    CompensatedSum weightAccumulator;
    CompensatedSum weightedTimeAccumulator;
    double minimumP = std::numeric_limits<double>::infinity();
    double maximumP = 0.0;
    for (std::int64_t row = 0; row < rows; ++row)
      {
        const double rowWeight = hasWeight ? weight->value(row) : 1.0;
        const double rowTime = time.value(row);
        const double rowX = x.value(row);
        const double rowXp = xp.value(row);
        const double rowY = y.value(row);
        const double rowYp = yp.value(row);
        const double rowP = momentum.value(row);
        if (!(rowWeight > 0.0) || !std::isfinite(rowWeight))
          throw std::runtime_error(
            "SDDS macro weight must be positive and finite at selected-page row " +
            std::to_string(row + 1));
        if (!std::isfinite(rowTime))
          throw std::runtime_error(
            "SDDS arrival time is not finite at selected-page row " +
            std::to_string(row + 1));
        if (!std::isfinite(rowX) || !std::isfinite(rowXp) ||
            !std::isfinite(rowY) || !std::isfinite(rowYp) ||
            !std::isfinite(rowP) || !(rowP > 0.0))
          throw std::runtime_error(
            "Elegant phase-space coordinates must be finite and p must be positive at selected-page row " +
            std::to_string(row + 1));
        if (hasId) id->unsignedInteger(row);
        weightAccumulator.add(static_cast<long double>(rowWeight));
        weightedTimeAccumulator.add(
          static_cast<long double>(rowWeight) *
          static_cast<long double>(rowTime));
        minimumP = std::min(minimumP, rowP);
        maximumP = std::max(maximumP, rowP);
      }
    const long double weightSum = weightAccumulator.value();
    if (!(weightSum > 0.0L))
      throw std::runtime_error("SDDS selected-page weight sum is invalid");
    const long double referenceTime =
      weightedTimeAccumulator.value() / weightSum;

    hid_t file = H5Fcreate(options.output.c_str(), H5F_ACC_TRUNC,
                           H5P_DEFAULT, H5P_DEFAULT);
    requireHandle(file, "Cannot create output file: " + options.output);
    hid_t group = H5Gcreate2(file, "/particles", H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    requireHandle(group, "Cannot create /particles group");
    writeIntAttribute(group, "format_version", 4);
    writeIntAttribute(group, "sdds_protocol_version", headerVersion);
    writeIntAttribute(group, "sdds_page", options.page);
    writeDoubleAttribute(group, "source_arrival_time_reference_s",
                         static_cast<double>(referenceTime));
    writeStringAttribute(group, "coordinate_frame",
                         "fixed-lab-plane-crossing-events");
    writeStringAttribute(group, "arrival_time_reference",
                         "macro-weighted-mean-of-selected-sdds-page");
    writeStringAttribute(group, "source_columns", "x,xp,y,yp,t,p");
    writeStringAttribute(group, "source_data_mode",
      dataset->layout.data_mode.mode == SDDS_BINARY ? "binary" : "ascii");
    writeStringAttribute(group, "plane_position_unit", "m");
    writeStringAttribute(group, "arrival_time_offset_unit", "s");
    writeStringAttribute(group, "proper_velocity_unit", "gamma*v/c");
    writeStringAttribute(group, "macro_weight_definition",
      hasWeight ? "positive relative weight copied from SDDS column " +
                  options.weightColumn :
                  "unit relative weight for every selected SDDS row");
    writeStringAttribute(group, "source_id_definition",
      hasId ? "copied exactly from SDDS column " + options.idColumn :
              "one-based row number within selected SDDS page");

    const hsize_t dimensions[1] = {static_cast<hsize_t>(rows)};
    hid_t space = H5Screate_simple(1, dimensions, NULL);
    requireHandle(space, "Cannot create particle dataset space");
    hid_t fileType = fileRecordType();
    hid_t datasetOutput = H5Dcreate2(group, "records", fileType, space,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    H5Tclose(fileType);
    requireHandle(datasetOutput, "Cannot create /particles/records dataset");
    hid_t memoryType = memoryRecordType();

    const std::size_t chunkRecords = 65536;
    std::vector<ParticleRecord> buffer;
    buffer.reserve(std::min<std::uint64_t>(
      chunkRecords, static_cast<std::uint64_t>(rows)));
    std::uint64_t written = 0;
    for (std::int64_t row = 0; row < rows; ++row)
      {
        ParticleRecord record = {};
        record.planePosition[0] = x.value(row);
        record.planePosition[1] = y.value(row);
        const double slopeX = xp.value(row);
        const double slopeY = yp.value(row);
        const double p = momentum.value(row);
        const double rowTime = time.value(row);
        record.macroWeight = hasWeight ? weight->value(row) : 1.0;
        record.arrivalTimeOffset = static_cast<double>(
          static_cast<long double>(rowTime) - referenceTime);
        const double slopeNorm = std::sqrt(
          1.0 + slopeX * slopeX + slopeY * slopeY);
        record.properVelocity[2] = p / slopeNorm;
        record.properVelocity[0] = slopeX * record.properVelocity[2];
        record.properVelocity[1] = slopeY * record.properVelocity[2];
        record.sourceId = hasId ? id->unsignedInteger(row) :
          static_cast<std::uint64_t>(row + 1);
        buffer.push_back(record);
        if (buffer.size() == chunkRecords)
          {
            writeChunk(datasetOutput, memoryType, buffer, written);
            written += buffer.size();
            buffer.clear();
          }
      }
    writeChunk(datasetOutput, memoryType, buffer, written);
    written += buffer.size();

    delete id;
    delete weight;
    H5Tclose(memoryType);
    H5Dclose(datasetOutput);
    H5Gclose(group);
    requireStatus(H5Fclose(file), "Cannot close output HDF5 file");
    if (written != static_cast<std::uint64_t>(rows))
      throw std::runtime_error("SDDS-to-HDF5 record count mismatch");

    std::cout << std::setprecision(16)
      << "Converted SDDS protocol " << headerVersion << " "
      << (dataset->layout.data_mode.mode == SDDS_BINARY ? "binary" : "ASCII")
      << " page " << options.page << ": " << written
      << " particles -> " << options.output << "\n"
      << "Elegant mapping: (x,xp,y,yp,t,p) -> "
         "(x_plane,y_plane,t-t_ref,ux,uy,uz)\n"
      << "t_ref=" << static_cast<double>(referenceTime)
      << " s, p range=[" << minimumP << ", " << maximumP << "]\n"
      << "source IDs: "
      << (hasId ? options.idColumn : "selected-page row number")
      << ", macro weights: "
      << (hasWeight ? options.weightColumn : "equal") << std::endl;
  }

  int positiveInt(const std::string& text, const std::string& option)
  {
    char* end = NULL;
    const long value = std::strtol(text.c_str(), &end, 10);
    if (!end || *end != '\0' || value <= 0 || value > INT_MAX)
      throw std::invalid_argument(option + " requires a positive integer");
    return static_cast<int>(value);
  }

  void usage(const char* executable)
  {
    std::cout
      << "Usage: " << executable
      << " --input beam.sdds --output particles.h5 [options]\n"
      << "Options:\n"
      << "  --page N                 SDDS page to convert (default: 1)\n"
      << "  --id-column NAME|none    Stable ID column (default: particleID;\n"
      << "                           missing/default falls back to row number)\n"
      << "  --weight-column NAME|none  Positive relative macro weight\n"
      << "                           (default: equal weights)\n"
      << "Required Elegant columns: x[m], xp, y[m], yp, t[s], p[m$be$nc].\n";
  }
}

int main(int argc, char** argv)
{
  Options options;
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
      const std::string value(argv[++index]);
      if (argument == "--input") options.input = value;
      else if (argument == "--output") options.output = value;
      else if (argument == "--page")
        {
          try { options.page = positiveInt(value, argument); }
          catch (const std::exception& error)
            {
              std::cerr << error.what() << std::endl;
              return 1;
            }
        }
      else if (argument == "--id-column")
        options.idColumn = value == "none" ? "" : value;
      else if (argument == "--weight-column")
        options.weightColumn = value == "none" ? "" : value;
      else
        {
          std::cerr << "Unknown option: " << argument << std::endl;
          return 1;
        }
    }
  if (options.input.empty() || options.output.empty())
    {
      usage(argv[0]);
      return 1;
    }
  try
    {
      convert(options);
    }
  catch (const std::exception& error)
    {
      std::cerr << "Conversion failed: " << error.what() << std::endl;
      return 1;
    }
  return 0;
}
