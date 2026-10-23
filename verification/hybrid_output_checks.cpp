#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>

namespace
{
  void require(bool condition, const std::string& message)
  {
    if (!condition) throw std::runtime_error(message);
  }

  class Handle
  {
  public:
    Handle(hid_t value, herr_t (*close)(hid_t)) : value_(value), close_(close)
    { require(value >= 0, "Cannot open HDF5 object"); }
    ~Handle() { close_(value_); }
    operator hid_t() const { return value_; }
  private:
    Handle(const Handle&) = delete;
    Handle& operator=(const Handle&) = delete;
    hid_t value_;
    herr_t (*close_)(hid_t);
  };

  herr_t collectDatasets(hid_t, const char* name, const H5O_info_t* info,
                         void* context)
  {
    if (info->type == H5O_TYPE_DATASET)
      static_cast<std::vector<std::string>*>(context)->push_back(name);
    return 0;
  }

  std::vector<std::string> datasets(hid_t file)
  {
    std::vector<std::string> result;
    require(H5Ovisit(file, H5_INDEX_NAME, H5_ITER_INC, collectDatasets,
                     &result
#if defined(H5Ovisit_vers) && H5Ovisit_vers > 1
                     // Newer HDF5 APIs require an explicit metadata selection.
                     , H5O_INFO_BASIC
#endif
                     ) >= 0, "Cannot enumerate HDF5 datasets");
    return result;
  }

  std::vector<hsize_t> dimensions(hid_t dataset)
  {
    Handle space(H5Dget_space(dataset), H5Sclose);
    const int rank = H5Sget_simple_extent_ndims(space);
    require(rank >= 0, "Cannot inspect HDF5 dimensions");
    std::vector<hsize_t> result(static_cast<std::size_t>(rank));
    if (rank)
      require(H5Sget_simple_extent_dims(space, result.data(), 0) >= 0,
        "Cannot read HDF5 dimensions");
    return result;
  }

  std::uint64_t scalar(hid_t file, const std::string& path)
  {
    Handle dataset(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose);
    require(dimensions(dataset).empty(), "Expected scalar: " + path);
    std::uint64_t value = 0;
    require(H5Dread(dataset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, &value) >= 0, "Cannot read scalar: " + path);
    return value;
  }

  struct FloatingMember
  {
    std::size_t offset;
    std::size_t width;
  };

  void floatingMembers(hid_t type, std::size_t offset,
                       std::vector<FloatingMember>& members)
  {
    const H5T_class_t kind = H5Tget_class(type);
    if (kind == H5T_FLOAT)
      {
        const std::size_t width = H5Tget_size(type);
        require(width == sizeof(float) || width == sizeof(double) ||
                  width == sizeof(long double), "Unsupported native float");
        members.push_back(FloatingMember{offset, width});
      }
    else if (kind == H5T_COMPOUND)
      for (int member = 0; member < H5Tget_nmembers(type); ++member)
        {
          Handle child(H5Tget_member_type(type, member), H5Tclose);
          floatingMembers(child, offset + H5Tget_member_offset(type, member),
                          members);
        }
    else if (kind == H5T_ARRAY)
      {
        const int rank = H5Tget_array_ndims(type);
        require(rank > 0, "Invalid HDF5 array type");
        std::vector<hsize_t> shape(static_cast<std::size_t>(rank));
        require(H5Tget_array_dims2(type, shape.data()) >= 0,
                "Cannot read array type dimensions");
        Handle child(H5Tget_super(type), H5Tclose);
        std::size_t count = 1;
        for (std::size_t axis = 0; axis < shape.size(); ++axis)
          count *= static_cast<std::size_t>(shape[axis]);
        for (std::size_t item = 0; item < count; ++item)
          floatingMembers(child, offset + item * H5Tget_size(child), members);
      }
  }

  template<class Number>
  bool finiteValue(const unsigned char* data)
  {
    Number value;
    std::memcpy(&value, data, sizeof(value));
    return std::isfinite(value);
  }

  void validateFinite(hid_t dataset, hid_t fileType, std::size_t count,
                      const std::string& name)
  {
    if (H5Tdetect_class(fileType, H5T_FLOAT) <= 0 || !count) return;
    Handle type(H5Tget_native_type(fileType, H5T_DIR_ASCEND), H5Tclose);
    std::vector<FloatingMember> members;
    floatingMembers(type, 0, members);
    const std::size_t stride = H5Tget_size(type);
    require(stride > 0 && count <=
              std::numeric_limits<std::size_t>::max() / stride,
            "Native dataset size overflows: " + name);
    std::vector<unsigned char> data(count * stride, 0);
    require(H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    data.data()) >= 0, "Cannot read native data: " + name);
    for (std::size_t item = 0; item < count; ++item)
      for (std::size_t member = 0; member < members.size(); ++member)
        {
          const FloatingMember& field = members[member];
          const unsigned char* value = data.data() + item * stride + field.offset;
          const bool finite = field.width == sizeof(float) ?
            finiteValue<float>(value) : field.width == sizeof(double) ?
            finiteValue<double>(value) : finiteValue<long double>(value);
          require(finite, "Nonfinite scientific value: " + name +
                    " record=" + std::to_string(item));
        }
  }

  void validateCommit(hid_t file, const std::string& group,
                      const std::string& countName,
                      const std::string& recordsName,
                      std::uint64_t expectedRecords)
  {
    require(scalar(file, group + "/complete") == 1,
      "Output is incomplete: " + group);
    const std::uint64_t count = scalar(file, group + "/" + countName);
    require(expectedRecords ? count == expectedRecords : count >= 2,
      "Unexpected committed record count: " + group + ": " +
      std::to_string(count));
    Handle records(H5Dopen2(file, (group + "/" + recordsName).c_str(),
                            H5P_DEFAULT), H5Dclose);
    const std::vector<hsize_t> shape = dimensions(records);
    require(!shape.empty() && shape[0] == count,
      "Committed extent does not match data: " + group);
  }

  void compareFile(const std::string& left, const std::string& right,
                   const std::string& group, const std::string& countName,
                   const std::string& recordsName,
                   std::uint64_t expectedRecords)
  {
    Handle first(H5Fopen(left.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    Handle second(H5Fopen(right.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    validateCommit(first, group, countName, recordsName, expectedRecords);
    validateCommit(second, group, countName, recordsName, expectedRecords);
    const std::vector<std::string> names = datasets(first);
    require(names == datasets(second), "Dataset names differ: " + left);
    require(!names.empty(), "Empty comparison file: " + left);
    for (std::size_t index = 0; index < names.size(); ++index)
      {
        const std::string& name = names[index];
        Handle a(H5Dopen2(first, name.c_str(), H5P_DEFAULT), H5Dclose);
        Handle b(H5Dopen2(second, name.c_str(), H5P_DEFAULT), H5Dclose);
        const std::vector<hsize_t> shape = dimensions(a);
        require(shape == dimensions(b), "Dataset shape differs: " + name);
        Handle typeA(H5Dget_type(a), H5Tclose);
        Handle typeB(H5Dget_type(b), H5Tclose);
        require(H5Tequal(typeA, typeB) > 0, "Dataset type differs: " + name);
        require(H5Tdetect_class(typeA, H5T_VLEN) == 0 &&
                H5Tis_variable_str(typeA) == 0,
          "Variable-size dataset is unsupported: " + name);
        std::size_t bytes = H5Tget_size(typeA);
        require(bytes > 0, "Zero-sized HDF5 type: " + name);
        for (std::size_t axis = 0; axis < shape.size(); ++axis)
          {
            require(shape[axis] <= std::numeric_limits<std::size_t>::max() /
                      std::max<std::size_t>(bytes, 1),
              "Dataset byte count overflows: " + name);
            bytes *= static_cast<std::size_t>(shape[axis]);
          }
        // Read using the packed on-disk type, with zeroed padding. Comparison
        // includes every member of compound particle, field and ledger data.
        std::vector<unsigned char> dataA(bytes, 0), dataB(bytes, 0);
        if (bytes)
          {
            require(H5Dread(a, typeA, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            dataA.data()) >= 0, "Cannot read: " + name);
            require(H5Dread(b, typeB, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            dataB.data()) >= 0, "Cannot read: " + name);
          }
        require(dataA == dataB,
          "Hybrid output differs bitwise: " + left + ": " + name);
        // Equal byte arrays need only one finite-value scan. Inspect nested
        // compound and array members in the native HDF5 representation.
        validateFinite(a, typeA, bytes / H5Tget_size(typeA), name);
      }
    std::cout << "bitwise_equal=" << left << " datasets=" << names.size()
              << std::endl;
  }
}

int main(int argc, char** argv)
{
  try
    {
      require(argc == 3, "Usage: aprl_hybrid_output_checks SERIAL_DIR HYBRID_DIR");
      const std::string first(argv[1]), second(argv[2]);
      const char* particlePlanes[2] = {"activation", "downstream"};
      for (unsigned int plane = 0; plane < 2; ++plane)
        {
          const std::string suffix = std::string("/detectors/") +
            particlePlanes[plane] + ".h5";
          compareFile(first + suffix, second + suffix, "/particle_plane",
            "committed_records", "records", 10410);
        }
      compareFile(first + "/detectors/field.h5", second + "/detectors/field.h5",
        "/field_plane", "committed_samples", "fields", 0);
      compareFile(first + "/ledger/energy-ledger.h5",
        second + "/ledger/energy-ledger.h5", "/energy_ledger",
        "committed_records", "records", 0);
      return 0;
    }
  catch (const std::exception& error)
    {
      std::cerr << "Hybrid verification failure: " << error.what() << std::endl;
      return 1;
    }
}
