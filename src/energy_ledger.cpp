#include "energy_ledger.h"

#include <cmath>
#include <stdexcept>

namespace fel
{
  namespace
  {
    void requireHandle(hid_t handle, const std::string& message)
    {
      if (handle < 0) throw std::runtime_error(message);
    }

    void requireStatus(herr_t status, const std::string& message)
    {
      if (status < 0) throw std::runtime_error(message);
    }

    struct DoubleField
    {
      const char* name;
      std::size_t memoryOffset;
    };

    const DoubleField DOUBLE_FIELDS[] = {
      {"time_box_s", HOFFSET(EnergyLedgerRecord, timeBox)},
      {"active_represented_electrons", HOFFSET(EnergyLedgerRecord, activeRepresentedElectrons)},
      {"particle_kinetic_active_J", HOFFSET(EnergyLedgerRecord, particleKineticActive)},
      {"particle_total_active_J", HOFFSET(EnergyLedgerRecord, particleTotalActive)},
      {"particle_kinetic_removed_J", HOFFSET(EnergyLedgerRecord, particleKineticRemoved)},
      {"particle_total_removed_J", HOFFSET(EnergyLedgerRecord, particleTotalRemoved)},
      {"field_energy_interior_J", HOFFSET(EnergyLedgerRecord, fieldEnergyInterior)},
      {"outward_field_energy_J", HOFFSET(EnergyLedgerRecord, outwardFieldEnergy)},
      {"prescribed_work_J", HOFFSET(EnergyLedgerRecord, prescribedWork)},
      {"balance_residual_J", HOFFSET(EnergyLedgerRecord, balanceResidual)},
      {"relative_balance_exchange", HOFFSET(EnergyLedgerRecord, relativeBalanceExchange)},
      {"relative_balance_initial", HOFFSET(EnergyLedgerRecord, relativeBalanceInitial)},
      {"field_fraction_kinetic_plus_field", HOFFSET(EnergyLedgerRecord, fieldFractionKinetic)},
      {"field_fraction_including_rest", HOFFSET(EnergyLedgerRecord, fieldFractionIncludingRest)},
      {"mean_gamma_lab", HOFFSET(EnergyLedgerRecord, meanGammaLab)},
      {"sigma_gamma_lab", HOFFSET(EnergyLedgerRecord, sigmaGammaLab)},
      {"relative_energy_spread_lab", HOFFSET(EnergyLedgerRecord, relativeEnergySpreadLab)},
      {"minimum_gamma_lab", HOFFSET(EnergyLedgerRecord, minimumGammaLab)},
      {"maximum_gamma_lab", HOFFSET(EnergyLedgerRecord, maximumGammaLab)},
      {"linear_chirp_gamma_per_m", HOFFSET(EnergyLedgerRecord, linearChirpGammaPerM)},
      {"uncorrelated_sigma_gamma_lab", HOFFSET(EnergyLedgerRecord, uncorrelatedSigmaGammaLab)},
      {"mean_gamma_box", HOFFSET(EnergyLedgerRecord, meanGammaBox)},
      {"sigma_gamma_box", HOFFSET(EnergyLedgerRecord, sigmaGammaBox)}
    };

    const std::size_t DOUBLE_FIELD_COUNT =
      sizeof(DOUBLE_FIELDS) / sizeof(DOUBLE_FIELDS[0]);
  }

  EnergyLedgerRecord::EnergyLedgerRecord()
    : step(0), activeMacroparticles(0), removedMacroparticles(0),
      closureValid(1), timeBox(0.0), activeRepresentedElectrons(0.0),
      particleKineticActive(0.0), particleTotalActive(0.0),
      particleKineticRemoved(0.0), particleTotalRemoved(0.0),
      fieldEnergyInterior(0.0), outwardFieldEnergy(0.0),
      prescribedWork(0.0), balanceResidual(0.0),
      relativeBalanceExchange(0.0), relativeBalanceInitial(0.0),
      fieldFractionKinetic(0.0), fieldFractionIncludingRest(0.0),
      meanGammaLab(0.0), sigmaGammaLab(0.0),
      relativeEnergySpreadLab(0.0), minimumGammaLab(0.0),
      maximumGammaLab(0.0), linearChirpGammaPerM(0.0),
      uncorrelatedSigmaGammaLab(0.0), meanGammaBox(0.0),
      sigmaGammaBox(0.0)
  {
    for (unsigned int face = 0; face < 6; ++face)
      outwardFieldEnergyFace[face] = 0.0;
  }

  EnergyLedgerWriter::EnergyLedgerWriter()
    : isOpen_(false), durableProgress_(false), bufferLimit_(0),
      committedRecords_(0), buffer_(), file_(-1), group_(-1),
      recordDataset_(-1), committedDataset_(-1), completeDataset_(-1),
      memoryRecordType_(-1), fileRecordType_(-1)
  {}

  EnergyLedgerWriter::~EnergyLedgerWriter()
  {
    try { close(false); }
    catch (...) { closeHandles(); }
  }

  hid_t EnergyLedgerWriter::createMemoryRecordType()
  {
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(EnergyLedgerRecord));
    requireHandle(type, "Cannot create energy-ledger memory datatype");
    try
      {
        requireStatus(H5Tinsert(type, "step",
          HOFFSET(EnergyLedgerRecord, step), H5T_NATIVE_UINT64),
          "Cannot add energy-ledger step");
        requireStatus(H5Tinsert(type, "active_macroparticles",
          HOFFSET(EnergyLedgerRecord, activeMacroparticles),
          H5T_NATIVE_UINT64), "Cannot add active particle count");
        requireStatus(H5Tinsert(type, "removed_macroparticles",
          HOFFSET(EnergyLedgerRecord, removedMacroparticles),
          H5T_NATIVE_UINT64), "Cannot add removed particle count");
        requireStatus(H5Tinsert(type, "closure_valid",
          HOFFSET(EnergyLedgerRecord, closureValid), H5T_NATIVE_UCHAR),
          "Cannot add energy-ledger validity flag");
        for (std::size_t field = 0; field < DOUBLE_FIELD_COUNT; ++field)
          requireStatus(H5Tinsert(type, DOUBLE_FIELDS[field].name,
            DOUBLE_FIELDS[field].memoryOffset, H5T_NATIVE_DOUBLE),
            std::string("Cannot add energy-ledger field ") +
              DOUBLE_FIELDS[field].name);
        hsize_t faceCount[1] = {6};
        hid_t faceType = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, faceCount);
        requireHandle(faceType, "Cannot create energy-ledger face type");
        const herr_t status = H5Tinsert(type,
          "outward_field_energy_face_J",
          HOFFSET(EnergyLedgerRecord, outwardFieldEnergyFace), faceType);
        H5Tclose(faceType);
        requireStatus(status, "Cannot add energy-ledger face energies");
      }
    catch (...)
      {
        H5Tclose(type);
        throw;
      }
    return type;
  }

  hid_t EnergyLedgerWriter::createFileRecordType() const
  {
    const std::size_t u64 = 8;
    const std::size_t f64 = 8;
    const std::size_t fileSize = 3 * u64 + 1 +
      (DOUBLE_FIELD_COUNT + 6) * f64;
    hid_t type = H5Tcreate(H5T_COMPOUND, fileSize);
    requireHandle(type, "Cannot create energy-ledger file datatype");
    try
      {
        std::size_t offset = 0;
        requireStatus(H5Tinsert(type, "step", offset, H5T_STD_U64LE),
          "Cannot add file step");
        offset += u64;
        requireStatus(H5Tinsert(type, "active_macroparticles", offset,
          H5T_STD_U64LE), "Cannot add file active particle count");
        offset += u64;
        requireStatus(H5Tinsert(type, "removed_macroparticles", offset,
          H5T_STD_U64LE), "Cannot add file removed particle count");
        offset += u64;
        requireStatus(H5Tinsert(type, "closure_valid", offset,
          H5T_STD_U8LE), "Cannot add file validity flag");
        offset += 1;
        for (std::size_t field = 0; field < DOUBLE_FIELD_COUNT; ++field)
          {
            requireStatus(H5Tinsert(type, DOUBLE_FIELDS[field].name,
              offset, H5T_IEEE_F64LE),
              std::string("Cannot add energy-ledger file field ") +
                DOUBLE_FIELDS[field].name);
            offset += f64;
          }
        hsize_t faceCount[1] = {6};
        hid_t faceType = H5Tarray_create2(H5T_IEEE_F64LE, 1, faceCount);
        requireHandle(faceType, "Cannot create file face type");
        const herr_t status = H5Tinsert(type,
          "outward_field_energy_face_J", offset, faceType);
        H5Tclose(faceType);
        requireStatus(status, "Cannot add file face energies");
      }
    catch (...)
      {
        H5Tclose(type);
        throw;
      }
    return type;
  }

  void EnergyLedgerWriter::open(const std::string& filename, int mpiSize,
      std::size_t bufferRecords, unsigned int compressionLevel,
      bool durableProgress, bool overwrite, const RunMetadata* metadata)
  {
    if (isOpen_) throw std::runtime_error("Energy ledger is already open");
    if (filename.empty())
      throw std::invalid_argument("Energy-ledger filename cannot be empty");
    if (mpiSize < 1 || bufferRecords == 0 || compressionLevel > 9)
      throw std::invalid_argument("Invalid energy-ledger output settings");

    durableProgress_ = durableProgress;
    bufferLimit_ = bufferRecords;
    committedRecords_ = 0;
    buffer_.clear();
    buffer_.reserve(bufferLimit_);

    hid_t access = H5Pcreate(H5P_FILE_ACCESS);
    requireHandle(access, "Cannot create energy-ledger file access");
#if H5_VERSION_GE(1, 10, 0)
    requireStatus(H5Pset_libver_bounds(access, H5F_LIBVER_LATEST,
      H5F_LIBVER_LATEST), "Cannot select crash-readable HDF5 format");
#endif
    file_ = H5Fcreate(filename.c_str(),
      overwrite ? H5F_ACC_TRUNC : H5F_ACC_EXCL,
      H5P_DEFAULT, access);
    H5Pclose(access);
    requireHandle(file_, "Cannot create energy-ledger file: " + filename);

    try
      {
        group_ = H5Gcreate2(file_, "/energy_ledger", H5P_DEFAULT,
          H5P_DEFAULT, H5P_DEFAULT);
        requireHandle(group_, "Cannot create energy-ledger group");
        memoryRecordType_ = createMemoryRecordType();
        fileRecordType_ = createFileRecordType();
        hsize_t initial[1] = {0};
        hsize_t maximum[1] = {H5S_UNLIMITED};
        hid_t space = H5Screate_simple(1, initial, maximum);
        requireHandle(space, "Cannot create energy-ledger dataspace");
        hid_t creation = H5Pcreate(H5P_DATASET_CREATE);
        requireHandle(creation, "Cannot create energy-ledger chunk settings");
        hsize_t chunk[1] = {static_cast<hsize_t>(bufferLimit_)};
        requireStatus(H5Pset_chunk(creation, 1, chunk),
          "Cannot set energy-ledger chunk size");
        requireStatus(H5Pset_shuffle(creation),
          "Cannot enable energy-ledger shuffle");
        if (compressionLevel > 0)
          requireStatus(H5Pset_deflate(creation, compressionLevel),
            "Cannot enable energy-ledger compression");
        recordDataset_ = H5Dcreate2(group_, "records", fileRecordType_,
          space, H5P_DEFAULT, creation, H5P_DEFAULT);
        H5Pclose(creation);
        H5Sclose(space);
        requireHandle(recordDataset_, "Cannot create energy-ledger records");

        hid_t scalar = H5Screate(H5S_SCALAR);
        requireHandle(scalar, "Cannot create energy-ledger scalar space");
        committedDataset_ = H5Dcreate2(group_, "committed_records",
          H5T_STD_U64LE, scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        completeDataset_ = H5Dcreate2(group_, "complete", H5T_STD_U8LE,
          scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(scalar);
        requireHandle(committedDataset_, "Cannot create ledger commit marker");
        requireHandle(completeDataset_, "Cannot create ledger completion marker");
        writeUnsignedScalar(committedDataset_, 0);
        writeByteScalar(completeDataset_, 0);

        writeIntAttribute(group_, "format_version", 1);
        writeIntAttribute(group_, "mpi_size", mpiSize);
        writeStringAttribute(group_, "frame", "boosted computational frame");
        writeStringAttribute(group_, "face_order", "x-,x+,y-,y+,z-,z+");
        writeStringAttribute(group_, "balance_equation",
          "delta(active_kinetic+removed_kinetic+interior_field)+outward_field-prescribed_work=residual");
        writeStringAttribute(group_, "field_surface",
          "first interior cell-centre layer next to CPML; grid outer layer for PEC");
        writeStringAttribute(group_, "staggering_note",
          "cell-centred E and half-step B give a convergence diagnostic, not an exact algebraic Yee invariant");
        writeStringAttribute(group_, "spread_note",
          "lab momentum statistics on an equal-box-time slice; particle planes are required for equal-location accelerator diagnostics");
        writeStringAttribute(group_, "reader_contract",
          "read only records[0:committed_records]");
        if (metadata)
          {
            writeStringAttribute(group_, "run_id", metadata->runId);
            writeStringAttribute(group_, "configuration_digest_fnv1a64",
              metadata->configurationDigest);
            writeStringAttribute(group_, "source_revision",
              metadata->sourceRevision);
            writeStringAttribute(group_, "run_manifest",
              metadata->manifestPath);
          }
        if (durableProgress_)
          requireStatus(H5Fflush(file_, H5F_SCOPE_GLOBAL),
            "Cannot initialize energy-ledger file on disk");
        isOpen_ = true;
      }
    catch (...)
      {
        closeHandles();
        buffer_.clear();
        throw;
      }
  }

  void EnergyLedgerWriter::append(const EnergyLedgerRecord& record)
  {
    if (!isOpen_)
      throw std::runtime_error("Energy-ledger append requires an open file");
    const double scalar[] = {
      record.timeBox, record.activeRepresentedElectrons,
      record.particleKineticActive, record.particleTotalActive,
      record.particleKineticRemoved, record.particleTotalRemoved,
      record.fieldEnergyInterior, record.outwardFieldEnergy,
      record.prescribedWork, record.balanceResidual,
      record.relativeBalanceExchange, record.relativeBalanceInitial,
      record.fieldFractionKinetic, record.fieldFractionIncludingRest,
      record.meanGammaLab, record.sigmaGammaLab,
      record.relativeEnergySpreadLab, record.minimumGammaLab,
      record.maximumGammaLab, record.linearChirpGammaPerM,
      record.uncorrelatedSigmaGammaLab, record.meanGammaBox,
      record.sigmaGammaBox
    };
    for (std::size_t field = 0;
         field < sizeof(scalar) / sizeof(scalar[0]); ++field)
      if (!std::isfinite(scalar[field]))
        throw std::invalid_argument("Energy-ledger record is not finite");
    for (unsigned int face = 0; face < 6; ++face)
      if (!std::isfinite(record.outwardFieldEnergyFace[face]))
        throw std::invalid_argument("Energy-ledger face energy is not finite");
    buffer_.push_back(record);
    if (buffer_.size() >= bufferLimit_) flush();
  }

  void EnergyLedgerWriter::flush()
  {
    if (!isOpen_)
      throw std::runtime_error("Energy-ledger flush requires an open file");
    if (buffer_.empty()) return;
    const hsize_t start[1] = {static_cast<hsize_t>(committedRecords_)};
    const hsize_t count[1] = {static_cast<hsize_t>(buffer_.size())};
    const hsize_t extent[1] = {start[0] + count[0]};
    requireStatus(H5Dset_extent(recordDataset_, extent),
      "Cannot extend energy-ledger dataset");
    hid_t fileSpace = H5Dget_space(recordDataset_);
    requireHandle(fileSpace, "Cannot select energy-ledger file space");
    try
      {
        requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
          start, NULL, count, NULL), "Cannot select energy-ledger batch");
        hid_t memorySpace = H5Screate_simple(1, count, NULL);
        requireHandle(memorySpace, "Cannot create energy-ledger memory space");
        const herr_t status = H5Dwrite(recordDataset_, memoryRecordType_,
          memorySpace, fileSpace, H5P_DEFAULT, &buffer_[0]);
        H5Sclose(memorySpace);
        requireStatus(status, "Cannot write energy-ledger batch");
      }
    catch (...)
      {
        H5Sclose(fileSpace);
        throw;
      }
    H5Sclose(fileSpace);
    if (durableProgress_)
      requireStatus(H5Fflush(file_, H5F_SCOPE_GLOBAL),
        "Cannot make energy-ledger records durable");
    committedRecords_ += static_cast<std::uint64_t>(buffer_.size());
    writeUnsignedScalar(committedDataset_, committedRecords_);
    if (durableProgress_)
      requireStatus(H5Fflush(file_, H5F_SCOPE_GLOBAL),
        "Cannot make energy-ledger marker durable");
    buffer_.clear();
  }

  void EnergyLedgerWriter::close(bool completed)
  {
    if (!isOpen_)
      {
        closeHandles();
        return;
      }
    flush();
    writeByteScalar(completeDataset_, completed ? 1 : 0);
    requireStatus(H5Fflush(file_, H5F_SCOPE_GLOBAL),
      "Cannot finalize energy-ledger file");
    isOpen_ = false;
    closeHandles();
    buffer_.clear();
  }

  void EnergyLedgerWriter::writeUnsignedScalar(hid_t dataset,
      std::uint64_t value) const
  {
    requireStatus(H5Dwrite(dataset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &value), "Cannot write energy-ledger counter");
  }

  void EnergyLedgerWriter::writeByteScalar(hid_t dataset,
      unsigned char value) const
  {
    requireStatus(H5Dwrite(dataset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &value), "Cannot write energy-ledger marker");
  }

  void EnergyLedgerWriter::writeIntAttribute(hid_t object,
      const char* name, int value) const
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create energy-ledger attribute space");
    hid_t attribute = H5Acreate2(object, name, H5T_STD_I32LE, space,
      H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, "Cannot create energy-ledger integer attribute");
    const herr_t status = H5Awrite(attribute, H5T_NATIVE_INT, &value);
    H5Aclose(attribute);
    requireStatus(status, "Cannot write energy-ledger integer attribute");
  }

  void EnergyLedgerWriter::writeStringAttribute(hid_t object,
      const char* name, const std::string& value) const
  {
    hid_t type = H5Tcopy(H5T_C_S1);
    requireHandle(type, "Cannot create energy-ledger string type");
    requireStatus(H5Tset_size(type, value.size() + 1),
      "Cannot size energy-ledger string type");
    requireStatus(H5Tset_strpad(type, H5T_STR_NULLTERM),
      "Cannot set energy-ledger string padding");
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create energy-ledger string space");
    hid_t attribute = H5Acreate2(object, name, type, space,
      H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, "Cannot create energy-ledger string attribute");
    const herr_t status = H5Awrite(attribute, type, value.c_str());
    H5Aclose(attribute);
    H5Tclose(type);
    requireStatus(status, "Cannot write energy-ledger string attribute");
  }

  void EnergyLedgerWriter::closeHandles() noexcept
  {
    if (recordDataset_ >= 0) H5Dclose(recordDataset_);
    if (committedDataset_ >= 0) H5Dclose(committedDataset_);
    if (completeDataset_ >= 0) H5Dclose(completeDataset_);
    if (memoryRecordType_ >= 0) H5Tclose(memoryRecordType_);
    if (fileRecordType_ >= 0) H5Tclose(fileRecordType_);
    if (group_ >= 0) H5Gclose(group_);
    if (file_ >= 0) H5Fclose(file_);
    recordDataset_ = committedDataset_ = completeDataset_ = -1;
    memoryRecordType_ = fileRecordType_ = -1;
    group_ = file_ = -1;
    isOpen_ = false;
  }

  bool EnergyLedgerWriter::isOpen() const { return isOpen_; }

  std::size_t EnergyLedgerWriter::memoryBytes() const
  {
    return buffer_.capacity() * sizeof(EnergyLedgerRecord);
  }
}
