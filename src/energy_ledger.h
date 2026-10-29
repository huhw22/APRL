#ifndef DIRECT_EB_ENERGY_LEDGER_H
#define DIRECT_EB_ENERGY_LEDGER_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "hdf5.h"
#include "run_metadata.h"

namespace aprl
{
  /* One global, boosted-frame energy audit sample.  The accelerator
   * statistics are laboratory momenta evaluated on the same box-time slice;
   * a particle plane remains the rigorous equal-lab-location diagnostic. */
  struct EnergyLedgerRecord
  {
    std::uint64_t step;
    std::uint64_t activeMacroparticles;
    std::uint64_t removedMacroparticles;
    std::uint64_t boundaryInjectionEntries;
    std::uint8_t closureValid;

    double timeBox;
    double activeRepresentedElectrons;
    double particleKineticActive;
    double particleTotalActive;
    double particleKineticRemoved;
    double particleTotalRemoved;
    double fieldEnergyInterior;
    double outwardFieldEnergyFace[6];
    double outwardFieldEnergy;
    double prescribedWork;
    double boundarySourceWork;
    double boundarySourceAbsWork;
    // Moving-volume terms. Positive remap/overwrite adds energy to the grid;
    // longitudinalBoundaryOutward combines propagation, remap and overwrite.
    double movingInsertedField;
    double movingCroppedField;
    double movingRemapFieldChange;
    double movingRemapInterfaceCorrection;
    double inletOverwriteWork;
    double longitudinalBoundaryOutward;
    double movingEnergyDiagnosticRankSeconds;
    double frozenConstraintWork;
    double carrierConstraintWork;
    double constraintEnergyDiagnosticRankSeconds;
    double balanceResidualRaw;
    double balanceResidual;
    double relativeBalanceExchange;
    double relativeBalanceInitial;
    double fieldFractionKinetic;
    double fieldFractionIncludingRest;

    double meanGammaLab;
    double sigmaGammaLab;
    double relativeEnergySpreadLab;
    double minimumGammaLab;
    double maximumGammaLab;
    double linearChirpGammaPerM;
    double uncorrelatedSigmaGammaLab;
    double meanGammaBox;
    double sigmaGammaBox;

    EnergyLedgerRecord();
  };

  /* Root-only, append-only HDF5 writer.  The commit marker makes interrupted
   * interactive output readable without forcing every sample to disk. */
  class EnergyLedgerWriter
  {
  public:
    EnergyLedgerWriter();
    ~EnergyLedgerWriter();

    EnergyLedgerWriter(const EnergyLedgerWriter&) = delete;
    EnergyLedgerWriter& operator=(const EnergyLedgerWriter&) = delete;

    void open(const std::string& filename, int mpiSize,
              std::size_t bufferRecords, unsigned int compressionLevel,
              bool durableProgress, bool overwrite = false,
              const RunMetadata* metadata = NULL,
              bool approximateBoundaryInjection = false,
              bool movingVolume = false);
    void append(const EnergyLedgerRecord& record);
    void flush();
    /* Explicit snapshot durability; ordinary throughput flush is unchanged. */
    void flushForCheckpoint();
    void close(bool completed = false);

    bool isOpen() const;
    std::uint64_t committedRecords() const;
    std::size_t memoryBytes() const;

    /* Exposed for the small C++ report utility so readers use exactly the
     * same native compound mapping as the writer.  The caller closes it. */
    static hid_t createMemoryRecordType();

  private:
    hid_t createFileRecordType() const;
    void writeUnsignedScalar(hid_t dataset, std::uint64_t value) const;
    void writeByteScalar(hid_t dataset, unsigned char value) const;
    void writeIntAttribute(hid_t object, const char* name, int value) const;
    void writeStringAttribute(hid_t object, const char* name,
                              const std::string& value) const;
    void closeHandles() noexcept;

    bool isOpen_;
    bool durableProgress_;
    std::size_t bufferLimit_;
    std::uint64_t committedRecords_;
    std::vector<EnergyLedgerRecord> buffer_;
    hid_t file_;
    hid_t group_;
    hid_t recordDataset_;
    hid_t committedDataset_;
    hid_t completeDataset_;
    hid_t memoryRecordType_;
    hid_t fileRecordType_;
  };
}

#endif
