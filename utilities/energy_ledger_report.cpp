#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "energy_ledger.h"

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

  std::uint64_t readUnsigned(hid_t group, const char* name)
  {
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
    requireHandle(dataset, std::string("Cannot open ") + name);
    std::uint64_t value = 0;
    const herr_t status = H5Dread(dataset, H5T_NATIVE_UINT64,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
    H5Dclose(dataset);
    requireStatus(status, std::string("Cannot read ") + name);
    return value;
  }

  unsigned char readByte(hid_t group, const char* name)
  {
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
    requireHandle(dataset, std::string("Cannot open ") + name);
    unsigned char value = 0;
    const herr_t status = H5Dread(dataset, H5T_NATIVE_UCHAR,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
    H5Dclose(dataset);
    requireStatus(status, std::string("Cannot read ") + name);
    return value;
  }

  aprl::EnergyLedgerRecord readRecord(hid_t dataset, hid_t memoryType,
                                     std::uint64_t index)
  {
    hid_t fileSpace = H5Dget_space(dataset);
    requireHandle(fileSpace, "Cannot open energy-ledger dataspace");
    const hsize_t start[1] = {static_cast<hsize_t>(index)};
    const hsize_t count[1] = {1};
    requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
      start, NULL, count, NULL), "Cannot select energy-ledger record");
    hid_t memorySpace = H5Screate_simple(1, count, NULL);
    requireHandle(memorySpace, "Cannot create record memory space");
    aprl::EnergyLedgerRecord record;
    const herr_t status = H5Dread(dataset, memoryType, memorySpace,
      fileSpace, H5P_DEFAULT, &record);
    H5Sclose(memorySpace);
    H5Sclose(fileSpace);
    requireStatus(status, "Cannot read energy-ledger record");
    return record;
  }

  void printRecord(const char* label, const aprl::EnergyLedgerRecord& record)
  {
    std::cout << label << ": step=" << record.step
      << " time_box_s=" << record.timeBox
      << " closure_valid=" << static_cast<unsigned int>(record.closureValid)
      << " active_macroparticles=" << record.activeMacroparticles
      << " removed_macroparticles=" << record.removedMacroparticles << "\n"
      << "  particle_kinetic_active_J=" << record.particleKineticActive
      << " particle_kinetic_removed_J=" << record.particleKineticRemoved
      << " field_energy_interior_J=" << record.fieldEnergyInterior << "\n"
      << "  outward_field_energy_J=" << record.outwardFieldEnergy
      << " prescribed_work_J=" << record.prescribedWork
      << " balance_residual_J=" << record.balanceResidual << "\n"
      << "  relative_balance_exchange=" << record.relativeBalanceExchange
      << " relative_balance_initial=" << record.relativeBalanceInitial
      << " field_fraction_kinetic_plus_field="
      << record.fieldFractionKinetic
      << " field_fraction_including_rest="
      << record.fieldFractionIncludingRest << "\n"
      << "  mean_gamma_lab=" << record.meanGammaLab
      << " sigma_gamma_lab=" << record.sigmaGammaLab
      << " relative_energy_spread_lab="
      << record.relativeEnergySpreadLab
      << " uncorrelated_sigma_gamma_lab="
      << record.uncorrelatedSigmaGammaLab
      << " linear_chirp_gamma_per_m="
      << record.linearChirpGammaPerM << "\n";
  }
}

int main(int argc, char** argv)
{
  if (argc != 2 && argc != 3)
    {
      std::cerr << "Usage: energy_ledger_report LEDGER.h5 [RECORD_INDEX]\n";
      return 2;
    }
  hid_t file = -1;
  hid_t group = -1;
  hid_t dataset = -1;
  hid_t memoryType = -1;
  try
    {
      file = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
      requireHandle(file, "Cannot open energy-ledger file");
      group = H5Gopen2(file, "/energy_ledger", H5P_DEFAULT);
      requireHandle(group, "Cannot open /energy_ledger");
      const std::uint64_t committed = readUnsigned(
        group, "committed_records");
      const unsigned char complete = readByte(group, "complete");
      if (committed == 0)
        throw std::runtime_error("Energy ledger has no committed records");
      dataset = H5Dopen2(group, "records", H5P_DEFAULT);
      requireHandle(dataset, "Cannot open energy-ledger records");
      memoryType = aprl::EnergyLedgerWriter::createMemoryRecordType();
      const aprl::EnergyLedgerRecord initial =
        readRecord(dataset, memoryType, 0);
      const aprl::EnergyLedgerRecord final =
        readRecord(dataset, memoryType, committed - 1);
      bool hasSelected = argc == 3;
      std::uint64_t selectedIndex = 0;
      aprl::EnergyLedgerRecord selected;
      if (hasSelected)
        {
          char* end = 0;
          const unsigned long long parsed = std::strtoull(argv[2], &end, 10);
          if (!argv[2][0] || !end || *end != '\0' || parsed >= committed)
            throw std::runtime_error(
              "RECORD_INDEX must be a committed zero-based record index");
          selectedIndex = static_cast<std::uint64_t>(parsed);
          selected = readRecord(dataset, memoryType, selectedIndex);
        }

      std::cout << std::setprecision(12)
        << "committed_records=" << committed
        << " complete=" << static_cast<unsigned int>(complete) << "\n";
      printRecord("initial", initial);
      if (hasSelected)
        {
          printRecord("selected", selected);
          std::cout
            << "initial_to_selected: particle_kinetic_accounted_J="
            << (selected.particleKineticActive +
                selected.particleKineticRemoved -
                initial.particleKineticActive)
            << " field_energy_J="
            << (selected.fieldEnergyInterior - initial.fieldEnergyInterior)
            << " outward_field_energy_J="
            << selected.outwardFieldEnergy
            << " residual_J=" << selected.balanceResidual
            << " mean_gamma_lab="
            << (selected.meanGammaLab - initial.meanGammaLab)
            << " sigma_gamma_lab="
            << (selected.sigmaGammaLab - initial.sigmaGammaLab) << "\n";
        }
      printRecord("final", final);
      std::cout
        << "change: particle_kinetic_accounted_J="
        << (final.particleKineticActive + final.particleKineticRemoved -
            initial.particleKineticActive)
        << " field_energy_J="
        << (final.fieldEnergyInterior - initial.fieldEnergyInterior)
        << " mean_gamma_lab="
        << (final.meanGammaLab - initial.meanGammaLab)
        << " sigma_gamma_lab="
        << (final.sigmaGammaLab - initial.sigmaGammaLab) << "\n";

      H5Tclose(memoryType);
      H5Dclose(dataset);
      H5Gclose(group);
      H5Fclose(file);
      return 0;
    }
  catch (const std::exception& error)
    {
      if (memoryType >= 0) H5Tclose(memoryType);
      if (dataset >= 0) H5Dclose(dataset);
      if (group >= 0) H5Gclose(group);
      if (file >= 0) H5Fclose(file);
      std::cerr << "energy_ledger_report: " << error.what() << "\n";
      return 1;
    }
}
