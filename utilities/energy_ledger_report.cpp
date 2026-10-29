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

  int readIntAttribute(hid_t object, const char* name)
  {
    hid_t attribute = H5Aopen(object, name, H5P_DEFAULT);
    requireHandle(attribute, std::string("Cannot open ") + name);
    int value = 0;
    const herr_t status = H5Aread(attribute, H5T_NATIVE_INT, &value);
    H5Aclose(attribute);
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

  void printRecord(const char* label, const aprl::EnergyLedgerRecord& record,
                   bool hasAbsoluteBoundaryWork, bool hasMovingWork, bool hasBorisWork)
  {
    std::cout << label << ": step=" << record.step
      << " time_box_s=" << record.timeBox
      << " closure_valid=" << static_cast<unsigned int>(record.closureValid)
      << " active_macroparticles=" << record.activeMacroparticles
      << " removed_macroparticles=" << record.removedMacroparticles
      << " boundary_injection_entries="
      << record.boundaryInjectionEntries << "\n"
      << "  particle_kinetic_active_J=" << record.particleKineticActive
      << " particle_kinetic_removed_J=" << record.particleKineticRemoved
      << " field_energy_interior_J=" << record.fieldEnergyInterior << "\n"
      << "  outward_field_energy_J=" << record.outwardFieldEnergy
      << " prescribed_work_J=" << record.prescribedWork
      << " boundary_source_work_J=" << record.boundarySourceWork
      << " balance_residual_raw_J=" << record.balanceResidualRaw
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
    if (hasMovingWork)
      std::cout << "  longitudinal_boundary_outward_J=" << record.longitudinalBoundaryOutward
        << " moving_inserted_field_J=" << record.movingInsertedField
        << " moving_cropped_field_J=" << record.movingCroppedField
        << " moving_remap_field_change_J=" << record.movingRemapFieldChange
        << " moving_remap_interface_correction_J=" << record.movingRemapInterfaceCorrection
        << " inlet_overwrite_work_J=" << record.inletOverwriteWork
        << " frozen_constraint_work_J=" << record.frozenConstraintWork
        << " carrier_constraint_work_J=" << record.carrierConstraintWork << "\n";
    if (hasBorisWork)
      std::cout << "  delta_particle_kinetic_box_J=" << record.deltaParticleKineticBox
        << " delta_particle_kinetic_lab_J=" << record.deltaParticleKineticLab
        << " delta_particle_cPz_box_J=" << record.deltaParticleCPzBox
        << " particle_lorentz_delta_residual_J=" << record.particleLorentzDeltaResidual
        << " delta_field_energy_J=" << record.deltaFieldEnergy << "\n"
        << "  particle_grid_work_J=" << record.particleGridWork
        << " prescribed_boris_work_J=" << record.prescribedBorisWork
        << " magnetic_rotation_delta_K_J=" << record.magneticRotationDeltaK
        << " boris_roundoff_residual_J=" << record.borisRoundoffResidual
        << " boris_absolute_work_J=" << record.borisAbsoluteWork << "\n"
        << "  particle_work_balance_residual_J=" << record.particleWorkBalanceResidual
        << " field_exchange_residual_J=" << record.fieldExchangeResidual
        << " balance_residual_boris_J=" << record.balanceResidualBoris
        << " partial_work_substeps=" << record.partialWorkSubsteps << "\n"
        << "  released_kinetic_box_J=" << record.releasedKineticBox
        << " released_kinetic_lab_J=" << record.releasedKineticLab
        << " released_charge_C=" << record.releasedCharge << "\n";
    if (hasAbsoluteBoundaryWork)
      {
        std::cout << "  boundary_source_abs_work_J="
                  << record.boundarySourceAbsWork
                  << " boundary_source_abs_work_over_field_energy=";
        if (record.fieldEnergyInterior > 0.0)
          std::cout << record.boundarySourceAbsWork /
            record.fieldEnergyInterior;
        else
          std::cout << "undefined_zero_field_energy";
        std::cout << "\n";
      }
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
      const int formatVersion = readIntAttribute(group, "format_version");
      if (formatVersion < 1 || formatVersion > 5)
        throw std::runtime_error("Unsupported energy-ledger format version");
      const std::uint64_t committed = readUnsigned(
        group, "committed_records");
      const unsigned char complete = readByte(group, "complete");
      if (committed == 0)
        throw std::runtime_error("Energy ledger has no committed records");
      dataset = H5Dopen2(group, "records", H5P_DEFAULT);
      requireHandle(dataset, "Cannot open energy-ledger records");
      memoryType = aprl::EnergyLedgerWriter::createMemoryRecordType();
      aprl::EnergyLedgerRecord initial =
        readRecord(dataset, memoryType, 0);
      aprl::EnergyLedgerRecord final =
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

      if (formatVersion == 1)
        {
          initial.balanceResidualRaw = initial.balanceResidual;
          final.balanceResidualRaw = final.balanceResidual;
          if (hasSelected)
            selected.balanceResidualRaw = selected.balanceResidual;
        }

      std::cout << std::setprecision(12)
        << "format_version=" << formatVersion
        << " committed_records=" << committed
        << " complete=" << static_cast<unsigned int>(complete) << "\n";
      printRecord("initial", initial, formatVersion >= 3, formatVersion >= 4, formatVersion >= 5);
      if (hasSelected)
        {
          printRecord("selected", selected, formatVersion >= 3, formatVersion >= 4, formatVersion >= 5);
          std::cout
            << "initial_to_selected: particle_kinetic_accounted_J="
            << (selected.particleKineticActive +
                selected.particleKineticRemoved -
                initial.particleKineticActive)
            << " field_energy_J="
            << (selected.fieldEnergyInterior - initial.fieldEnergyInterior)
            << " outward_field_energy_J="
            << selected.outwardFieldEnergy
            << " boundary_source_work_J="
            << selected.boundarySourceWork
            << " residual_raw_J=" << selected.balanceResidualRaw
            << " residual_J=" << selected.balanceResidual
            << " mean_gamma_lab="
            << (selected.meanGammaLab - initial.meanGammaLab)
            << " sigma_gamma_lab="
            << (selected.sigmaGammaLab - initial.sigmaGammaLab) << "\n";
        }
      printRecord("final", final, formatVersion >= 3, formatVersion >= 4, formatVersion >= 5);
      if (formatVersion >= 3 && final.boundaryInjectionEntries > 0)
        std::cout << "note: boundary-source work / stored field energy is "
          "a source-activity scale, not the fraction of missing external "
          "near field or a bounded percentage.\n";
      std::cout
        << "change: particle_kinetic_accounted_J="
        << (final.particleKineticActive + final.particleKineticRemoved -
            initial.particleKineticActive)
        << " field_energy_J="
        << (final.fieldEnergyInterior - initial.fieldEnergyInterior)
        << " boundary_source_work_J=" << final.boundarySourceWork
        << " residual_raw_J=" << final.balanceResidualRaw
        << " residual_J=" << final.balanceResidual
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
