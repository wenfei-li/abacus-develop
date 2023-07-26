#ifndef LCAO_DFTU_NEW_INTERFACE_H
#define LCAO_DFTU_NEW_INTERFACE_H

#include "LCAO_dftu_new.h"
#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"
#include <memory>

class LCAO_DftU_New_Interface
{
  public:
    /// @brief Constructor for LCAO_DftU_New_Interface
    /// @param ld_in
    LCAO_DftU_New_Interface(std::shared_ptr<LCAO_DftU_New> ld_in);
    /// @brief output dftu-related labels, and energy corrections
    /// @param[in] etot
    /// @param[in] nks
    /// @param[in] nat
    /// @param[in] ekb
    /// @param[in] kvec_d
    /// @param[in] ucell
    /// @param[in] orb
    /// @param[in] GridD
    /// @param[in] ParaV
    /// @param[in] psi
    /// @param[in] psid
    /// @param[in] dm_gamma
    /// @param[in] dm_k
    void out_deepks_labels(double etot,
                           int nks,
                           int nat,
                           const ModuleBase::matrix& ekb,
                           const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                           const UnitCell& ucell,
                           const LCAO_Orbitals& orb,
                           Grid_Driver& GridD,
                           const Parallel_Orbitals* ParaV,
                           const psi::Psi<std::complex<double>>& psi,
                           const psi::Psi<double>& psid,
                           const std::vector<ModuleBase::matrix>& dm_gamma,
                           const std::vector<ModuleBase::ComplexMatrix>& dm_k);
  private:
    std::shared_ptr<LCAO_DftU_New> ld;
};

#endif