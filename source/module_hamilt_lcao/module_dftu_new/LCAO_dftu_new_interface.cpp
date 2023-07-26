#include "LCAO_dftu_new_interface.h"

#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_elecstate/cal_dm.h"

LCAO_DftU_New_Interface::LCAO_DftU_New_Interface(std::shared_ptr<LCAO_DftU_New> ld_in) : ld(ld_in)
{
}

void LCAO_DftU_New_Interface::out_deepks_labels(double etot,
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
                                              const std::vector<ModuleBase::ComplexMatrix>& dm_k)
{
    ModuleBase::TITLE("LCAO_DftU_New_Interface", "out_deepks_labels");                                                    // end deepks_out_labels

    if (GlobalV::deepks_scf)
    {
        // this part is for integrated test of deepks
        // so it is printed no matter even if deepks_out_labels is not used
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            ld->cal_projected_DM(dm_gamma, ucell, orb, GridD);
        }
        else
        {
            ld->cal_projected_DM_k(dm_k, ucell, orb, GridD, nks, kvec_d);
        }
        ld->check_projected_dm(); // print out the projected dm for NSCF calculaiton

        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            ld->cal_e_delta_band(dm_gamma);
        }
        else
        {
            ld->cal_e_delta_band_k(dm_k, nks);
        }
        std::cout << "E_delta_band = " << std::setprecision(8) << ld->e_delta_band << " Ry"
                  << " = " << std::setprecision(8) << ld->e_delta_band * ModuleBase::Ry_to_eV << " eV"
                  << std::endl;
        std::cout << "E_delta_NN= " << std::setprecision(8) << ld->E_delta << " Ry"
                  << " = " << std::setprecision(8) << ld->E_delta * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
}