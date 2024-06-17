#include "elecstate_lcao_tddft.h"

#include "cal_dm.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_base/timer.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
namespace elecstate
{

// multi-k case
void ElecStateLCAO_TDDFT::psiToRho_td(const psi::Psi<std::complex<double>>& psi)
{
    ModuleBase::TITLE("ElecStateLCAO", "psiToRho");
    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");

    this->calculate_weights_td();
    this->calEBand();

    ModuleBase::GlobalFunc::NOTE("Calculate the density matrix.");

    // this part for calculating DMK in 2d-block format, not used for charge now
    //    psi::Psi<std::complex<double>> dm_k_2d();

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx"
        || GlobalV::KS_SOLVER == "lapack") // Peize Lin test 2019-05-15
    {
        elecstate::cal_dm_psi(this->DM->get_paraV_pointer(), this->wg, psi, *(this->DM));
        this->DM->cal_DMR();

// interface for RI-related calculation, which needs loc.dm_k 
#ifdef __EXX
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            const K_Vectors* kv = this->DM->get_kv_pointer();
            this->loc->dm_k.resize(kv->get_nks());
            for (int ik = 0; ik < kv->get_nks(); ++ik)
            {
                this->loc->set_dm_k(ik, this->DM->get_DMK_pointer(ik));         
            }
        }
#endif
    }


    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx); // mohan 2009-11-10
    }

    //------------------------------------------------------------
    // calculate the charge density on real space grid.
    //------------------------------------------------------------

    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
    this->gint_k->transfer_DM2DtoGrid(this->DM->get_DMR_vector()); // transfer DM2D to DM_grid in gint
    Gint_inout inout(this->loc->DM_R, this->charge->rho, Gint_Tools::job_type::rho); // rho calculation
    this->gint_k->cal_gint(&inout);

    this->charge->renormalize_rho();

    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");
    return;
}

void ElecStateLCAO_TDDFT::calculate_weights_td()
{
    ModuleBase::TITLE("ElecState", "calculate_weights");

    if (GlobalV::ocp == 1)
    {
        int num = 0;
        num = this->klist->get_nks() * GlobalV::NBANDS;
        if (num != GlobalV::ocp_kb.size())
        {
            ModuleBase::WARNING_QUIT("ElecStateLCAO_TDDFT::calculate_weights_td",
                                     "size of occupation array is wrong , please check ocp_set");
        }

        double num_elec = 0.0;
        for (int i = 0; i < GlobalV::ocp_kb.size(); i++)
        {
            num_elec += GlobalV::ocp_kb[i];
        }
        if (std::abs(num_elec - GlobalV::nelec) > 1.0e-5)
        {
            ModuleBase::WARNING_QUIT("ElecStateLCAO_TDDFT::calculate_weights_td",
                                     "total number of occupations is wrong , please check ocp_set");
        }

        for (int ik = 0; ik < this->klist->get_nks(); ik++)
        {
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                this->wg(ik, ib) = GlobalV::ocp_kb[ik * GlobalV::NBANDS + ib];
            }
        }
    }
    return;
}


} // namespace elecstate
