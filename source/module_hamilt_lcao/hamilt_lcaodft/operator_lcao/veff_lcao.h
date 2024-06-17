#ifndef VEFFLCAO_H
#define VEFFLCAO_H
#include "module_base/timer.h"
#include "module_elecstate/potentials/potential_new.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "operator_lcao.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"

namespace hamilt
{

#ifndef __VEFFTEMPLATE
#define __VEFFTEMPLATE

template <class T>
class Veff : public T
{
};

#endif

/// @brief Effective potential class, used for calculating Hamiltonian with grid integration tools
/// If user want to separate the contribution of V_{eff} into V_{H} and V_{XC} and V_{local pseudopotential} and so on,
/// the user can separate the Potential class into different parts, and construct different Veff class for each part.
/// @tparam TK 
/// @tparam TR 
template <typename TK, typename TR>
class Veff<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    /**
     * @brief Construct a new Veff object for multi-kpoint calculation
     * @param GK_in: the pointer of Gint_k object, used for grid integration
    */
    Veff<OperatorLCAO<TK, TR>>(Gint_k* GK_in,
                          Local_Orbital_Charge* loc_in,
                          LCAO_Matrix* LM_in,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                          elecstate::Potential* pot_in,
                          hamilt::HContainer<TR>* hR_in,
                          std::vector<TK>* hK_in,
                          const UnitCell* ucell_in,
                          Grid_Driver* GridD_in,
                          const Parallel_Orbitals* paraV)
        : GK(GK_in),
          loc(loc_in),
          pot(pot_in),
          ucell(ucell_in),
          gd(GridD_in),
          OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in)
    {
        this->cal_type = lcao_gint;

        this->initialize_HR(ucell_in, GridD_in, paraV);
        GK_in->initialize_pvpR(*ucell_in, GridD_in);
    }
    /**
     * @brief Construct a new Veff object for Gamma-only calculation
     * @param GG_in: the pointer of Gint_Gamma object, used for grid integration
    */
    Veff<OperatorLCAO<TK, TR>>(Gint_Gamma* GG_in,
                          Local_Orbital_Charge* loc_in,
                          LCAO_Matrix* LM_in,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                          elecstate::Potential* pot_in,
                          hamilt::HContainer<TR>* hR_in,
                          std::vector<TK>* hK_in,
                          const UnitCell* ucell_in,
                          Grid_Driver* GridD_in,
                          const Parallel_Orbitals* paraV
                          )
        : GG(GG_in), loc(loc_in), pot(pot_in),
        OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in)
    {
        this->cal_type = lcao_gint;
        this->initialize_HR(ucell_in, GridD_in, paraV);

        GG_in->initialize_pvpR(*ucell_in, GridD_in);
    }

    ~Veff<OperatorLCAO<TK, TR>>(){};

    /**
     * @brief contributeHR() is used to calculate the HR matrix
     * <phi_{\mu, 0}|V_{eff}|phi_{\nu, R}>
     * the contribution of V_{eff} is calculated by the contribution of V_{H} and V_{XC} and V_{local pseudopotential} and so on.
     * grid integration is used to calculate the contribution Hamiltonian of effective potential
     */
    virtual void contributeHR() override;
  
  const UnitCell* ucell;
  Grid_Driver* gd;
  private:
    // used for k-dependent grid integration.
    Gint_k* GK = nullptr;

    // used for gamma only algorithms.
    Gint_Gamma* GG = nullptr;

    // Charge calculating method in LCAO base and contained grid base calculation: DM_R, DM, pvpR_reduced
    Local_Orbital_Charge* loc = nullptr;

    elecstate::Potential* pot = nullptr;

    int nspin = 1;
    int current_spin = 0;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the electronic kinetic matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(const UnitCell* ucell_in, Grid_Driver* GridD_in, const Parallel_Orbitals* paraV);

};

} // namespace hamilt
#endif