#ifndef ESOLVER_KS_LCAO_H
#define ESOLVER_KS_LCAO_H
#include "esolver_ks.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"
// for grid integration
#include "module_basis/module_ao/ORB_control.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#ifdef __EXX
#include "module_ri/Exx_LRI_interface.h"
#include "module_ri/Mix_DMk_2D.h"
#endif
#include "module_basis/module_nao/two_center_bundle.h"
#include "module_io/output_dm.h"
#include "module_io/output_mat_sparse.h"

#include <memory>

namespace ModuleESolver
{
template <typename TK, typename TR>
class ESolver_KS_LCAO : public ESolver_KS<TK>
{
  public:
    ESolver_KS_LCAO();
    ~ESolver_KS_LCAO();

    void before_all_runners(Input& inp, UnitCell& cell) override;

    void init_after_vc(Input& inp, UnitCell& cell) override;

    double cal_energy() override;

    void cal_force(ModuleBase::matrix& force) override;

    void cal_stress(ModuleBase::matrix& stress) override;

    void after_all_runners() override;

    void nscf() override;

    void get_S();

    void cal_mag(const int istep, const bool print = false);

  protected:
    virtual void before_scf(const int istep) override;

    virtual void iter_init(const int istep, const int iter) override;

    virtual void hamilt2density(const int istep, const int iter, const double ethr) override;

    virtual void update_pot(const int istep, const int iter) override;

    virtual void iter_finish(const int iter) override;

    virtual void after_scf(const int istep) override;

    virtual bool do_after_converge(int& iter) override;

    virtual void others(const int istep) override;

    // we will get rid of this class soon, don't use it, mohan 2024-03-28
    ORB_control orb_con; // Basis_LCAO

    // we will get rid of this class soon, don't use it, mohan 2024-03-28
    Record_adj RA;

    // we will get rid of this class soon, don't use it, mohan 2024-03-28
    Local_Orbital_wfc LOWF;

    // we will get rid of this class soon, don't use it, mohan 2024-03-28
    Local_Orbital_Charge LOC;

    // used for k-dependent grid integration.
    Gint_k GK;

    // used for gamma only algorithms.
    Gint_Gamma GG;

    // we will get rid of this class soon, don't use it, mohan 2024-03-28
    LCAO_Matrix LM;

    Grid_Technique GridT;

    TwoCenterBundle two_center_bundle_;

    // Temporarily store the stress to unify the interface with PW,
    // because it's hard to seperate force and stress calculation in LCAO.
    // The copy costs memory and time !
    // Are there any better way to solve this problem ?
    ModuleBase::matrix scs;
    bool have_force = false;

    void init_basis_lcao(ORB_control& orb_con, Input& inp, UnitCell& ucell);

    //--------------common for all calculation, not only scf-------------
    // set matrix and grid integral
    void set_matrix_grid(Record_adj& ra);

    void beforesolver(const int istep);
    //----------------------------------------------------------------------

    /// @brief create ModuleIO::Output_DM object to output density matrix
    ModuleIO::Output_DM create_Output_DM(int is, int iter);

    /// @brief create ModuleIO::Output_Mat_Sparse object to output sparse density matrix of H, S, T, r
    ModuleIO::Output_Mat_Sparse<TK> create_Output_Mat_Sparse(int istep);

    void read_mat_npz(std::string& zipname, hamilt::HContainer<double>& hR);
    void output_mat_npz(std::string& zipname, const hamilt::HContainer<double>& hR);

    /// @brief check if skip the corresponding output in md calculation
    bool md_skip_out(std::string calculation, int istep, int interval);

#ifdef __EXX
    std::shared_ptr<Exx_LRI_Interface<TK, double>> exd = nullptr;
    std::shared_ptr<Exx_LRI_Interface<TK, std::complex<double>>> exc = nullptr;
    std::shared_ptr<Exx_LRI<double>> exx_lri_double = nullptr;
    std::shared_ptr<Exx_LRI<std::complex<double>>> exx_lri_complex = nullptr;
#endif

  private:
    // tmp interfaces  before sub-modules are refactored
    void dftu_cal_occup_m(const int& iter, const std::vector<std::vector<TK>>& dm) const;

#ifdef __DEEPKS
    void dpks_cal_e_delta_band(const std::vector<std::vector<TK>>& dm) const;

    void dpks_cal_projected_DM(const elecstate::DensityMatrix<TK, double>* dm) const;
#endif
};
} // namespace ModuleESolver
#endif
