#ifndef ESOLVER_SDFT_PW_H
#define ESOLVER_SDFT_PW_H

#include "esolver_ks_pw.h"
#include "module_hamilt_pw/hamilt_stodft/sto_hchi.h"
#include "module_hamilt_pw/hamilt_stodft/sto_iter.h"
#include "module_hamilt_pw/hamilt_stodft/sto_wf.h"

namespace ModuleESolver
{

class ESolver_SDFT_PW : public ESolver_KS_PW<std::complex<double>>
{
  public:
    ESolver_SDFT_PW();
    ~ESolver_SDFT_PW();

    void before_all_runners(Input& inp, UnitCell& cell) override;

    double cal_energy() override;

    void cal_force(ModuleBase::matrix& force) override;

    void cal_stress(ModuleBase::matrix& stress) override;

  public:
    Stochastic_WF stowf;

  protected:
    virtual void before_scf(const int istep) override;

    virtual void hamilt2density(const int istep, const int iter, const double ethr) override;

    virtual void nscf() override;

    virtual void others(const int istep) override;

    virtual void iter_finish(const int iter) override;

    virtual void after_scf(const int istep) override;

    virtual void after_all_runners() override;

  private:
    int nche_sto;   ///< norder of Chebyshev
    int method_sto; ///< method of SDFT
};

} // namespace ModuleESolver

// temporary setting: removed GlobalC but not breaking design philosophy
namespace GlobalTemp
{

extern const ModuleBase::matrix* veff;

}

#endif
