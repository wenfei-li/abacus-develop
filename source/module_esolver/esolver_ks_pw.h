#ifndef ESOLVER_KS_PW_H
#define ESOLVER_KS_PW_H
#include "./esolver_ks.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/velocity_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/wfinit.h"

#include <memory>
#include <module_base/macros.h>

namespace ModuleESolver
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class ESolver_KS_PW : public ESolver_KS<T, Device>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    ESolver_KS_PW();

    ~ESolver_KS_PW();

    void before_all_runners(Input& inp, UnitCell& cell) override;

    void init_after_vc(Input& inp, UnitCell& cell) override;

    double cal_energy() override;

    void cal_force(ModuleBase::matrix& force) override;

    void cal_stress(ModuleBase::matrix& stress) override;

    virtual void hamilt2density(const int istep, const int iter, const double ethr) override;

    virtual void hamilt2estates(const double ethr) override;

    virtual void nscf() override;

    void after_all_runners() override;

  protected:
    virtual void before_scf(const int istep) override;

    virtual void iter_init(const int istep, const int iter) override;

    virtual void update_pot(const int istep, const int iter) override;

    virtual void iter_finish(const int iter) override;

    virtual void after_scf(const int istep) override;

    virtual void others(const int istep) override;

    // temporary, this will be removed in the future;
    // Init Global class
    void Init_GlobalC(Input& inp, UnitCell& ucell, pseudopot_cell_vnl& ppcell);

  protected:
    //! hide the psi in ESolver_KS for tmp use
    psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* psi = nullptr;

    // psi_initializer controller
    psi::WFInit<T, Device>* p_wf_init = nullptr;

  private:
    Device* ctx = {};

    base_device::AbacusDevice_t device = {};

    psi::Psi<T, Device>* kspw_psi = nullptr;

    psi::Psi<std::complex<double>, Device>* __kspw_psi = nullptr;

    using castmem_2d_d2h_op
        = base_device::memory::cast_memory_op<std::complex<double>, T, base_device::DEVICE_CPU, Device>;
};
} // namespace ModuleESolver
#endif
