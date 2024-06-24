#include "hsolver_pw.h"

#include "diago_bpcg.h"
#include "diago_cg.h"
#include "diago_dav_subspace.h"
#include "diago_david.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_global.h" // for MPI
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_hsolver/diagh.h"
#include "module_hsolver/diago_iter_assist.h"

#include <algorithm>
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif
namespace hsolver
{

template <typename T, typename Device>
HSolverPW<T, Device>::HSolverPW(ModulePW::PW_Basis_K* wfc_basis_in, wavefunc* pwf_in)
{
    this->classname = "HSolverPW";
    this->wfc_basis = wfc_basis_in;
    this->pwf = pwf_in;
    this->diag_ethr = GlobalV::PW_DIAG_THR;
    /*this->init(pbas_in);*/
}

template <typename T, typename Device>
void HSolverPW<T, Device>::initDiagh(const psi::Psi<T, Device>& psi)
{
    if (this->method == "cg")
    {
        // if (this->pdiagh != nullptr)
        // {
        //     if (this->pdiagh->method != this->method)
        //     {
        //         delete reinterpret_cast<DiagoCG<T, Device>*>(this->pdiagh);
        //     }
        //     else
        //     {
        //         return;
        //     }
        // }

        // this->pdiagh = new DiagoCG<T, Device>(precondition.data());

        // // warp the subspace_func into a lambda function
        // auto ngk_pointer = psi.get_ngk_pointer();
        // auto subspace_func = [this, ngk_pointer](const ct::Tensor& psi_in, ct::Tensor& psi_out) {
        //     // psi_in should be a 2D tensor:
        //     // psi_in.shape() = [nbands, nbasis]
        //     const auto ndim = psi_in.shape().ndim();
        //     REQUIRES_OK(ndim == 2, "dims of psi_in should be less than or equal to 2");
        //     // Convert a Tensor object to a psi::Psi object
        //     auto psi_in_wrapper = psi::Psi<T, Device>(psi_in.data<T>(),
        //                                               1,
        //                                               psi_in.shape().dim_size(0),
        //                                               psi_in.shape().dim_size(1),
        //                                               ngk_pointer);
        //     auto psi_out_wrapper = psi::Psi<T, Device>(psi_out.data<T>(),
        //                                                1,
        //                                                psi_out.shape().dim_size(0),
        //                                                psi_out.shape().dim_size(1),
        //                                                ngk_pointer);
        //     auto eigen = ct::Tensor(ct::DataTypeToEnum<Real>::value,
        //                             ct::DeviceType::CpuDevice,
        //                             ct::TensorShape({psi_in.shape().dim_size(0)}));

        //     DiagoIterAssist<T, Device>::diagH_subspace(hamilt_, psi_in_wrapper, psi_out_wrapper, eigen.data<Real>());
        // };
        // this->pdiagh = new DiagoCG<T, Device>(GlobalV::BASIS_TYPE,
        //                                       GlobalV::CALCULATION,
        //                                       DiagoIterAssist<T, Device>::need_subspace,
        //                                       subspace_func,
        //                                       DiagoIterAssist<T, Device>::PW_DIAG_THR,
        //                                       DiagoIterAssist<T, Device>::PW_DIAG_NMAX,
        //                                       GlobalV::NPROC_IN_POOL);
        // this->pdiagh->method = this->method;
    }
    else if (this->method == "dav")
    {
        // #ifdef __MPI
        //         const diag_comm_info comm_info = {POOL_WORLD, GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL};
        // #else
        //         const diag_comm_info comm_info = {GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL};
        // #endif

        //         if (this->pdiagh != nullptr)
        //         {
        //             if (this->pdiagh->method != this->method)
        //             {
        //                 delete (DiagoDavid<T, Device>*)this->pdiagh;

        //                 this->pdiagh = new DiagoDavid<T, Device>(precondition.data(),
        //                                                          GlobalV::PW_DIAG_NDIM,
        //                                                          GlobalV::use_paw,
        //                                                          comm_info);

        //                 this->pdiagh->method = this->method;
        //             }
        //         }
        //         else
        //         {
        //             this->pdiagh
        //                 = new DiagoDavid<T, Device>(precondition.data(), GlobalV::PW_DIAG_NDIM, GlobalV::use_paw,
        //                 comm_info);

        //             this->pdiagh->method = this->method;
        //         }
    }
    else if (this->method == "dav_subspace")
    {
        // #ifdef __MPI
        //         const diag_comm_info comm_info = {POOL_WORLD, GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL};
        // #else
        //         const diag_comm_info comm_info = {GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL};
        // #endif
        //         if (this->pdiagh != nullptr)
        //         {
        //             if (this->pdiagh->method != this->method)
        //             {
        //                 delete (Diago_DavSubspace<T, Device>*)this->pdiagh;

        //                 this->pdiagh = new Diago_DavSubspace<T, Device>(precondition.data(),
        //                                                                 GlobalV::PW_DIAG_NDIM,
        //                                                                 DiagoIterAssist<T, Device>::PW_DIAG_THR,
        //                                                                 DiagoIterAssist<T, Device>::PW_DIAG_NMAX,
        //                                                                 DiagoIterAssist<T, Device>::need_subspace,
        //                                                                 comm_info);

        //                 this->pdiagh->method = this->method;
        //             }
        //         }
        //         else
        //         {
        //             this->pdiagh = new Diago_DavSubspace<T, Device>(precondition.data(),
        //                                                             GlobalV::PW_DIAG_NDIM,
        //                                                             DiagoIterAssist<T, Device>::PW_DIAG_THR,
        //                                                             DiagoIterAssist<T, Device>::PW_DIAG_NMAX,
        //                                                             DiagoIterAssist<T, Device>::need_subspace,
        //                                                             comm_info);
        //             this->pdiagh->method = this->method;
        //         }
    }
    else if (this->method == "bpcg")
    {
        // if (this->pdiagh != nullptr)
        // {
        //     if (this->pdiagh->method != this->method)
        //     {
        //         delete (DiagoBPCG<T, Device>*)this->pdiagh;
        //         this->pdiagh = new DiagoBPCG<T, Device>(precondition.data());
        //         this->pdiagh->method = this->method;
        //         reinterpret_cast<DiagoBPCG<T, Device>*>(this->pdiagh)->init_iter(psi);
        //     }
        // }
        // else
        // {
        //     this->pdiagh = new DiagoBPCG<T, Device>(precondition.data());
        //     this->pdiagh->method = this->method;
        //     reinterpret_cast<DiagoBPCG<T, Device>*>(this->pdiagh)->init_iter(psi);
        // }
    }
    else
    {
        ModuleBase::WARNING_QUIT("HSolverPW::solve", "This method of DiagH is not supported!");
    }
}

template <typename T, typename Device>
void HSolverPW<T, Device>::solve(hamilt::Hamilt<T, Device>* pHamilt,
                                 psi::Psi<T, Device>& psi,
                                 elecstate::ElecState* pes,
                                 const std::string method_in,
                                 const bool skip_charge)
{
    ModuleBase::TITLE("HSolverPW", "solve");
    ModuleBase::timer::tick("HSolverPW", "solve");
    // prepare for the precondition of diagonalization
    this->precondition.resize(psi.get_nbasis());
    this->hamilt_ = pHamilt;
    // select the method of diagonalization
    this->method = method_in;
    this->initDiagh(psi);

    std::vector<Real> eigenvalues(pes->ekb.nr * pes->ekb.nc, 0);

    if (this->is_first_scf)
    {
        is_occupied.resize(psi.get_nk() * psi.get_nbands(), true);
    }
    else
    {
        if (this->diago_full_acc)
        {
            is_occupied.assign(is_occupied.size(), true);
        }
        else
        {
            for (int i = 0; i < psi.get_nk(); i++)
            {
                if (pes->klist->wk[i] > 0.0)
                {
                    for (int j = 0; j < psi.get_nbands(); j++)
                    {
                        if (pes->wg(i, j) / pes->klist->wk[i] < 0.01)
                        {
                            is_occupied[i * psi.get_nbands() + j] = false;
                        }
                    }
                }
            }
        }
    }

    /// Loop over k points for solve Hamiltonian to charge density
    for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
    {
        /// update H(k) for each k point
        pHamilt->updateHk(ik);

#ifdef USE_PAW
        if (GlobalV::use_paw)
        {
            const int npw = this->wfc_basis->npwk[ik];
            ModuleBase::Vector3<double>* _gk = new ModuleBase::Vector3<double>[npw];
            for (int ig = 0; ig < npw; ig++)
            {
                _gk[ig] = this->wfc_basis->getgpluskcar(ik, ig);
            }

            std::vector<double> kpt(3, 0);
            kpt[0] = this->wfc_basis->kvec_c[ik].x;
            kpt[1] = this->wfc_basis->kvec_c[ik].y;
            kpt[2] = this->wfc_basis->kvec_c[ik].z;

            double** kpg;
            double** gcar;
            kpg = new double*[npw];
            gcar = new double*[npw];
            for (int ipw = 0; ipw < npw; ipw++)
            {
                kpg[ipw] = new double[3];
                kpg[ipw][0] = _gk[ipw].x;
                kpg[ipw][1] = _gk[ipw].y;
                kpg[ipw][2] = _gk[ipw].z;

                gcar[ipw] = new double[3];
                gcar[ipw][0] = this->wfc_basis->getgcar(ik, ipw).x;
                gcar[ipw][1] = this->wfc_basis->getgcar(ik, ipw).y;
                gcar[ipw][2] = this->wfc_basis->getgcar(ik, ipw).z;
            }

            GlobalC::paw_cell.set_paw_k(npw,
                                        wfc_basis->npwk_max,
                                        kpt.data(),
                                        this->wfc_basis->get_ig2ix(ik).data(),
                                        this->wfc_basis->get_ig2iy(ik).data(),
                                        this->wfc_basis->get_ig2iz(ik).data(),
                                        (const double**)kpg,
                                        GlobalC::ucell.tpiba,
                                        (const double**)gcar);

            std::vector<double>().swap(kpt);
            for (int ipw = 0; ipw < npw; ipw++)
            {
                delete[] kpg[ipw];
                delete[] gcar[ipw];
            }
            delete[] kpg;
            delete[] gcar;

            GlobalC::paw_cell.get_vkb();

            GlobalC::paw_cell.set_currentk(ik);
        }
#endif

        this->updatePsiK(pHamilt, psi, ik);

        // template add precondition calculating here
        update_precondition(precondition, ik, this->wfc_basis->npwk[ik]);

#ifdef USE_PAW
        GlobalC::paw_cell.set_currentk(ik);
#endif

        /// solve eigenvector and eigenvalue for H(k)
        this->hamiltSolvePsiK(pHamilt, psi, eigenvalues.data() + ik * pes->ekb.nc);

        if (skip_charge)
        {
            GlobalV::ofs_running << "Average iterative diagonalization steps for k-points " << ik
                                 << " is: " << DiagoIterAssist<T, Device>::avg_iter
                                 << " ; where current threshold is: " << DiagoIterAssist<T, Device>::PW_DIAG_THR
                                 << " . " << std::endl;
            DiagoIterAssist<T, Device>::avg_iter = 0.0;
        }
        /// calculate the contribution of Psi for charge density rho
    }
    // END Loop over k points

    castmem_2d_2h_op()(cpu_ctx, cpu_ctx, pes->ekb.c, eigenvalues.data(), pes->ekb.nr * pes->ekb.nc);

    this->is_first_scf = false;

    this->endDiagh();

    if (skip_charge)
    {
        ModuleBase::timer::tick("HSolverPW", "solve");
        return;
    }
    reinterpret_cast<elecstate::ElecStatePW<T, Device>*>(pes)->psiToRho(psi);

#ifdef USE_PAW
    if (GlobalV::use_paw)
    {
        if (typeid(Real) != typeid(double))
        {
            ModuleBase::WARNING_QUIT("HSolverPW::solve", "PAW is only supported for double precision!");
        }

        GlobalC::paw_cell.reset_rhoij();
        for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
        {
            const int npw = this->wfc_basis->npwk[ik];
            ModuleBase::Vector3<double>* _gk = new ModuleBase::Vector3<double>[npw];
            for (int ig = 0; ig < npw; ig++)
            {
                _gk[ig] = this->wfc_basis->getgpluskcar(ik, ig);
            }

            std::vector<double> kpt(3, 0);
            kpt[0] = this->wfc_basis->kvec_c[ik].x;
            kpt[1] = this->wfc_basis->kvec_c[ik].y;
            kpt[2] = this->wfc_basis->kvec_c[ik].z;

            double** kpg;
            double** gcar;
            kpg = new double*[npw];
            gcar = new double*[npw];
            for (int ipw = 0; ipw < npw; ipw++)
            {
                kpg[ipw] = new double[3];
                kpg[ipw][0] = _gk[ipw].x;
                kpg[ipw][1] = _gk[ipw].y;
                kpg[ipw][2] = _gk[ipw].z;

                gcar[ipw] = new double[3];
                gcar[ipw][0] = this->wfc_basis->getgcar(ik, ipw).x;
                gcar[ipw][1] = this->wfc_basis->getgcar(ik, ipw).y;
                gcar[ipw][2] = this->wfc_basis->getgcar(ik, ipw).z;
            }

            GlobalC::paw_cell.set_paw_k(npw,
                                        wfc_basis->npwk_max,
                                        kpt.data(),
                                        this->wfc_basis->get_ig2ix(ik).data(),
                                        this->wfc_basis->get_ig2iy(ik).data(),
                                        this->wfc_basis->get_ig2iz(ik).data(),
                                        (const double**)kpg,
                                        GlobalC::ucell.tpiba,
                                        (const double**)gcar);

            std::vector<double>().swap(kpt);
            for (int ipw = 0; ipw < npw; ipw++)
            {
                delete[] kpg[ipw];
                delete[] gcar[ipw];
            }
            delete[] kpg;
            delete[] gcar;

            GlobalC::paw_cell.get_vkb();

            psi.fix_k(ik);
            GlobalC::paw_cell.set_currentk(ik);
            int nbands = psi.get_nbands();
            for (int ib = 0; ib < nbands; ib++)
            {
                GlobalC::paw_cell.accumulate_rhoij(reinterpret_cast<std::complex<double>*>(psi.get_pointer(ib)),
                                                   pes->wg(ik, ib));
            }
        }

        std::vector<std::vector<double>> rhoijp;
        std::vector<std::vector<int>> rhoijselect;
        std::vector<int> nrhoijsel;

#ifdef __MPI
        if (GlobalV::RANK_IN_POOL == 0)
        {
            GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

            for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
            {
                GlobalC::paw_cell.set_rhoij(iat,
                                            nrhoijsel[iat],
                                            rhoijselect[iat].size(),
                                            rhoijselect[iat].data(),
                                            rhoijp[iat].data());
            }
        }
#else
        GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

        for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
        {
            GlobalC::paw_cell.set_rhoij(iat,
                                        nrhoijsel[iat],
                                        rhoijselect[iat].size(),
                                        rhoijselect[iat].data(),
                                        rhoijp[iat].data());
        }

#endif
        double* nhatgr;
        GlobalC::paw_cell.get_nhat(pes->charge->nhat, nhatgr);
    }
#endif
    ModuleBase::timer::tick("HSolverPW", "solve");
    return;
}

/*
    lcao_in_pw
*/
template <typename T, typename Device>
void HSolverPW<T, Device>::solve(hamilt::Hamilt<T, Device>* pHamilt, // ESolver_KS_PW::p_hamilt
                                 psi::Psi<T, Device>& psi,           // ESolver_KS_PW::kspw_psi
                                 elecstate::ElecState* pes,          // ESolver_KS_PW::pes
                                 psi::Psi<T, Device>& transform,
                                 const bool skip_charge)
{
    ModuleBase::TITLE("HSolverPW", "solve");
    ModuleBase::timer::tick("HSolverPW", "solve");
    std::vector<Real> eigenvalues(pes->ekb.nr * pes->ekb.nc, 0);
    for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
    {
        /// update H(k) for each k point
        pHamilt->updateHk(ik);
#ifdef USE_PAW
        if (GlobalV::use_paw)
        {
            const int npw = this->wfc_basis->npwk[ik];
            ModuleBase::Vector3<double>* _gk = new ModuleBase::Vector3<double>[npw];
            for (int ig = 0; ig < npw; ig++)
            {
                _gk[ig] = this->wfc_basis->getgpluskcar(ik, ig);
            }

            std::vector<double> kpt(3, 0);
            kpt[0] = this->wfc_basis->kvec_c[ik].x;
            kpt[1] = this->wfc_basis->kvec_c[ik].y;
            kpt[2] = this->wfc_basis->kvec_c[ik].z;

            double** kpg;
            double** gcar;
            kpg = new double*[npw];
            gcar = new double*[npw];
            for (int ipw = 0; ipw < npw; ipw++)
            {
                kpg[ipw] = new double[3];
                kpg[ipw][0] = _gk[ipw].x;
                kpg[ipw][1] = _gk[ipw].y;
                kpg[ipw][2] = _gk[ipw].z;

                gcar[ipw] = new double[3];
                gcar[ipw][0] = this->wfc_basis->getgcar(ik, ipw).x;
                gcar[ipw][1] = this->wfc_basis->getgcar(ik, ipw).y;
                gcar[ipw][2] = this->wfc_basis->getgcar(ik, ipw).z;
            }

            GlobalC::paw_cell.set_paw_k(npw,
                                        wfc_basis->npwk_max,
                                        kpt.data(),
                                        this->wfc_basis->get_ig2ix(ik).data(),
                                        this->wfc_basis->get_ig2iy(ik).data(),
                                        this->wfc_basis->get_ig2iz(ik).data(),
                                        (const double**)kpg,
                                        GlobalC::ucell.tpiba,
                                        (const double**)gcar);

            std::vector<double>().swap(kpt);
            for (int ipw = 0; ipw < npw; ipw++)
            {
                delete[] kpg[ipw];
                delete[] gcar[ipw];
            }
            delete[] kpg;
            delete[] gcar;

            GlobalC::paw_cell.get_vkb();

            GlobalC::paw_cell.set_currentk(ik);
        }
#endif
        psi.fix_k(ik);
        transform.fix_k(ik);
        /// solve eigenvector and eigenvalue for H(k)

        // hsolver::DiagoIterAssist<T, Device>::diagH_subspace(
        //     pHamilt,                                // interface to hamilt
        //     transform,                              // transform matrix between lcao and pw
        //     psi,                                    // psi in pw basis
        //     eigenvalues.data() + ik * pes->ekb.nc,  // eigenvalues
        //     psi.get_nbands()                        // number of the lowest energies bands
        //     );

        hsolver::DiagoIterAssist<T, Device>::diagH_subspace_init(
            pHamilt,                 // interface to hamilt
            transform.get_pointer(), // transform matrix between lcao and pw
            transform.get_nbands(),
            transform.get_nbasis(),
            psi,                                  // psi in pw basis
            eigenvalues.data() + ik * pes->ekb.nc // eigenvalues
        );

        if (skip_charge)
        {
            GlobalV::ofs_running << "Average iterative diagonalization steps for k-points " << ik
                                 << " is: " << DiagoIterAssist<T, Device>::avg_iter
                                 << " ; where current threshold is: " << DiagoIterAssist<T, Device>::PW_DIAG_THR
                                 << " . " << std::endl;
            DiagoIterAssist<T, Device>::avg_iter = 0.0;
        }
        /// calculate the contribution of Psi for charge density rho
    }
    castmem_2d_2h_op()(cpu_ctx, cpu_ctx, pes->ekb.c, eigenvalues.data(), pes->ekb.nr * pes->ekb.nc);

    if (skip_charge)
    {
        ModuleBase::timer::tick("HSolverPW", "solve");
        return;
    }
    reinterpret_cast<elecstate::ElecStatePW<T, Device>*>(pes)->psiToRho(psi);

#ifdef USE_PAW
    if (GlobalV::use_paw)
    {
        if (typeid(Real) != typeid(double))
        {
            ModuleBase::WARNING_QUIT("HSolverPW::solve", "PAW is only supported for double precision!");
        }

        GlobalC::paw_cell.reset_rhoij();
        for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
        {
            const int npw = this->wfc_basis->npwk[ik];
            ModuleBase::Vector3<double>* _gk = new ModuleBase::Vector3<double>[npw];
            for (int ig = 0; ig < npw; ig++)
            {
                _gk[ig] = this->wfc_basis->getgpluskcar(ik, ig);
            }

            std::vector<double> kpt(3, 0);
            kpt[0] = this->wfc_basis->kvec_c[ik].x;
            kpt[1] = this->wfc_basis->kvec_c[ik].y;
            kpt[2] = this->wfc_basis->kvec_c[ik].z;

            double** kpg;
            double** gcar;
            kpg = new double*[npw];
            gcar = new double*[npw];
            for (int ipw = 0; ipw < npw; ipw++)
            {
                kpg[ipw] = new double[3];
                kpg[ipw][0] = _gk[ipw].x;
                kpg[ipw][1] = _gk[ipw].y;
                kpg[ipw][2] = _gk[ipw].z;

                gcar[ipw] = new double[3];
                gcar[ipw][0] = this->wfc_basis->getgcar(ik, ipw).x;
                gcar[ipw][1] = this->wfc_basis->getgcar(ik, ipw).y;
                gcar[ipw][2] = this->wfc_basis->getgcar(ik, ipw).z;
            }

            GlobalC::paw_cell.set_paw_k(npw,
                                        wfc_basis->npwk_max,
                                        kpt.data(),
                                        this->wfc_basis->get_ig2ix(ik).data(),
                                        this->wfc_basis->get_ig2iy(ik).data(),
                                        this->wfc_basis->get_ig2iz(ik).data(),
                                        (const double**)kpg,
                                        GlobalC::ucell.tpiba,
                                        (const double**)gcar);

            std::vector<double>().swap(kpt);
            for (int ipw = 0; ipw < npw; ipw++)
            {
                delete[] kpg[ipw];
                delete[] gcar[ipw];
            }
            delete[] kpg;
            delete[] gcar;

            GlobalC::paw_cell.get_vkb();

            psi.fix_k(ik);
            GlobalC::paw_cell.set_currentk(ik);
            int nbands = psi.get_nbands();
            for (int ib = 0; ib < nbands; ib++)
            {
                GlobalC::paw_cell.accumulate_rhoij(reinterpret_cast<std::complex<double>*>(psi.get_pointer(ib)),
                                                   pes->wg(ik, ib));
            }
        }

        std::vector<std::vector<double>> rhoijp;
        std::vector<std::vector<int>> rhoijselect;
        std::vector<int> nrhoijsel;

#ifdef __MPI
        if (GlobalV::RANK_IN_POOL == 0)
        {
            GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

            for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
            {
                GlobalC::paw_cell.set_rhoij(iat,
                                            nrhoijsel[iat],
                                            rhoijselect[iat].size(),
                                            rhoijselect[iat].data(),
                                            rhoijp[iat].data());
            }
        }
#else
        GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

        for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
        {
            GlobalC::paw_cell.set_rhoij(iat,
                                        nrhoijsel[iat],
                                        rhoijselect[iat].size(),
                                        rhoijselect[iat].data(),
                                        rhoijp[iat].data());
        }

#endif
        double* nhatgr;
        GlobalC::paw_cell.get_nhat(pes->charge->nhat, nhatgr);
    }
#endif
    ModuleBase::timer::tick("HSolverPW", "solve");
    return;
}

template <typename T, typename Device>
void HSolverPW<T, Device>::endDiagh()
{
    // DiagoCG would keep 9*nbasis memory in cache during loop-k
    // it should be deleted before calculating charge
    // if (this->method == "cg")
    // {
    //     delete reinterpret_cast<DiagoCG<T, Device>*>(this->pdiagh);
    //     this->pdiagh = nullptr;
    // }
    // if (this->method == "dav")
    // {
    //     delete reinterpret_cast<DiagoDavid<T, Device>*>(this->pdiagh);
    //     this->pdiagh = nullptr;
    // }
    // if (this->method == "dav_subspace")
    // {
    //     delete reinterpret_cast<Diago_DavSubspace<T, Device>*>(this->pdiagh);
    //     this->pdiagh = nullptr;
    // }
    // if (this->method == "bpcg")
    // {
    //     delete reinterpret_cast<DiagoBPCG<T, Device>*>(this->pdiagh);
    //     this->pdiagh = nullptr;
    // }

    // in PW base, average iteration steps for each band and k-point should be printing
    if (DiagoIterAssist<T, Device>::avg_iter > 0.0)
    {
        GlobalV::ofs_running << "Average iterative diagonalization steps: "
                             << DiagoIterAssist<T, Device>::avg_iter / this->wfc_basis->nks
                             << " ; where current threshold is: " << DiagoIterAssist<T, Device>::PW_DIAG_THR << " . "
                             << std::endl;

        // std::cout << "avg_iter == " << DiagoIterAssist<T, Device>::avg_iter << std::endl;

        // reset avg_iter
        DiagoIterAssist<T, Device>::avg_iter = 0.0;
    }
    // psi only should be initialed once for PW
    if (!this->initialed_psi)
    {
        this->initialed_psi = true;
    }
}

template <typename T, typename Device>
void HSolverPW<T, Device>::updatePsiK(hamilt::Hamilt<T, Device>* pHamilt, psi::Psi<T, Device>& psi, const int ik)
{
    psi.fix_k(ik);
    if (GlobalV::psi_initializer) // new psi initialization method branch
    {
        // do nothing here, because we have already initialize, allocate and make initial guess
        // basis_type lcao_in_pw function may be inserted here
    }
    else if (!this->initialed_psi) // old psi initialization method branch
    {
        if (GlobalV::BASIS_TYPE == "pw")
        {
            hamilt::diago_PAO_in_pw_k2(this->ctx, ik, psi, this->wfc_basis, this->pwf, pHamilt);
        }
        /* lcao_in_pw now is based on newly implemented psi initializer, so it does not appear here*/
    }
}

template <typename T, typename Device>
void HSolverPW<T, Device>::hamiltSolvePsiK(hamilt::Hamilt<T, Device>* hm, psi::Psi<T, Device>& psi, Real* eigenvalue)
{
    if (this->method == "cg")
    {
        // warp the subspace_func into a lambda function
        auto ngk_pointer = psi.get_ngk_pointer();
        auto subspace_func = [this, ngk_pointer](const ct::Tensor& psi_in, ct::Tensor& psi_out) {
            // psi_in should be a 2D tensor:
            // psi_in.shape() = [nbands, nbasis]
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim == 2, "dims of psi_in should be less than or equal to 2");
            // Convert a Tensor object to a psi::Psi object
            auto psi_in_wrapper = psi::Psi<T, Device>(psi_in.data<T>(),
                                                      1,
                                                      psi_in.shape().dim_size(0),
                                                      psi_in.shape().dim_size(1),
                                                      ngk_pointer);
            auto psi_out_wrapper = psi::Psi<T, Device>(psi_out.data<T>(),
                                                       1,
                                                       psi_out.shape().dim_size(0),
                                                       psi_out.shape().dim_size(1),
                                                       ngk_pointer);
            auto eigen = ct::Tensor(ct::DataTypeToEnum<Real>::value,
                                    ct::DeviceType::CpuDevice,
                                    ct::TensorShape({psi_in.shape().dim_size(0)}));

            DiagoIterAssist<T, Device>::diagH_subspace(hamilt_, psi_in_wrapper, psi_out_wrapper, eigen.data<Real>());
        };
        DiagoCG<T, Device> cg(GlobalV::BASIS_TYPE,
                              GlobalV::CALCULATION,
                              DiagoIterAssist<T, Device>::need_subspace,
                              subspace_func,
                              DiagoIterAssist<T, Device>::PW_DIAG_THR,
                              DiagoIterAssist<T, Device>::PW_DIAG_NMAX,
                              GlobalV::NPROC_IN_POOL);

        // warp the hpsi_func and spsi_func into a lambda function
        using ct_Device = typename ct::PsiToContainer<Device>::type;

        // warp the hpsi_func and spsi_func into a lambda function
        auto hpsi_func = [hm, ngk_pointer](const ct::Tensor& psi_in, ct::Tensor& hpsi_out) {
            ModuleBase::timer::tick("DiagoCG_New", "hpsi_func");
            // psi_in should be a 2D tensor:
            // psi_in.shape() = [nbands, nbasis]
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");
            // Convert a Tensor object to a psi::Psi object
            auto psi_wrapper = psi::Psi<T, Device>(psi_in.data<T>(),
                                                   1,
                                                   ndim == 1 ? 1 : psi_in.shape().dim_size(0),
                                                   ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1),
                                                   ngk_pointer);
            psi::Range all_bands_range(true, psi_wrapper.get_current_k(), 0, psi_wrapper.get_nbands() - 1);
            using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;
            hpsi_info info(&psi_wrapper, all_bands_range, hpsi_out.data<T>());
            hm->ops->hPsi(info);
            ModuleBase::timer::tick("DiagoCG_New", "hpsi_func");
        };
        auto spsi_func = [this, hm](const ct::Tensor& psi_in, ct::Tensor& spsi_out) {
            ModuleBase::timer::tick("DiagoCG_New", "spsi_func");
            // psi_in should be a 2D tensor:
            // psi_in.shape() = [nbands, nbasis]
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");

            if (GlobalV::use_uspp)
            {
                // Convert a Tensor object to a psi::Psi object
                hm->sPsi(psi_in.data<T>(),
                         spsi_out.data<T>(),
                         ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1),
                         ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1),
                         ndim == 1 ? 1 : psi_in.shape().dim_size(0));
            }
            else
            {
                base_device::memory::synchronize_memory_op<T, Device, Device>()(
                    this->ctx,
                    this->ctx,
                    spsi_out.data<T>(),
                    psi_in.data<T>(),
                    static_cast<size_t>((ndim == 1 ? 1 : psi_in.shape().dim_size(0))
                                        * (ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1))));
            }

            ModuleBase::timer::tick("DiagoCG_New", "spsi_func");
        };
        auto psi_tensor = ct::TensorMap(psi.get_pointer(),
                                        ct::DataTypeToEnum<T>::value,
                                        ct::DeviceTypeToEnum<ct_Device>::value,
                                        ct::TensorShape({psi.get_nbands(), psi.get_nbasis()}))
                              .slice({0, 0}, {psi.get_nbands(), psi.get_current_nbas()});
        auto eigen_tensor = ct::TensorMap(eigenvalue,
                                          ct::DataTypeToEnum<Real>::value,
                                          ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                                          ct::TensorShape({psi.get_nbands()}));
        auto prec_tensor = ct::TensorMap(precondition.data(),
                                         ct::DataTypeToEnum<Real>::value,
                                         ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                                         ct::TensorShape({static_cast<int>(precondition.size())}))
                               .to_device<ct_Device>()
                               .slice({0}, {psi.get_current_nbas()});

        cg.diag(hpsi_func, spsi_func, psi_tensor, eigen_tensor, prec_tensor);
        // TODO: Double check tensormap's potential problem
        ct::TensorMap(psi.get_pointer(), psi_tensor, {psi.get_nbands(), psi.get_nbasis()}).sync(psi_tensor);
    }
    else if (this->method == "dav_subspace")
    {
#ifdef __MPI
        const diag_comm_info comm_info = {POOL_WORLD, GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL};
#else
        const diag_comm_info comm_info = {GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL};
#endif
        Diago_DavSubspace<T, Device> dav_subspace(this->precondition,
                                                  psi.get_nbands(),
                                                  psi.get_k_first() ? psi.get_current_nbas()
                                                                    : psi.get_nk() * psi.get_nbasis(),
                                                  GlobalV::PW_DIAG_NDIM,
                                                  DiagoIterAssist<T, Device>::PW_DIAG_THR,
                                                  DiagoIterAssist<T, Device>::PW_DIAG_NMAX,
                                                  DiagoIterAssist<T, Device>::need_subspace,
                                                  comm_info);

        bool scf;
        if (GlobalV::CALCULATION == "nscf")
        {
            scf = false;
        }
        else
        {
            scf = true;
        }

        auto ngk_pointer = psi.get_ngk_pointer();

        auto hpsi_func = [hm, ngk_pointer](T* hpsi_out,
                                           T* psi_in,
                                           const int nband_in,
                                           const int nbasis_in,
                                           const int band_index1,
                                           const int band_index2) {
            ModuleBase::timer::tick("DavSubspace", "hpsi_func");

            // Convert "pointer data stucture" to a psi::Psi object
            auto psi_iter_wrapper = psi::Psi<T, Device>(psi_in, 1, nband_in, nbasis_in, ngk_pointer);

            psi::Range bands_range(1, 0, band_index1, band_index2);

            using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;
            hpsi_info info(&psi_iter_wrapper, bands_range, hpsi_out);
            hm->ops->hPsi(info);

            ModuleBase::timer::tick("DavSubspace", "hpsi_func");
        };

        auto subspace_func = [hm, ngk_pointer](T* psi_out,
                                               T* psi_in,
                                               Real* eigenvalue_in_hsolver,
                                               const int nband_in,
                                               const int nbasis_max_in) {
            // Convert "pointer data stucture" to a psi::Psi object
            auto psi_in_wrapper = psi::Psi<T, Device>(psi_in, 1, nband_in, nbasis_max_in, ngk_pointer);
            auto psi_out_wrapper = psi::Psi<T, Device>(psi_out, 1, nband_in, nbasis_max_in, ngk_pointer);

            DiagoIterAssist<T, Device>::diagH_subspace(hm,
                                                       psi_in_wrapper,
                                                       psi_out_wrapper,
                                                       eigenvalue_in_hsolver,
                                                       nband_in);
        };

        DiagoIterAssist<T, Device>::avg_iter += static_cast<double>(
            dav_subspace
                .diag(hpsi_func, subspace_func, psi.get_pointer(), psi.get_nbasis(), eigenvalue, is_occupied, scf));

        this->pdiagh = nullptr;
    }
    else if (this->method == "bpcg")
    {
        // this->pdiagh->diag(hm, psi, eigenvalue);
        DiagoBPCG<T, Device> bpcg(precondition.data());
        bpcg.init_iter(psi);
        bpcg.diag(hm, psi, eigenvalue);
    }
    else if (this->method == "dav")
    {
#ifdef __MPI
        const diag_comm_info comm_info = {POOL_WORLD, GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL};
#else
        const diag_comm_info comm_info = {GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL};
#endif
        // Allow 5 tries at most. If ntry > ntry_max = 5, exit diag loop.
        const int ntry_max = 5;
        // In non-self consistent calculation, do until totally converged. Else allow 5 eigenvecs to be NOT
        // converged.
        const int notconv_max = ("nscf" == GlobalV::CALCULATION) ? 0 : 5;
        // do diag and add davidson iteration counts up to avg_iter
        const Real david_diag_thr = DiagoIterAssist<T, Device>::PW_DIAG_THR;
        const int david_maxiter = DiagoIterAssist<T, Device>::PW_DIAG_NMAX;

        DiagoDavid<T, Device> david(precondition.data(), GlobalV::PW_DIAG_NDIM, GlobalV::use_paw, comm_info);
        DiagoIterAssist<T, Device>::avg_iter += static_cast<double>(
            david.diag(hm, psi, eigenvalue, david_diag_thr, david_maxiter, ntry_max, notconv_max));
    }
    return;
}

template <typename T, typename Device>
void HSolverPW<T, Device>::update_precondition(std::vector<Real>& h_diag, const int ik, const int npw)
{
    h_diag.assign(h_diag.size(), 1.0);
    int precondition_type = 2;
    const auto tpiba2 = static_cast<Real>(this->wfc_basis->tpiba2);

    //===========================================
    // Conjugate-Gradient diagonalization
    // h_diag is the precondition matrix
    // h_diag(1:npw) = MAX( 1.0, g2kin(1:npw) );
    //===========================================
    if (precondition_type == 1)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            Real g2kin = static_cast<Real>(this->wfc_basis->getgk2(ik, ig)) * tpiba2;
            h_diag[ig] = std::max(static_cast<Real>(1.0), g2kin);
        }
    }
    else if (precondition_type == 2)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            Real g2kin = static_cast<Real>(this->wfc_basis->getgk2(ik, ig)) * tpiba2;

            if (this->method == "dav_subspace")
            {
                h_diag[ig] = g2kin;
            }
            else
            {
                h_diag[ig] = 1 + g2kin + sqrt(1 + (g2kin - 1) * (g2kin - 1));
            }
        }
    }
    if (GlobalV::NSPIN == 4)
    {
        const int size = h_diag.size();
        for (int ig = 0; ig < npw; ig++)
        {
            h_diag[ig + size / 2] = h_diag[ig];
        }
    }
}

template <typename T, typename Device>
typename HSolverPW<T, Device>::Real HSolverPW<T, Device>::cal_hsolerror()
{
    return this->diag_ethr * static_cast<Real>(std::max(1.0, GlobalV::nelec));
}

template <typename T, typename Device>
typename HSolverPW<T, Device>::Real HSolverPW<T, Device>::set_diagethr(const int istep, const int iter, const Real drho)
{
    // It is too complex now and should be modified.
    if (iter == 1)
    {
        if (std::abs(this->diag_ethr - 1.0e-2) < 1.0e-6)
        {
            if (GlobalV::init_chg == "file")
            {
                //======================================================
                // if you think that the starting potential is good
                // do not spoil it with a louly first diagonalization:
                // set a strict this->diag_ethr in the input file ()diago_the_init
                //======================================================
                this->diag_ethr = 1.0e-5;
            }
            else
            {
                //=======================================================
                // starting atomic potential is probably far from scf
                // don't waste iterations in the first diagonalization
                //=======================================================
                this->diag_ethr = 1.0e-2;
            }
        }
        // if (GlobalV::FINAL_SCF) this->diag_ethr = 1.0e-2;
        if (GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
        {
            this->diag_ethr = std::max(this->diag_ethr, static_cast<Real>(GlobalV::PW_DIAG_THR));
        }
    }
    else
    {
        if (iter == 2)
        {
            this->diag_ethr = 1.e-2;
        }
        this->diag_ethr = std::min(this->diag_ethr,
                                   static_cast<Real>(0.1) * drho
                                       / std::max(static_cast<Real>(1.0), static_cast<Real>(GlobalV::nelec)));
    }
    // It is essential for single precision implementation to keep the diag_ethr value
    // less or equal to the single-precision limit of convergence(0.5e-4).
    // modified by denghuilu at 2023-05-15
    if (GlobalV::precision_flag == "single")
    {
        this->diag_ethr = std::max(this->diag_ethr, static_cast<Real>(0.5e-4));
    }
    return this->diag_ethr;
}

template <typename T, typename Device>
typename HSolverPW<T, Device>::Real HSolverPW<T, Device>::reset_diagethr(std::ofstream& ofs_running,
                                                                         const Real hsover_error,
                                                                         const Real drho)
{
    ofs_running << " Notice: Threshold on eigenvalues was too large.\n";
    ModuleBase::WARNING("scf", "Threshold on eigenvalues was too large.");
    ofs_running << " hsover_error=" << hsover_error << " > DRHO=" << drho << std::endl;
    ofs_running << " Origin diag_ethr = " << this->diag_ethr << std::endl;
    this->diag_ethr = 0.1 * drho / GlobalV::nelec;
    // It is essential for single precision implementation to keep the diag_ethr value
    // less or equal to the single-precision limit of convergence(0.5e-4).
    // modified by denghuilu at 2023-05-15
    if (GlobalV::precision_flag == "single")
    {
        this->diag_ethr = std::max(this->diag_ethr, static_cast<Real>(0.5e-4));
    }
    ofs_running << " New    diag_ethr = " << this->diag_ethr << std::endl;
    return this->diag_ethr;
}

template class HSolverPW<std::complex<float>, base_device::DEVICE_CPU>;
template class HSolverPW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class HSolverPW<std::complex<float>, base_device::DEVICE_GPU>;
template class HSolverPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver
