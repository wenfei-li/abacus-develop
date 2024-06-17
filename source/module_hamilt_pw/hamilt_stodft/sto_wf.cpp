#include "sto_wf.h"

#include <cassert>

#include "module_base/memory.h"
#include "time.h"

//---------Temporary------------------------------------
#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//------------------------------------------------------

Stochastic_WF::Stochastic_WF()
{
}

Stochastic_WF::~Stochastic_WF()
{
    delete chi0;
    delete shchi;
    delete chiortho;
    delete[] nchip;
}

void Stochastic_WF::init(K_Vectors* p_kv, const int npwx_in)
{
    this->nks = p_kv->get_nks();
    this->ngk = p_kv->ngk.data();
    this->npwx = npwx_in;
    nchip = new int[nks];

    if (nks <= 0)
    {
        ModuleBase::WARNING_QUIT("Stochastic_WF", "nks <=0!");
    }
}

void Init_Sto_Orbitals(Stochastic_WF& stowf, const int seed_in)
{
    if (seed_in == 0 || seed_in == -1)
    {
        srand((unsigned)time(nullptr) + GlobalV::MY_RANK * 10000); // GlobalV global variables are reserved
    }
    else
    {
        srand((unsigned)std::abs(seed_in) + GlobalV::MY_RANK * 10000);
    }

    Allocate_Chi0(stowf);
    Update_Sto_Orbitals(stowf, seed_in);
}

void Allocate_Chi0(Stochastic_WF& stowf)
{
    bool firstrankmore = false;
    int igroup = 0;
    // I am not sure which is better.
    // former processor calculate more bands
    if (firstrankmore)
    {
        igroup = GlobalV::MY_STOGROUP;
    }
    // latter processor calculate more bands
    else
    {
        igroup = GlobalV::NSTOGROUP - GlobalV::MY_STOGROUP - 1;
    }
    const int nchi = INPUT.nbands_sto;
    const int npwx = stowf.npwx;
    const int nks = stowf.nks;
    const int ngroup = GlobalV::NSTOGROUP;
    if (ngroup <= 0)
    {
        ModuleBase::WARNING_QUIT("Init_Sto_Orbitals", "ngroup <= 0!");
    }
    int tmpnchip = int(nchi / ngroup);
    if (igroup < nchi % ngroup)
    {
        ++tmpnchip;
    }

    stowf.nchip_max = tmpnchip;
    size_t size = stowf.nchip_max * npwx * nks;
    stowf.chi0 = new psi::Psi<std::complex<double>>(nks, stowf.nchip_max, npwx, stowf.ngk);
    ModuleBase::Memory::record("SDFT::chi0", size * sizeof(std::complex<double>));

    for (int ik = 0; ik < nks; ++ik)
    {
        stowf.nchip[ik] = tmpnchip;
    }
}

void Update_Sto_Orbitals(Stochastic_WF& stowf, const int seed_in)
{
    const int nchi = INPUT.nbands_sto;
    stowf.chi0->fix_k(0);
    if (seed_in >= 0)
    {
        for (int i = 0; i < stowf.chi0->size(); ++i)
        {
            const double phi = 2 * ModuleBase::PI * rand() / double(RAND_MAX);
            stowf.chi0->get_pointer()[i] = std::complex<double>(cos(phi), sin(phi)) / sqrt(double(nchi));
        }
    }
    else
    {
        for (int i = 0; i < stowf.chi0->size(); ++i)
        {
            if (rand() / double(RAND_MAX) < 0.5)
            {
                stowf.chi0->get_pointer()[i] = -1.0 / sqrt(double(nchi));
            }
            else
            {
                stowf.chi0->get_pointer()[i] = 1.0 / sqrt(double(nchi));
            }
        }
    }
}

#ifdef __MPI
void Init_Com_Orbitals(Stochastic_WF& stowf)
{
    const bool firstrankmore = false;
    const int npwx = stowf.npwx;
    const int nks = stowf.nks;
    int igroup;
    // former processor calculate more bands
    if (firstrankmore)
    {
        igroup = GlobalV::MY_STOGROUP;
    }
    // latter processor calculate more bands
    else
    {
        igroup = GlobalV::NSTOGROUP - GlobalV::MY_STOGROUP - 1;
    }
    const int ngroup = GlobalV::NSTOGROUP;
    const int n_in_pool = GlobalV::NPROC_IN_POOL;
    const int i_in_group = GlobalV::RANK_IN_STOGROUP;
    const int i_in_pool = GlobalV::RANK_IN_POOL;

    int* totnpw = new int[nks];
    for (int ik = 0; ik < nks; ++ik)
    {
        int* npwip = new int[n_in_pool];
        const int npw = stowf.ngk[ik];
        totnpw[ik] = 0;

        MPI_Allgather(&npw, 1, MPI_INT, npwip, 1, MPI_INT, POOL_WORLD);
        for (int i_in_p = 0; i_in_p < n_in_pool; ++i_in_p)
        {
            totnpw[ik] += npwip[i_in_p];
        }

        int tmpnchip = int(totnpw[ik] / ngroup);
        if (igroup < totnpw[ik] % ngroup)
        {
            ++tmpnchip;
        }
        stowf.nchip[ik] = tmpnchip;
        stowf.nchip_max = std::max(tmpnchip, stowf.nchip_max);
        delete[] npwip;
    }
    size_t size = stowf.nchip_max * npwx * nks;
    stowf.chi0 = new psi::Psi<std::complex<double>>(nks, stowf.nchip_max, npwx, stowf.ngk);
    stowf.chi0->zero_out();
    ModuleBase::Memory::record("SDFT::chi0", size * sizeof(std::complex<double>));
    for (int ik = 0; ik < nks; ++ik)
    {
        int* npwip = new int[n_in_pool];
        const int npw = stowf.ngk[ik];
        MPI_Allgather(&npw, 1, MPI_INT, npwip, 1, MPI_INT, POOL_WORLD);
        const int re = totnpw[ik] % ngroup;
        int ip = 0, ig0 = 0;
        const int nchipk = stowf.nchip[ik];
        // give value to orbitals in one parallel group one by one.
        for (int ichi = 0; ichi < nchipk; ++ichi)
        {
            int ig;
            if (igroup < re)
            {
                // It has more nchip.
                ig = igroup * nchipk + ichi - ig0;
            }
            else
            {
                // It has less nchip and should add re.
                ig = igroup * nchipk + re + ichi - ig0;
            }
            // get which ip stores this ig.
            while (ig >= npwip[ip])
            {
                ig -= npwip[ip];
                ig0 += npwip[ip];
                ++ip;
            }
            if (i_in_pool == ip)
            {
                stowf.chi0->operator()(ik, ichi, ig) = 1;
            }
        }

        delete[] npwip;
    }
    delete[] totnpw;
}
#else
void Init_Com_Orbitals(Stochastic_WF& stowf)
{
    const int npwx = stowf.npwx;
    const int nks = stowf.nks;
    size_t size = stowf.nchip_max * npwx * nks;
    stowf.chi0 = new psi::Psi<std::complex<double>>(nks, npwx, npwx, stowf.ngk);
    stowf.chi0->zero_out();
    ModuleBase::Memory::record("SDFT::chi0", size * sizeof(std::complex<double>));
    for (int ik = 0; ik < nks; ++ik)
    {
        const int npw = stowf.ngk[ik];
        stowf.nchip[ik] = npwx;
        stowf.nchip_max = npwx;
        for (int ichi = 0; ichi < npw; ++ichi)
        {
            stowf.chi0->operator()(ik, ichi, ichi) = 1;
        }
    }
}
#endif
void Init_Sto_Orbitals_Ecut(Stochastic_WF& stowf,
                            const int seed_in,
                            const K_Vectors& kv,
                            const ModulePW::PW_Basis_K& wfcpw,
                            const int max_ecut)
{
    Allocate_Chi0(stowf);

    ModulePW::PW_Basis pwmax;
#ifdef __MPI
    pwmax.initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
    pwmax.initgrids(wfcpw.lat0, wfcpw.latvec, max_ecut);
    const int nx = pwmax.nx;
    const int ny = pwmax.ny;
    const int nz = pwmax.nz;
    const int nkstot = kv.get_nkstot();
    const int nks = kv.get_nks();
    const int nchitot = INPUT.nbands_sto;
    bool* updown = new bool[nx * ny * nz];
    int* nrecv = new int[GlobalV::NSTOGROUP];
    const int nchiper = stowf.nchip[0];
#ifdef __MPI
    MPI_Allgather(&nchiper, 1, MPI_INT, nrecv, 1, MPI_INT, PARAPW_WORLD);
#endif
    int ichi_start = 0;
    for (int i = 0; i < GlobalV::MY_STOGROUP; ++i)
    {
        ichi_start += nrecv[i];
    }

    for (int ik = 0; ik < nks; ++ik)
    {
        const int iktot = kv.getik_global(ik);
        const int npw = wfcpw.npwk[ik];
        int* ig2ixyz = new int[npw];

        for (int ig = 0; ig < npw; ++ig)
        {
            ModuleBase::Vector3<double> gdirect = wfcpw.getgdirect(ik, ig);
            int ix = static_cast<int>(gdirect.x);
            int iy = static_cast<int>(gdirect.y);
            int iz = static_cast<int>(gdirect.z);
            ix = (ix + nx) % nx;
            iy = (iy + ny) % ny;
            iz = (iz + nz) % nz;
            ig2ixyz[ig] = ix * ny * nz + iy * nz + iz;
        }

        for (int ichi = 0; ichi < nchiper; ++ichi)
        {
            unsigned int seed = std::abs(seed_in) * (nkstot * nchitot) + iktot * nchitot + ichi_start + ichi;
            srand(seed);
            for (int i = 0; i < nx * ny * nz; ++i)
            {
                updown[i] = (rand() / double(RAND_MAX) < 0.5);
            }

            for (int ig = 0; ig < npw; ++ig)
            {
                if (updown[ig2ixyz[ig]])
                {
                    stowf.chi0->operator()(ik, ichi, ig) = -1.0 / sqrt(double(nchitot));
                }
                else
                {
                    stowf.chi0->operator()(ik, ichi, ig) = 1.0 / sqrt(double(nchitot));
                }
            }
        }
        delete[] ig2ixyz;
    }
    delete[] nrecv;
    delete[] updown;
}
