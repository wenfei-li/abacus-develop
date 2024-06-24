#include "grid_technique.h"

#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hsolver/kernels/cuda/helper_cuda.h"

Grid_Technique::Grid_Technique()
{
    allocate_find_R2 = false;
#if ((defined __CUDA) /* || (defined __ROCM) */)
    if (GlobalV::device_flag == "gpu")
    {
        is_malloced = false;
    }
#endif
}

Grid_Technique::~Grid_Technique()
{

#if ((defined __CUDA) /* || (defined __ROCM) */)
    if (GlobalV::device_flag == "gpu")
    {
        free_gpu_gint_variables(this->nat);
    }
#endif
}

// This function is called in esolver_ks_lcao_elec.cpp
// after the orbital information has been read,
// this function control the routinue to generate
// grid technique parameters.
void Grid_Technique::set_pbc_grid(const int& ncx_in,
                                  const int& ncy_in,
                                  const int& ncz_in,
                                  const int& bx_in,
                                  const int& by_in,
                                  const int& bz_in,
                                  const int& nbx_in,
                                  const int& nby_in,
                                  const int& nbz_in,
                                  const int& nbxx_in,
                                  const int& nbzp_start_in,
                                  const int& nbzp_in,
                                  const int& ny,
                                  const int& nplane,
                                  const int& startz_current,
                                  const UnitCell& ucell,
                                  const double& dr_uniform,
                                  const std::vector<double>& rcuts,
                                  const std::vector<std::vector<double>>& psi_u,
                                  const std::vector<std::vector<double>>& dpsi_u,
                                  const std::vector<std::vector<double>>& d2psi_u,
                                  const int& num_stream)
{
    ModuleBase::TITLE("Grid_Technique", "init");
    ModuleBase::timer::tick("Grid_Technique", "init");

    if (GlobalV::OUT_LEVEL != "m")
    {
        GlobalV::ofs_running << "\n SETUP EXTENDED REAL SPACE GRID FOR GRID INTEGRATION" << std::endl;
    }

    // copy ucell and orb parameters
    this->ucell = &ucell;
    this->dr_uniform = dr_uniform;

    this->nwmax = ucell.nwmax;
    this->ntype = ucell.ntype;

    this->rcuts = rcuts;
    double max_cut = *std::max_element(this->rcuts.begin(), this->rcuts.end());
    this->nr_max = static_cast<int>(1 / this->dr_uniform * max_cut) + 10;
    this->psi_u = psi_u;
    this->dpsi_u = dpsi_u;
    this->d2psi_u = d2psi_u;

    // (1) init_meshcell cell and big cell.
    this->set_grid_dim(ncx_in,
                       ncy_in,
                       ncz_in,
                       bx_in,
                       by_in,
                       bz_in,
                       nbx_in,
                       nby_in,
                       nbz_in,
                       nbxx_in,
                       nbzp_start_in,
                       nbzp_in);
    this->init_latvec(ucell);

    this->init_big_latvec(ucell);

    this->init_meshcell_pos();

    // (2) expand the grid

    this->init_grid_expansion(ucell, this->rcuts.data());

    // (3) calculate the extended grid.
    this->cal_extended_cell(this->dxe, this->dye, this->dze, this->nbx, this->nby, this->nbz);

    this->init_tau_in_bigcell(ucell);

    this->init_meshball();

    this->init_atoms_on_grid(ny, nplane, startz_current, ucell);

    this->cal_trace_lo(ucell);
#if ((defined __CUDA) /* || (defined __ROCM) */)
    if (GlobalV::device_flag == "gpu")
    {
        this->init_gpu_gint_variables(ucell, num_stream);
    }
#endif

    ModuleBase::timer::tick("Grid_Technique", "init");
    return;
}

void Grid_Technique::get_startind(const int& ny, const int& nplane, const int& startz_current)
{
    ModuleBase::TITLE("Grid_Technique", "get_startind");

    assert(nbxx >= 0);

    // calculates start_ind, which stores the
    // starting index of each bigcell
    this->start_ind = std::vector<int>(nbxx, 0);
    ModuleBase::Memory::record("GT::start_ind", sizeof(int) * nbxx);

    for (int i = 0; i < nbxx; i++)
    {
        int ibx = 0;
        int iby = 0;
        int ibz = 0;

        int ix = 0;
        int iy = 0;
        int iz = 0;

        ibx = i / (nby * nbzp);
        iby = (i - ibx * nby * nbzp) / nbzp;
        ibz = i % nbzp;

        ix = ibx * this->bx;
        iy = iby * this->by;
        iz = (ibz + nbzp_start) * this->bz - startz_current;

        int ind = iz + iy * nplane + ix * ny * nplane;

        start_ind[i] = ind;
    }

    return;
}

// PLEASE update this 'init_atoms_on_grid' to make
// it adapted to 'cuboid' shape of grid
// mohan add 2021-04-06
void Grid_Technique::init_atoms_on_grid(const int& ny,
                                        const int& nplane,
                                        const int& startz_current,
                                        const UnitCell& ucell)
{
    ModuleBase::TITLE("Grid_Technique", "init_atoms_on_grid");

    assert(nbxx >= 0);
    this->get_startind(ny, nplane, startz_current);

    // (1) prepare data.
    // counting the number of atoms whose orbitals have
    // values on the bigcell.
    this->how_many_atoms = std::vector<int>(nbxx, 0);
    ModuleBase::Memory::record("GT::how_many_atoms", sizeof(int) * nbxx);

    // (2) information about gloabl grid
    // and local grid.
    // mohan add 2010-07-02
    std::vector<int> ind_bigcell = std::vector<int>(nbxyz, 0);
    ModuleBase::Memory::record("GT::ind_bigcell", sizeof(int) * this->nxyze);
    std::vector<char> bigcell_on_processor = std::vector<char>(nbxyz, 0);
    ModuleBase::Memory::record("GT::bigcell_on_processor", sizeof(char) * this->nxyze);
    this->check_bigcell(ind_bigcell.data(), bigcell_on_processor.data());

    // (3) Find the atoms using
    // when doing grid integration.
    this->in_this_processor = std::vector<bool>(ucell.nat, false);
    ModuleBase::Memory::record("GT::in_this_processor", sizeof(int) * this->nxyze);

    // (4) init atoms on grid
    std::vector<int> index2normal = std::vector<int>(this->nxyze, 0);
    ModuleBase::Memory::record("GT::index2normal", sizeof(int) * this->nxyze);
    this->grid_expansion_index(true, index2normal.data());

    // (5) record how many atoms on
    // each local grid point (ix,iy,iz)
    int iat = 0;
    int normal = 0;
    this->total_atoms_on_grid = 0;
    int nat_local = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            for (int im = 0; im < this->meshball_ncells; im++)
            {
                // bcell[iat]: which bcell iat atom is in.
                // ball[im]: relative position of adjacent bcell.
                normal = index2normal[this->index_atom[iat] + this->index_ball[im]];

                if (normal >= nbxyz)
                {
                    std::cout << " index_atom=" << index_atom[iat] << std::endl;
                    std::cout << " index_ball=" << index_ball[im] << std::endl;
                    std::cout << " normal=" << normal << std::endl;
                    std::cout << " nbxyz=" << nbxyz << std::endl;
                    ModuleBase::WARNING_QUIT("Grid_Technique::init_atoms_on_grid", "normal >= nbxyz");
                }

                assert(normal >= 0);

                int f = ind_bigcell[normal];
                if (!bigcell_on_processor[normal])
                    continue;

                ++how_many_atoms[f];
                ++total_atoms_on_grid;

                this->in_this_processor[iat] = true;
            }
            if (this->in_this_processor[iat])
                ++nat_local;
            ++iat;
        }
    }

    if (GlobalV::test_gridt)
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Total_atoms_on_grid", total_atoms_on_grid);

    int stop = 0;
    if (total_atoms_on_grid == 0)
    {
        GlobalV::ofs_running << " No atoms on this sub-FFT-mesh." << std::endl;
        stop = 1;
    }
    Parallel_Reduce::reduce_all(stop);
    if (stop)
    {
        ModuleBase::WARNING("Grid_Technique::init_atoms_on_grid", "No atom on this sub-FFT-mesh.");
    }

    // calculate the trach of local ia to global iat
    if (nat_local > 0)
    {
        this->trace_iat.resize(nat_local);
        for (int iat = ucell.nat - 1; iat >= 0; iat--)
        {
            if (this->in_this_processor[iat])
            {
                this->trace_iat[--nat_local] = iat;
            }
        }
    }

    // need how_many_atoms first.
    this->cal_grid_integration_index();
    // bcell_start is needed.
    this->init_atoms_on_grid2(index2normal.data(), ucell);
    return;
}

void Grid_Technique::check_bigcell(int* ind_bigcell, char* bigcell_on_processor)
{
    // check if a given bigcell is treated on this processor
    const int zstart = nbzp_start;
    const int zend = nbzp + zstart;
    const int nbyz = nby * nbz;
    const int nz = nbzp;

    int iz_now = 0;
    int ix = 0;
    int iy = 0;
    int iz = 0;
    int ind = 0;
    bool flag = false;

    for (int i = 0; i < nbxyz; i++)
    {
        int iz_now = i % nbz;
        if (iz_now < zstart || iz_now >= zend)
        {
            flag = false;
        }
        else
        {
            flag = true;
            ix = i / nbyz;
            iy = (i - ix * nbyz) / nbz;
            iz = iz_now - zstart;
            ind = ix * nby * nz + iy * nz + iz;
            // no need to calculate index if bigcell is
            // not on this processor
        }

        ind_bigcell[i] = ind;
        bigcell_on_processor[i] = flag;
    }
    return;
}

void Grid_Technique::init_atoms_on_grid2(const int* index2normal, const UnitCell& ucell)
{
    ModuleBase::TITLE("Grid_Techinique", "init_atoms_on_grid2");

    if (total_atoms_on_grid == 0)
    {
        ModuleBase::WARNING("Grid_Technique::init_atoms_on_grid2", "no atom on this sub FFT grid.");
        return;
    }

    std::vector<int> index2ucell = std::vector<int>(this->nxyze, 0);
    ModuleBase::Memory::record("GT::index2ucell", sizeof(int) * this->nxyze);
    this->grid_expansion_index(false, index2ucell.data());

    std::vector<int> ind_bigcell = std::vector<int>(nbxyz, 0);
    ModuleBase::Memory::record("GT::ind_bigcell", sizeof(int) * nbxyz);
    std::vector<char> bigcell_on_processor = std::vector<char>(nbxyz, 0);
    this->check_bigcell(ind_bigcell.data(), bigcell_on_processor.data());

    //--------------------------------------
    // save which atom is in the bigcell,unitcell
    //--------------------------------------
    assert(total_atoms_on_grid != 0);
    this->which_atom = std::vector<int>(total_atoms_on_grid, 0);
    ModuleBase::Memory::record("GT::which_atom", sizeof(int) * total_atoms_on_grid);

    this->which_bigcell = std::vector<int>(total_atoms_on_grid, 0);
    ModuleBase::Memory::record("GT::which_bigcell", sizeof(int) * total_atoms_on_grid);

    this->which_unitcell = std::vector<int>(total_atoms_on_grid, 0);
    ModuleBase::Memory::record("GT::which_unitcell", sizeof(int) * total_atoms_on_grid);

    // for each atom, first we need to locate which cell
    // the atom is in, then we search meshball aroung this
    // grid, and record each grid's atom position.
    int count = 0;
    int iat = 0;
    this->how_many_atoms = std::vector<int>(nbxx, 0);
    ModuleBase::Memory::record("GT::how many atoms", sizeof(int) * nbxx);
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            // zero bigcell of meshball indicate ?
            for (int im = 0; im < this->meshball_ncells; im++)
            {
                const int extgrid = this->index_atom[iat] + this->index_ball[im];
                const int normal = index2normal[extgrid];

                // mohan add 2010-07-01
                int f = ind_bigcell[normal];
                if (!bigcell_on_processor[normal])
                {
                    continue;
                }

                // it's not the normal order to calculate which_atom
                // and which_bigcell, especailly in 1D array.
                // Each grid's adjacent atom number is different,
                // so, first we need to locate which grid, using
                // bcell_start, then we need to count which adjacent atom.
                // using how_many_atoms.
                int index = this->bcell_start[f] + this->how_many_atoms[f];

                // we save which_atom and which_bigcell in 1D array,
                // once you want to use this in grid integration,
                // the only information you got is the 'normal' index,
                // so you need to use bcell_start
                // to get the 'mesh_index', then you can you this mesh_index
                // to use which_atom or which_bigcell.
                this->which_atom[index] = iat;
                this->which_bigcell[index] = im;
                this->which_unitcell[index] = index2ucell[extgrid];

                ++how_many_atoms[f];
                ++count;
            }
            ++iat;
        }
    }
    assert(count == total_atoms_on_grid);
    return;
}

void Grid_Technique::cal_grid_integration_index()
{
    // save the start
    this->bcell_start = std::vector<int>(nbxx, 0);
    ModuleBase::Memory::record("GT::bcell_start", sizeof(int) * nbxx);
    for (int i = 1; i < nbxx; i++)
    {
        this->bcell_start[i] = this->bcell_start[i - 1] + this->how_many_atoms[i - 1];
    }

    // calculate which grid has the largest number of atoms,
    // and how many atoms.
    this->max_atom = 0;
    for (int i = 0; i < nbxx; i++)
    {
        this->max_atom = std::max(this->max_atom, this->how_many_atoms[i]);
    }

#ifdef __MPI
    int* all = new int[GlobalV::NPROC];
    ModuleBase::GlobalFunc::ZEROS(all, GlobalV::NPROC);
    Parallel_Reduce::gather_int_all(max_atom, all);
    if (GlobalV::MY_RANK == 0)
    {
        GlobalV::ofs_warning << std::setw(15) << "Processor" << std::setw(15) << "Atom" << std::endl;
        for (int i = 0; i < GlobalV::NPROC; i++)
        {
            GlobalV::ofs_warning << std::setw(15) << i + 1 << std::setw(15) << all[i] << std::endl;
        }
    }
    delete[] all;
#endif

    if (GlobalV::test_gridt)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Max atom on bigcell", max_atom);
    }
    return;
}

// set 'lgd' variable
void Grid_Technique::cal_trace_lo(const UnitCell& ucell)
{
    ModuleBase::TITLE("Grid_Technique", "cal_trace_lo");
    // save the atom information in trace_lo,
    // in fact the trace_lo dimension can be reduced
    // to ucell.nat, but I think this is another way.
    this->trace_lo = std::vector<int>(GlobalV::NLOCAL, -1);
    ModuleBase::Memory::record("GT::trace_lo", sizeof(int) * GlobalV::NLOCAL);

    this->lnat = 0;
    this->lgd = 0;
    int iat = 0;
    int iw_all = 0;
    int iw_local = 0;

    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            if (this->in_this_processor[iat])
            {
                ++lnat;
                int nw0 = ucell.atoms[it].nw;
                if (GlobalV::NSPIN == 4)
                { // added by zhengdy-soc, need to be double in soc
                    nw0 *= 2;
                    this->lgd += nw0;
                }
                else
                {
                    this->lgd += ucell.atoms[it].nw;
                }

                for (int iw = 0; iw < nw0; iw++)
                {
                    this->trace_lo[iw_all] = iw_local;
                    ++iw_local;
                    ++iw_all;
                }
            }
            else
            {
                // global index of atomic orbitals
                iw_all += ucell.atoms[it].nw;
                if (GlobalV::NSPIN == 4)
                    iw_all += ucell.atoms[it].nw;
            }
            ++iat;
        }
    }

    if (GlobalV::OUT_LEVEL != "m")
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Atom number in sub-FFT-grid", lnat);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Local orbitals number in sub-FFT-grid", lgd);
    }

    assert(iw_local == lgd);
    assert(iw_all == GlobalV::NLOCAL);
    return;
}

#if ((defined __CUDA) /* || (defined __ROCM) */)

void Grid_Technique::init_gpu_gint_variables(const UnitCell& ucell, const int num_stream)
{

    int dev_id = base_device::information::set_device_by_rank();
    if (is_malloced)
    {
        free_gpu_gint_variables(this->nat);
    }
    nstreams = num_stream;
    double ylmcoef[100];
    ModuleBase::GlobalFunc::ZEROS(ylmcoef, 100);
    for (int i = 0; i < 100; i++)
    {
        ylmcoef[i] = ModuleBase::Ylm::ylmcoef[i];
    }
    checkCudaErrors(cudaMalloc((void**)&ylmcoef_g, 100 * sizeof(double)));
    checkCudaErrors(cudaMemcpy(ylmcoef_g, ylmcoef, 100 * sizeof(double), cudaMemcpyHostToDevice));

    double max_cut = *std::max_element(this->rcuts.begin(), this->rcuts.end());

    int atom_nw_now[ucell.ntype];
    int ucell_atom_nwl_now[ucell.ntype];
    for (int i = 0; i < ucell.ntype; i++)
    {
        atom_nw_now[i] = ucell.atoms[i].nw;
        ucell_atom_nwl_now[i] = ucell.atoms[i].nwl;
    }

    nr_max = static_cast<int>(1000 * max_cut) + 10;
    // double psi_u_now[ucell.ntype * ucell.nwmax * nr_max *
    // 2];
    double* psi_u_now = (double*)malloc(ucell.ntype * ucell.nwmax * nr_max * 2 * sizeof(double));
    memset(psi_u_now, 0, ucell.ntype * ucell.nwmax * nr_max * 2 * sizeof(double));
    bool* atom_iw2_new_now = (bool*)malloc(ucell.ntype * ucell.nwmax * sizeof(bool));
    memset(atom_iw2_new_now, 0, ucell.ntype * ucell.nwmax * sizeof(bool));
    int* atom_iw2_ylm_now = (int*)malloc(ucell.ntype * ucell.nwmax * sizeof(int));
    memset(atom_iw2_ylm_now, 0, ucell.ntype * ucell.nwmax * sizeof(int));
    int* atom_iw2_l_now = (int*)malloc(ucell.ntype * ucell.nwmax * sizeof(int));
    memset(atom_iw2_l_now, 0, ucell.ntype * ucell.nwmax * sizeof(int));

    Atom* atomx;
    for (int i = 0; i < ucell.ntype; i++)
    {
        atomx = &ucell.atoms[i];
        for (int j = 0; j < ucell.nwmax; j++)
        {
            if (j < atomx->nw)
            {
                atom_iw2_new_now[i * ucell.nwmax + j] = atomx->iw2_new[j];
                atom_iw2_ylm_now[i * ucell.nwmax + j] = atomx->iw2_ylm[j];
                atom_iw2_l_now[i * ucell.nwmax + j] = atomx->iw2l[j];
                for (int k = 0; k < this->nr_max; k++)
                {
                    int index_temp = (i * ucell.nwmax * this->nr_max + j * this->nr_max + k) * 2;
                    if (k < this->psi_u[i * this->nwmax + j].size())
                    {
                        psi_u_now[index_temp] = this->psi_u[i * this->nwmax + j].data()[k];
                        psi_u_now[index_temp + 1] = this->dpsi_u[i * this->nwmax + j].data()[k];
                    }
                }
            }
        }
    }

    checkCudaErrors(cudaMalloc((void**)&atom_nw_g, ucell.ntype * sizeof(int)));
    checkCudaErrors(cudaMemcpy(atom_nw_g, atom_nw_now, ucell.ntype * sizeof(int), cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMalloc((void**)&atom_nwl_g, ucell.ntype * sizeof(int)));
    checkCudaErrors(cudaMemcpy(atom_nwl_g, ucell_atom_nwl_now, ucell.ntype * sizeof(int), cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMalloc((void**)&psi_u_g, ucell.ntype * ucell.nwmax * nr_max * sizeof(double) * 2));
    checkCudaErrors(cudaMemcpy(psi_u_g,
                               psi_u_now,
                               ucell.ntype * ucell.nwmax * nr_max * sizeof(double) * 2,
                               cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMalloc((void**)&atom_new_g, ucell.ntype * ucell.nwmax * sizeof(bool)));
    checkCudaErrors(cudaMalloc((void**)&atom_ylm_g, ucell.ntype * ucell.nwmax * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&atom_l_g, ucell.ntype * ucell.nwmax * sizeof(int)));

    checkCudaErrors(
        cudaMemcpy(atom_new_g, atom_iw2_new_now, ucell.ntype * ucell.nwmax * sizeof(bool), cudaMemcpyHostToDevice));
    checkCudaErrors(
        cudaMemcpy(atom_ylm_g, atom_iw2_ylm_now, ucell.ntype * ucell.nwmax * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(
        cudaMemcpy(atom_l_g, atom_iw2_l_now, ucell.ntype * ucell.nwmax * sizeof(int), cudaMemcpyHostToDevice));

    gemm_algo_selector(bxyz, fastest_matrix_mul, ucell);

    is_malloced = true;

    free(psi_u_now);
    free(atom_iw2_new_now);
    free(atom_iw2_ylm_now);
}

void Grid_Technique::free_gpu_gint_variables(int nat)
{
    if (!is_malloced)
    {
        return;
    }

    checkCudaErrors(cudaFree(ylmcoef_g));
    checkCudaErrors(cudaFree(atom_nwl_g));
    checkCudaErrors(cudaFree(psi_u_g));
    checkCudaErrors(cudaFree(atom_new_g));
    checkCudaErrors(cudaFree(atom_ylm_g));
    checkCudaErrors(cudaFree(atom_nw_g));
    checkCudaErrors(cudaFree(atom_l_g));

    is_malloced = false;
}
#endif
