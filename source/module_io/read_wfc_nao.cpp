#include "read_wfc_nao.h"

#include "module_base/parallel_common.h"
#include "module_base/timer.h"
#include "module_io/write_wfc_nao.h"

/**
 * @brief Extracts a submatrix from a global orbital coefficient matrix and stores it in a linear array.
 *        This function is designed for real-valued matrices.
 *
 * @param myid The rank of the current MPI process.
 * @param naroc Array with two elements: number of rows and columns this process is responsible for.
 * @param nb Block size used in the global index calculation.
 * @param dim0 Global number of rows in the matrix.
 * @param dim1 Global number of columns in the matrix.
 * @param iprow The row index in the process grid.
 * @param ipcol The column index in the process grid.
 * @param nbands Total number of bands (relevant matrix columns) involved in the computation.
 * @param work Output array where the submatrix will be stored linearly.
 * @param CTOT Global matrix from which the submatrix is extracted.
 * @return int Always returns 0 as a success indicator.
 */
inline int CTOT2q(const int myid, const int naroc[2], const int nb, const int dim0, const int dim1, const int iprow,
                  const int ipcol, const int nbands, double* work, double** const CTOT)
{
    for (int j = 0; j < naroc[1]; ++j)
    {
        int igcol = Local_Orbital_wfc::globalIndex(j, nb, dim1, ipcol);
        if (igcol >= nbands)
        {
            continue;
        }
        for (int i = 0; i < naroc[0]; ++i)
        {
            int igrow = Local_Orbital_wfc::globalIndex(i, nb, dim0, iprow);
            // GlobalV::ofs_running << "i,j,igcol,igrow" << i << " " << j << " " << igcol << " " << igrow << std::endl;
            if (myid == 0)
            {
                work[j * naroc[0] + i] = CTOT[igcol][igrow];
            }
        }
    }
    return 0;
}

/**
 * @brief Extracts a submatrix from a global orbital coefficient matrix and stores it in a linear array.
 *        This function is designed for complex-valued matrices.
 *
 * @param myid The rank of the current MPI process.
 * @param naroc Array with two elements: number of rows and columns this process is responsible for.
 * @param nb Block size used in the global index calculation.
 * @param dim0 Global number of rows in the matrix.
 * @param dim1 Global number of columns in the matrix.
 * @param iprow The row index in the process grid.
 * @param ipcol The column index in the process grid.
 * @param nbands Total number of bands (relevant matrix columns) involved in the computation.
 * @param work Output array where the submatrix will be stored linearly.
 * @param CTOT Global matrix from which the submatrix is extracted.
 * @return int Always returns 0 as a success indicator.
 */
inline int CTOT2q_c(const int myid, const int naroc[2], const int nb, const int dim0, const int dim1, const int iprow,
                    const int ipcol, const int nbands, std::complex<double>* work, std::complex<double>** const CTOT)
{
    for (int j = 0; j < naroc[1]; ++j)
    {
        int igcol = Local_Orbital_wfc::globalIndex(j, nb, dim1, ipcol);
        if (igcol >= nbands)
        {
            continue;
        }
        for (int i = 0; i < naroc[0]; ++i)
        {
            int igrow = Local_Orbital_wfc::globalIndex(i, nb, dim0, iprow);
            // ofs_running << "i,j,igcol,igrow" << i << " " << j << " " << igcol << " " << igrow << std::endl;
            if (myid == 0)
            {
                work[j * naroc[0] + i] = CTOT[igcol][igrow];
            }
        }
    }
    return 0;
}

// be called in local_orbital_wfc::allocate_k
int ModuleIO::read_wfc_nao_complex(std::complex<double>** ctot, const int& ik, const int& nb2d, const int& nbands_g,
                                   const int& nlocal_g, const std::string& global_readin_dir,
                                   const ModuleBase::Vector3<double> kvec_c, const Parallel_Orbitals* ParaV,
                                   psi::Psi<std::complex<double>>* psi, elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("ModuleIO", "read_wfc_nao_complex");
    ModuleBase::timer::tick("ModuleIO", "read_wfc_nao_complex");

    std::string ss = global_readin_dir + ModuleIO::wfc_nao_gen_fname(1, false, true, ik);
    std::ifstream ifs;
    int error = 0;

    if (GlobalV::DRANK == 0)
    {
        ifs.open(ss.c_str());
        if (!ifs)
        {
            GlobalV::ofs_warning << " Can't open file:" << ss << std::endl;
            error = 1;
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_int(error);
#endif

    if (error == 1)
    {
        return 1;
    }

    // otherwise, find the file.

    if (GlobalV::MY_RANK == 0)
    {
        int ikr = 0;
        double kx = 0.0, ky = 0.0, kz = 0.0;
        int nbands = 0, nlocal = 0;
        ModuleBase::GlobalFunc::READ_VALUE(ifs, ikr);
        ifs >> kx >> ky >> kz;
        ModuleBase::GlobalFunc::READ_VALUE(ifs, nbands);
        ModuleBase::GlobalFunc::READ_VALUE(ifs, nlocal);

        if (ikr != ik + 1)
        {
            GlobalV::ofs_warning << " ikr=" << ikr << " ik=" << ik << std::endl;
            GlobalV::ofs_warning << " k index is not correct" << std::endl;
            error = 4;
        }
        else if (std::abs(kx - kvec_c.x) > 1.0e-5 || std::abs(ky - kvec_c.y) > 1.0e-5
                 || std::abs(kz - kvec_c.z) > 1.0e-5)
        {
            GlobalV::ofs_warning << " k std::vector is not correct" << std::endl;
            GlobalV::ofs_warning << " Read in kx=" << kx << " ky=" << ky << " kz=" << kz << std::endl;
            GlobalV::ofs_warning << " In fact, kx=" << kvec_c.x << " ky=" << kvec_c.y << " kz=" << kvec_c.z
                                 << std::endl;
            error = 4;
        }
        else if (nbands != nbands_g)
        {
            GlobalV::ofs_warning << " read in nbands=" << nbands;
            GlobalV::ofs_warning << " NBANDS=" << nbands_g << std::endl;
            error = 2;
        }
        else if (nlocal != nlocal_g)
        {
            GlobalV::ofs_warning << " read in nlocal=" << nlocal;
            GlobalV::ofs_warning << " NLOCAL=" << nlocal_g << std::endl;
            error = 3;
        }

        ctot = new std::complex<double>*[nbands_g];
        for (int i = 0; i < nbands_g; i++)
        {
            ctot[i] = new std::complex<double>[nlocal_g];
        }

        for (int i = 0; i < nbands_g; ++i)
        {
            int ib = 0;
            ModuleBase::GlobalFunc::READ_VALUE(ifs, ib);
            ib -= 1; // because in C++, ib should start from 0
            //------------------------------------------------
            // read the eigenvalues!
            // very important to determine the occupations.
            //------------------------------------------------
            ModuleBase::GlobalFunc::READ_VALUE(ifs, pelec->ekb(ik, ib));
            ModuleBase::GlobalFunc::READ_VALUE(ifs, pelec->wg(ik, ib));
            assert(i == ib);
            double a = 0.0, b = 0.0;
            for (int j = 0; j < nlocal_g; ++j)
            {
                ifs >> a >> b;
                ctot[i][j] = std::complex<double>(a, b);
                // std::cout << ctot[i][j] << " " << std::endl;
            }
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_int(error);
    // then broadcast the "weigths" of present k-point, with length pelec->wg.nc
    Parallel_Common::bcast_double(&pelec->wg.c[ik * pelec->wg.nc], pelec->wg.nc);
#endif

    if (error == 2)
    {
        return 2;
    }
    if (error == 3)
    {
        return 3;
    }
    if (error == 4)
    {
        return 4;
    }

    ModuleIO::distri_wfc_nao_complex(ctot, ik, nb2d, nbands_g, ParaV, psi);

    // mohan add 2012-02-15,
    // still have bugs, but can solve it later.
    // distribute the wave functions again.

    if (GlobalV::DRANK == 0)
    {
        // delte the ctot
        for (int i = 0; i < nbands_g; i++)
        {
            delete[] ctot[i];
        }
        delete[] ctot;
    }

    //---------
    // TEST
    //---------
    /*
    for(int i=0; i<nbands_g; ++i)
    {
        std::cout << " c band i=" << i+1 << std::endl;
        for(int j=0; j<nlocal_g; ++j)
        {
            std::cout << " " << c[i][j];
        }
        std::cout << std::endl;
    }
    */

    ModuleBase::timer::tick("ModuleIO", "read_wfc_nao_complex");
    return 0;
}

int ModuleIO::read_wfc_nao(double** ctot, const int& is, const bool& gamma_only_local, const int& nb2d,
                           const int& nbands_g, const int& nlocal_g, const std::string& global_readin_dir,
                           const Parallel_Orbitals* ParaV, psi::Psi<double>* psid, elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("ModuleIO", "read_wfc_nao");
    ModuleBase::timer::tick("ModuleIO", "read_wfc_nao");

    std::string ss = global_readin_dir + ModuleIO::wfc_nao_gen_fname(1, GlobalV::GAMMA_ONLY_LOCAL, false, is);
    std::ifstream ifs;

    int error = 0;

    if (GlobalV::DRANK == 0)
    {
        ifs.open(ss.c_str());
        if (!ifs)
        {
            GlobalV::ofs_warning << " Can't open file:" << ss << std::endl;
            error = 1;
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_int(error);
#endif

    if (error == 1)
    {
        return 1;
    }

    // otherwise, find the file.

    if (GlobalV::MY_RANK == 0)
    {
        int nbands = 0, nlocal = 0;
        ModuleBase::GlobalFunc::READ_VALUE(ifs, nbands);
        ModuleBase::GlobalFunc::READ_VALUE(ifs, nlocal);

        if (nbands != nbands_g)
        {
            GlobalV::ofs_warning << " read in nbands=" << nbands;
            GlobalV::ofs_warning << " NBANDS=" << nbands_g << std::endl;
            error = 2;
            ModuleBase::WARNING_QUIT("ModuleIO::read_wfc_nao",
                                     "The nbands read in from file do not match the global parameter NBANDS!");
        }
        else if (nlocal != nlocal_g)
        {
            GlobalV::ofs_warning << " read in nlocal=" << nlocal;
            GlobalV::ofs_warning << " NLOCAL=" << nlocal_g << std::endl;
            error = 3;
            ModuleBase::WARNING_QUIT("ModuleIO::read_wfc_nao",
                                     "The nlocal read in from file do not match the global parameter NLOCAL!");
        }

        ctot = new double*[nbands_g];
        for (int i = 0; i < nbands_g; i++)
        {
            ctot[i] = new double[nlocal_g];
        }

        for (int i = 0; i < nbands_g; i++)
        {
            int ib = 0;
            ModuleBase::GlobalFunc::READ_VALUE(ifs, ib);
            ModuleBase::GlobalFunc::READ_VALUE(ifs, pelec->ekb(is, i));
            ModuleBase::GlobalFunc::READ_VALUE(ifs, pelec->wg(is, i));
            assert((i + 1) == ib);
            // std::cout << " ib=" << ib << std::endl;
            for (int j = 0; j < nlocal_g; j++)
            {
                ifs >> ctot[i][j];
                // std::cout << ctot[i][j] << " ";
            }
            // std::cout << std::endl;
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_int(error);
    Parallel_Common::bcast_double(&(pelec->ekb(is, 0)), nbands_g);
    Parallel_Common::bcast_double(&(pelec->wg(is, 0)), nbands_g);
#endif
    if (error == 2)
    {
        return 2;
    }
    if (error == 3)
    {
        return 3;
    }

    ModuleIO::distri_wfc_nao(ctot, is, nb2d, nbands_g, nlocal_g, ParaV, psid);

    // mohan add 2012-02-15,
    // still have bugs, but can solve it later.
    // distribute the wave functions again.
    // mohan comment out 2021-02-09

    if (GlobalV::MY_RANK == 0)
    {
        // delte the ctot
        for (int i = 0; i < nbands_g; i++)
        {
            delete[] ctot[i];
        }
        delete[] ctot;
    }

    ModuleBase::timer::tick("ModuleIO", "read_wfc_nao");
    return 0;
}

void ModuleIO::distri_wfc_nao(double** ctot, const int& is, const int& nb2d, const int& nbands_g, const int& nlocal_g,
                              const Parallel_Orbitals* ParaV, psi::Psi<double>* psid)
{
    ModuleBase::TITLE("ModuleIO", "distri_wfc_nao");
#ifdef __MPI

    // 1. alloc work array; set some parameters

    long maxnloc = 0; // maximum number of elements in local matrix
    MPI_Reduce(&ParaV->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, ParaV->comm_2D);
    MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, ParaV->comm_2D);
    // reduce and bcast could be replaced by allreduce

    int nprocs = 0, myid = 0;
    MPI_Comm_size(ParaV->comm_2D, &nprocs);
    MPI_Comm_rank(ParaV->comm_2D, &myid);

    double* work = new double[maxnloc]; // work/buffer matrix
#ifdef __DEBUG
    assert(nb2d > 0);
#endif
    const int nb = nb2d;
    int info = 0;
    int naroc[2]; // maximum number of row or column

    // 2. copy from ctot to psi
    for (int iprow = 0; iprow < ParaV->dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < ParaV->dim1; ++ipcol)
        {
            // 2.1 get and bcast local 2d matrix info
            const int coord[2] = {iprow, ipcol};
            int src_rank = 0;
            MPI_Cart_rank(ParaV->comm_2D, coord, &src_rank);
            if (myid == src_rank)
            {
                naroc[0] = ParaV->nrow;
                naroc[1] = ParaV->ncol;
            }
            info = MPI_Bcast(naroc, 2, MPI_INT, src_rank, ParaV->comm_2D);

            // 2.2 copy from nctot to work, then bcast work
            info = CTOT2q(myid, naroc, nb, ParaV->dim0, ParaV->dim1, iprow, ipcol, nbands_g, work, ctot);
            info = MPI_Bcast(work, maxnloc, MPI_DOUBLE, 0, ParaV->comm_2D);
            // GlobalV::ofs_running << "iprow, ipcow : " << iprow << ipcol << std::endl;
            // for (int i=0; i<maxnloc; ++i)
            //{
            // GlobalV::ofs_running << *(work+i)<<" ";
            //}
            // GlobalV::ofs_running << std::endl;
            // 2.3 copy from work to psi
            const int inc = 1;
            if (myid == src_rank)
            {
                BlasConnector::copy(ParaV->nloc_wfc, work, inc, &(psid[0](is, 0, 0)), inc);
            }
        } // loop ipcol
    }     // loop	iprow

    delete[] work;
#else
    for (int i = 0; i < nbands_g; i++)
    {
        for (int j = 0; j < nlocal_g; j++)
        {
            psid[0](is, i, j) = ctot[i][j];
        }
    }
#endif
    return;
}

void ModuleIO::distri_wfc_nao_complex(std::complex<double>** ctot, const int& ik, const int& nb2d, const int& nbands_g,
                                      const Parallel_Orbitals* ParaV, psi::Psi<std::complex<double>>* psi)
{
    ModuleBase::TITLE("ModuleIO", "distri_wfc_nao_complex");
#ifdef __MPI

    // 1. alloc work array; set some parameters

    long maxnloc = 0; // maximum number of elements in local matrix
    MPI_Reduce(&ParaV->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, ParaV->comm_2D);
    MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, ParaV->comm_2D);
    // reduce and bcast could be replaced by allreduce

    int nprocs = 0, myid = 0;
    MPI_Comm_size(ParaV->comm_2D, &nprocs);
    MPI_Comm_rank(ParaV->comm_2D, &myid);

    std::complex<double>* work = new std::complex<double>[maxnloc]; // work/buffer matrix
#ifdef __DEBUG
    assert(nb2d > 0);
#endif
    const int nb = nb2d;
    int info = 0;
    int naroc[2]; // maximum number of row or column

    // 2. copy from ctot to psi
    for (int iprow = 0; iprow < ParaV->dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < ParaV->dim1; ++ipcol)
        {
            // 2.1 get and bcast local 2d matrix info
            const int coord[2] = {iprow, ipcol};
            int src_rank = 0;
            MPI_Cart_rank(ParaV->comm_2D, coord, &src_rank);
            if (myid == src_rank)
            {
                naroc[0] = ParaV->nrow;
                naroc[1] = ParaV->ncol_bands;
            }
            info = MPI_Bcast(naroc, 2, MPI_INT, src_rank, ParaV->comm_2D);

            // 2.2 copy from ctot to work, then bcast work
            info = CTOT2q_c(myid, naroc, nb, ParaV->dim0, ParaV->dim1, iprow, ipcol, nbands_g, work, ctot);
            info = MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, 0, ParaV->comm_2D);
            // ofs_running << "iprow, ipcow : " << iprow << ipcol << std::endl;
            // for (int i=0; i<maxnloc; ++i)
            //{
            // ofs_running << *(work+i)<<" ";
            //}
            // ofs_running << std::endl;
            // 2.3 copy from work to psi_k
            const int inc = 1;
            if (myid == src_rank)
            {
                BlasConnector::copy(ParaV->nloc_wfc, work, inc, &(psi[0](ik, 0, 0)), inc);
            }
        } // loop ipcol
    }     // loop	iprow

    delete[] work;
#else
    ModuleBase::WARNING_QUIT("ModuleIO::distri_wfc_nao", "check the code without MPI.");
#endif
    return;
}
