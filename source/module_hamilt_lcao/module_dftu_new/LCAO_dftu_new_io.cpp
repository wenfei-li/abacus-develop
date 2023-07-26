//wenfei 2022-1-11
//This file contains subroutines that contains interface with libnpy
//since many arrays must be saved in numpy format
//It also contains subroutines for printing density matrices
//which is used in unit tests

//There are 2 subroutines for printing density matrices:
//1. print_dm : for gamma only
//2. print_dm_k : for multi-k

//There are 4 subroutines in this file that prints to npy file:
//1. save_npy_d : descriptor ->dm_eig.npy
//2. save_npy_gvx : gvx ->grad_vx.npy
//3. save_npy_e : energy
//4. save_npy_f : force
//5. save_npy_s : stress
//6. save_npy_o : bandgap
//7. save_npy_orbital_precalc : orbital_precalc

#include "LCAO_dftu_new.h"

void LCAO_Deepks::print_dm(const ModuleBase::matrix &dm)
{
    std::ofstream ofs("dm");
    ofs << std::setprecision(15);
    for (int mu=0;mu<GlobalV::NLOCAL;mu++)
    {
        for (int nu=0;nu<GlobalV::NLOCAL;nu++)
        {
            ofs << dm(mu,nu) << " ";
        }
        ofs << std::endl;
    }
}

void LCAO_Deepks::print_dm_k(const int nks, const std::vector<ModuleBase::ComplexMatrix>& dm)
{
    std::stringstream ss;
    for(int ik=0;ik<nks;ik++)
    {
        ss.str("");
        ss<<"dm_"<<ik;
        std::ofstream ofs(ss.str().c_str());
        ofs << std::setprecision(15);

        for (int mu=0;mu<GlobalV::NLOCAL;mu++)
        {
            for (int nu=0;nu<GlobalV::NLOCAL;nu++)
            {
                ofs << dm[ik](mu,nu) << " ";
            }
            ofs << std::endl;
        }
    }
}