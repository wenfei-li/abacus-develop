//wenfei 2022-1-11
//This file contains subroutines that contains subroutines for printing density matrices

//There are 2 subroutines for printing density matrices:
//1. print_dm : for gamma only
//2. print_dm_k : for multi-k

#include "LCAO_dftu_new.h"

void LCAO_DftU_New::print_dm(const ModuleBase::matrix &dm)
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

void LCAO_DftU_New::print_dm_k(const int nks, const std::vector<ModuleBase::ComplexMatrix>& dm)
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