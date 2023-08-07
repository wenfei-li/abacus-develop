#include "LCAO_dftu_new.h"

void LCAO_DftU_New::cal_gedm(const int nat)
{
    ModuleBase::TITLE("LCAO_DftU_New", "cal_gedm");
    ModuleBase::timer::tick ("LCAO_DftU_New","cal_gedm");

    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    for(int is = 0; is < GlobalV::NSPIN; is ++)
    {
        for (int inl = 0;inl < this->inlmax;inl++)
        {
            ModuleBase::GlobalFunc::ZEROS(this->gedm[is][inl], pdm_size);
        }
    }

    this->E_delta = 0.0;

    for(int iat = 0; iat < nat; iat ++)
    {
        if(!if_has_u[iat]) continue;
        for(int is = 0; is < GlobalV::NSPIN; is ++)
        {
            for(int inl = 0; inl < this->inlmax; inl ++)
            {
                const int l = this->inl_l[inl];
                if(std::abs(uvalue[iat][l]) < 1e-8) continue;

                double noccup = 0.0;
                for(int m1 = 0; m1 < 2*l+1; m1 ++)
                {
                    int ind = m1 * (2*l+1) + m1;
                    noccup += this->pdm[is][inl][ind];
                }

                for(int m1 = 0; m1 < 2*l+1; m1 ++)
                {
                    for(int m2 = 0; m2 < 2*l+1; m2 ++)
                    {
                        int ind = m1 * (2*l+1) + m2;
                        if(m1 == m2)
                        {
                            this->gedm[is][inl][ind] = uvalue[iat][l] * (0.5 - noccup);
                        }
                        else
                        {
                            this->gedm[is][inl][ind] = uvalue[iat][l] * (- noccup);
                        }
                    }
                }

                // Tr (nn) = \sum_mm' n_{mm'} n_{m'm}
                double trace = 0.0;
                for(int m1 = 0; m1 < 2*l+1; m1 ++)
                {
                    for(int m2 = 0; m2 < 2*l+1; m2 ++)
                    {
                        int ind1 = m1 * (2*l+1) + m2;
                        int ind2 = m2 * (2*l+1) + m1;
                        trace += this->pdm[is][inl][ind1] * this->pdm[is][inl][ind2];
                    }
                }

                E_delta += 0.5 * uvalue[iat][l] * (noccup - trace);
            }
        }
    }

    ModuleBase::timer::tick ("LCAO_DftU_New","cal_gedm");
}

void LCAO_DftU_New::check_gedm()
{
    std::ofstream ofs("gedm.dat");
    for(int is = 0; is < GlobalV::NSPIN; is ++)
    {
        for(int inl=0;inl<inlmax;inl++)
        {
            int nm = 2 * inl_l[inl] + 1;
            for (int m1 = 0;m1 < nm;++m1)
            {
                for (int m2 = 0;m2 < nm;++m2)
                {
                    int index = m1 * nm + m2;
                    //*2 is for Hartree to Ry
                    ofs << this->gedm[is][inl][index] << " ";
                }
            }   
            ofs << std::endl;     
        }
    }
}