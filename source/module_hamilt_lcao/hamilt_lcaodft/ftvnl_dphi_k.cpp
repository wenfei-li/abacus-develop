#include "FORCE.h"

#include <map>
#include <unordered_map>

#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/cal_dm.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/write_HS.h"

#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


template<>
void Force_LCAO<std::complex<double>>::cal_ftvnl_dphi(
    const elecstate::DensityMatrix<std::complex<double>, double>* dm,
    const Parallel_Orbitals& pv,
    const UnitCell& ucell,
    ForceStressArrays& fsr,
    const bool isforce,
    const bool isstress,
    ModuleBase::matrix& ftvnl_dphi,
    ModuleBase::matrix& stvnl_dphi,
    Record_adj* ra)
{
    ModuleBase::TITLE("Force_LCAO", "cal_ftvnl_dphi");
    ModuleBase::timer::tick("Force_LCAO", "cal_ftvnl_dphi");

    const int nspin_DMR = (GlobalV::NSPIN == 2) ? 2 : 1;

    int total_irr = 0;
    // get the adjacent atom's information.

#ifdef _OPENMP
#pragma omp parallel
    {
        int num_threads = omp_get_num_threads();
        ModuleBase::matrix local_stvnl_dphi(3, 3);
        int local_total_irr = 0;
#pragma omp for schedule(dynamic)
#else
		ModuleBase::matrix& local_stvnl_dphi = stvnl_dphi;
		int& local_total_irr = total_irr;
#endif
		for (int iat = 0; iat < ucell.nat; iat++)
        {
            const int T1 = ucell.iat2it[iat];
            Atom* atom1 = &ucell.atoms[T1];
            const int I1 = ucell.iat2ia[iat];
            // get iat1
            int iat1 = ucell.itia2iat(T1, I1);
            //
            int irr = pv.nlocstart[iat];
            const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
            double* ftvnl_dphi_iat;
			if (isforce)
			{
				ftvnl_dphi_iat = &ftvnl_dphi(iat, 0);
			}
#ifdef _OPENMP
            // using local stack to avoid false sharing in multi-threaded case
            double ftvnl_dphi_temp[3] = {0.0, 0.0, 0.0};
            if (num_threads > 1)
            {
                ftvnl_dphi_iat = ftvnl_dphi_temp;
            }
#endif
            for (int cb = 0; cb < ra->na_each[iat]; ++cb)
            {
                const int T2 = ra->info[iat][cb][3];
                const int I2 = ra->info[iat][cb][4];
                const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                Atom* atom2 = &ucell.atoms[T2];
                // get iat2
                int iat2 = ucell.itia2iat(T2, I2);
                double Rx = ra->info[iat][cb][0];
                double Ry = ra->info[iat][cb][1];
                double Rz = ra->info[iat][cb][2];
                // get BaseMatrix
                if (pv.get_row_size(iat1) <= 0 || pv.get_col_size(iat2) <= 0)
                {
                    continue;
                }
                std::vector<hamilt::BaseMatrix<double>*> tmp_matrix;
                for (int is = 0; is < nspin_DMR; ++is)
                {
                    tmp_matrix.push_back(dm->get_DMR_pointer(is+1)->find_matrix(iat1, iat2, Rx, Ry, Rz));
                }

                for (int mu = 0; mu < pv.get_row_size(iat1); ++mu)
                {
                    for (int nu = 0; nu < pv.get_col_size(iat2); ++nu)
                    {
                        // get value from dm
                        double dm2d1 = 0.0;
                        for (int is = 0; is < nspin_DMR; ++is)
                        {
                            dm2d1 += tmp_matrix[is]->get_value(mu, nu);
                        }
                        double dm2d2 = 2.0 * dm2d1;
                        //
                        if (isforce)
                        {
                            ftvnl_dphi_iat[0] += dm2d2 * fsr.DHloc_fixedR_x[irr];
                            ftvnl_dphi_iat[1] += dm2d2 * fsr.DHloc_fixedR_y[irr];
                            ftvnl_dphi_iat[2] += dm2d2 * fsr.DHloc_fixedR_z[irr];
                        }
                        if (isstress)
                        {
                            local_stvnl_dphi(0, 0) -= dm2d1 * fsr.stvnl11[irr];
                            local_stvnl_dphi(0, 1) -= dm2d1 * fsr.stvnl12[irr];
                            local_stvnl_dphi(0, 2) -= dm2d1 * fsr.stvnl13[irr];
                            local_stvnl_dphi(1, 1) -= dm2d1 * fsr.stvnl22[irr];
                            local_stvnl_dphi(1, 2) -= dm2d1 * fsr.stvnl23[irr];
                            local_stvnl_dphi(2, 2) -= dm2d1 * fsr.stvnl33[irr];
                        }
                        //}
                        ++local_total_irr;
                        ++irr;
                    } // end kk
                }     // end jj
            }         // end cb
#ifdef _OPENMP
            if (isforce && num_threads > 1)
            {
                ftvnl_dphi(iat, 0) += ftvnl_dphi_iat[0];
                ftvnl_dphi(iat, 1) += ftvnl_dphi_iat[1];
                ftvnl_dphi(iat, 2) += ftvnl_dphi_iat[2];
            }
#endif
        } // end iat
#ifdef _OPENMP
#pragma omp critical(cal_ftvnl_dphi_k_reduce)
        {
            total_irr += local_total_irr;
            if (isstress)
            {
                stvnl_dphi(0, 0) += local_stvnl_dphi(0, 0);
                stvnl_dphi(0, 1) += local_stvnl_dphi(0, 1);
                stvnl_dphi(0, 2) += local_stvnl_dphi(0, 2);
                stvnl_dphi(1, 1) += local_stvnl_dphi(1, 1);
                stvnl_dphi(1, 2) += local_stvnl_dphi(1, 2);
                stvnl_dphi(2, 2) += local_stvnl_dphi(2, 2);
            }
        }
    }
#endif
    assert(total_irr == pv.nnr);

    if (isstress)
    {
        StressTools::stress_fill(ucell.lat0, ucell.omega, stvnl_dphi);
    }

    ModuleBase::timer::tick("Force_LCAO", "cal_ftvnl_dphi");
    return;
}
