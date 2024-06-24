#ifndef SPARSE_FORMAT_ST_H
#define SPARSE_FORMAT_ST_H

#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_HS_arrays.hpp"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"

namespace sparse_format
{
//! calculate overlap matrix with lattice vector R
void cal_SR(const Parallel_Orbitals& pv,
            std::set<Abfs::Vector3_Order<int>>& all_R_coor,
            std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>& SR_sparse,
            std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>& SR_soc_sparse,
            Grid_Driver& grid,
            const double& sparse_thr,
            hamilt::Hamilt<std::complex<double>>* p_ham);

//! calculate kinetic matrix with lattice vector R
void cal_TR(const UnitCell& ucell,
            const Parallel_Orbitals& pv,
            LCAO_Matrix& lm,
            LCAO_HS_Arrays& HS_arrays,
            Grid_Driver& grid,
            const TwoCenterBundle& two_center_bundle,
            const double& sparse_thr);

//! cal_STN_R_for_T is only called by cal_TR
void cal_STN_R_for_T(const UnitCell& ucell,
                     const Parallel_Orbitals& pv,
                     LCAO_Matrix& lm,
                     LCAO_HS_Arrays& HS_arrays,
                     Grid_Driver& grid,
                     const double& sparse_thr);
} // namespace sparse_format

#endif
