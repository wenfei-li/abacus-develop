#ifndef GRID_TECHNIQUE_H
#define GRID_TECHNIQUE_H

#include "grid_index.h"
#include "grid_meshball.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include  "module_cell/unitcell.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#if ((defined __CUDA) /* || (defined __ROCM) */)
#include <cuda_runtime.h>

#include "kernels/cuda/cuda_tools.cuh"
#include "kernels/cuda/gemm_selector.cuh"
#endif

// Author: mohan
// Date: 2009-10-17
class Grid_Technique : public Grid_MeshBall
{
    // public variables.
  public:
    Grid_Technique();
    ~Grid_Technique();
    //------------------------------------
    // 1: Info about atom number on grid.
    //------------------------------------
    // record how many atoms on each grid.
    std::vector<int> how_many_atoms;
    // max atom on grid
    int max_atom;
    // sum of how_many_atoms
    int total_atoms_on_grid;
    std::vector<int> start_ind;

    //------------------------------------
    // 2: Info about which atom on grid.
    //------------------------------------
    // save the start position of each big cell's adjacent
    // atoms in 1D grid.
    std::vector<int> bcell_start;
    // save the 'iat' atom.
    // dim: total_atoms_on_grid.
    std::vector<int> which_atom;

    //--------------------------------------
    // save the bigcell index in meshball.
    // dim: total_atoms_on_grid.
    //--------------------------------------
    std::vector<int> which_bigcell;
    std::vector<int> which_unitcell;

    //------------------------------------
    // 3: which atom on local grid.
    //------------------------------------
    int lnat;      // local nat.  
    int lgd;       // local grid dimension.  lgd * lgd symmetry matrix.
    std::vector<bool> in_this_processor;
    std::vector<int> trace_iat;
    std::vector<int> trace_lo; // trace local orbital.

    //---------------------------------------
    // nnrg: number of matrix elements on
    // each processor's real space grid.
    // use: GridT.in_this_processor
    //---------------------------------------
    int nnrg;
    bool allocate_find_R2;
    std::vector<int> nlocdimg;
    std::vector<int> nlocstartg;
    std::vector<int> nad; // number of adjacent atoms for each atom.
    std::vector<std::vector<int>> find_R2;
    std::vector<std::vector<int>> find_R2_sorted_index;
    std::vector<std::vector<int>> find_R2st;

    int binary_search_find_R2_offset(int val, int iat) const;

    //UnitCell and LCAO_Obrbitals
    const UnitCell* ucell;
    const LCAO_Orbitals* orb;

    // UnitCell parameters
    int nwmax;
    int nr_max;
    int ntype;

    // LCAO Orbitals
    double dr_uniform;
		std::vector<double> rcuts;
    std::vector<std::vector<double>> psi_u;
    std::vector<std::vector<double>> dpsi_u;
    std::vector<std::vector<double>> d2psi_u;

    // indexes for nnrg -> orbital index + R index
    std::vector<gridIntegral::gridIndex> nnrg_index;

    void set_pbc_grid(const int& ncx_in, const int& ncy_in, const int& ncz_in, const int& bx_in, const int& by_in,
                      const int& bz_in, const int& nbx_in, const int& nby_in, const int& nbz_in, const int& nbxx_in,
                      const int& nbzp_start_in, const int& nbzp_in, const int& ny, const int& nplane,
                      const int& startz_current, const UnitCell& ucell, const double& dr_uniform, 
                      const std::vector<double>& rcuts,const std::vector<std::vector<double>>& psi_u,
                      const std::vector<std::vector<double>>& dpsi_u,const std::vector<std::vector<double>>& d2psi_u,
                      const int& num_stream);

    /// number of elements(basis-pairs) in this processon
    /// on all adjacent atoms-pairs(Grid division)
    void cal_nnrg(Parallel_Orbitals* pv);
    int cal_RindexAtom(const int& u1,
                       const int& u2,
                       const int& u3,
                       const int& iat2) const;

  private:
    void cal_max_box_index(void);

    int maxB1;
    int maxB2;
    int maxB3;

    int minB1;
    int minB2;
    int minB3;

    int nB1;
    int nB2;
    int nB3;

    int nbox;

    // atoms on meshball
    void init_atoms_on_grid(const int& ny,
                            const int& nplane,
                            const int& startz_current,
                            const UnitCell& ucell);
    void init_atoms_on_grid2(const int* index2normal,
                            const UnitCell& ucell);
    void cal_grid_integration_index(void);
    void cal_trace_lo(const UnitCell& ucell);
    void check_bigcell(int* ind_bigcell, char* bigcell_on_processor);
    void get_startind(const int& ny,
                      const int& nplane,
                      const int& startz_current);

#if ((defined __CUDA) /* || (defined __ROCM) */)
  public:
    double* ylmcoef_g;
    bool is_malloced;

    int* atom_nw_g;
    int* atom_nwl_g;
    double* psi_u_g;
    bool* atom_new_g;
    int* atom_ylm_g;
    int* atom_l_g;

    int nstreams = 4;
    // streams[nstreams]
    // TODO it needs to be implemented through configuration files
    matrix_multiple_func_type fastest_matrix_mul;

  private:
    void init_gpu_gint_variables(const UnitCell& ucell, const int num_stream);
    void free_gpu_gint_variables(int nat);

#endif
};
#endif
