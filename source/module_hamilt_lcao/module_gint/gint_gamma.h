//=========================================================
//AUTHOR : mohan
//DATE : 2009-09-16
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#ifndef GINT_GAMMA_H
#define GINT_GAMMA_H
#include "gint.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "grid_technique.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#ifdef _OPENMP
#include <omp.h>
#endif

//=========================================================
// ModuleBase::Integral On 3D Grids, different from Grid_Integral
// Feature : Matrix Elements Of Local Potential For 
// Numerical Orbitals
//=========================================================

class Gint_Gamma : public Gint
{
	public:

	// Gint_Gamma();
	// ~Gint_Gamma();

    //------------------------------------------------------
    // in gint_gamma_vl.cpp 
    //------------------------------------------------------
    // there is an additional step in calculating vlocal for gamma point
    // namely the redistribution of Hamiltonian from grid to 2D block format
    // hence we have an additional layer outside the unified interface
    void cal_vlocal(Gint_inout* inout, const bool new_e_iteration);

    //------------------------------------------------------
    // in gint_gamma_env.cpp 
    //------------------------------------------------------
	// calcualte the envelope function
	void cal_env(const double* wfc, double* rho,UnitCell &ucell);

    //------------------------------------------------------
    // in veff_lcao.cpp
    //------------------------------------------------------
    /// transfer this->hRGint to Veff::hR
    void transfer_pvpR(hamilt::HContainer<double>* hR,const UnitCell* ucell);

private:

    double*** DM = nullptr;   //pointer to LOC.DM

};

#endif
