#ifndef LCAO_DFTU_NEW_H
#define LCAO_DFTU_NEW_H

#include "module_base/intarray.h"
#include "module_base/complexmatrix.h"
#include "module_basis/module_ao/ORB_gen_tables.h"
#include <unordered_map>

#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_base/intarray.h"
#include "module_base/complexmatrix.h"
#include "module_io/winput.h"
#include "module_base/matrix.h"
#include "module_base/timer.h"

// Wenfei 2023/07/25
// After noting the form of DeePKS operator is the same as that of DFT+U,
// I am going to add a mode to use DeePKS as the DFT+U correction term.
// In fact the only major change I need to make is on cal_gedm, which calculates
// the correction potential, but of course there will be some minor changes as well.
// The values of gedm are obtained as the derivative of energy w.r.t. descriptors in DeePKS,
// while in DFT+U calculated as: (U-J)(1/2-n), where n is the occupation matrix.
// It should be noted that in DeePKS, the same projectors are used on both sides,
// corresponding to the so-called 'full' representation of occupation matrix calculation,
// while the implementation in module_dftu uses the 'dual' representation.
// Also, two potential perils: 1. DeePKS is applied on every l channel of every atom, so there will be
// some unnecessary operations in the loop (but can be improved by adding some skipping conditions later)
// 2. the DeePKS projectors are the same for all atoms. To modify it will be lots of work,
// so I intend to let it be and see how it goes.

class LCAO_Deepks
{

//-------------------
// public variables
//-------------------
public:
    
    ///(Unit: Ry) Correction energy provided by NN
    double E_delta = 0.0;
    ///(Unit: Ry)  \f$tr(\rho H_\delta), \rho = \sum_i{c_{i, \mu}c_{i,\nu}} \f$ (for gamma_only)
    double e_delta_band = 0.0;

    ///Correction term to the Hamiltonian matrix: \f$\langle\psi|V_\delta|\psi\rangle\f$ (for gamma only)
    double* H_V_delta;
    ///Correction term to Hamiltonian, for multi-k
    ///In R space:
    double* H_V_deltaR;
    ///In k space:
    std::complex<double>** H_V_delta_k;

//-------------------
// private variables
//-------------------
private:

	int lmaxd = 0; //max l of descirptors
	int nmaxd = 0; //#. descriptors per l
	int inlmax = 0; //tot. number {i,n,l} - atom, n, l

    bool init_pdm = false; //for DeePKS NSCF calculation

    // saves <psi(0)|alpha(R)>, for gamma only
    std::vector<std::vector<std::unordered_map<int,std::vector<std::vector<double>>>>> nlm_save;
    
    // saves <psi(0)|alpha(R)>, for multi_k
    typedef std::tuple<int,int,int,int> key_tuple;
    std::vector<std::map<key_tuple,std::unordered_map<int,std::vector<std::vector<double>>>>> nlm_save_k;

    // projected density matrix
	double** pdm;	//[tot_Inl][2l+1][2l+1]	caoyu modified 2021-05-07

	// +U correction strength
	double** gedm;	//[tot_Inl][2l+1][2l+1]

	ModuleBase::IntArray* inl_index;	//caoyu add 2021-05-07
	int* inl_l;	//inl_l[inl_index] = l of descriptor with inl_index

    const Parallel_Orbitals* pv;

//-------------------
// subroutines, grouped according to the file they are in:
//-------------------

//-------------------
// LCAO_dftu_new.cpp
//-------------------

//This file contains constructor and destructor of the class LCAO_dftu_new, 
//as well as subroutines for initializing and releasing relevant data structures 

//Other than the constructor and the destructor, it contains 3 types of subroutines:
//1. subroutines that are related to calculating descriptors:
//  - init : allocates some arrays
//  - init_index : records the index (inl)
//  - allocate_nlm : allocates data structures (nlm_save) which is used to store <chi|alpha>
//2. subroutines that are related to calculating force label:
//  - init_gdmx : allocates gdmx; it is a private subroutine
//  - del_gdmx : releases gdmx
//3. subroutines that are related to V_delta:
//  - allocate_V_delta : allocates H_V_delta; if calculating force, it also calls
//      init_gdmx, as well as allocating F_delta
//  - allocate_V_deltaR : allcoates H_V_deltaR, for multi-k calculations

public:

    explicit LCAO_Deepks();
    ~LCAO_Deepks();

    ///Allocate memory and calculate the index of descriptor in all atoms. 
    ///(only for descriptor part, not including scf)
    void init(const LCAO_Orbitals &orb,
        const int nat,
        const int ntype,
        const Parallel_Orbitals& pv_in,
        std::vector<int> na);

    ///Allocate memory for correction to Hamiltonian
    void allocate_V_delta(const int nat, const int nks = 1);
    void allocate_V_deltaR(const int nnr);

    // array for storing gdmx, used for calculating gvx
	void init_gdmx(const int nat);
	void del_gdmx(const int nat);

    // array for storing gdm_epsl, used for calculating gvx
	void init_gdmepsl();
	void del_gdmepsl();

private:

    // arrange index of descriptor in all atoms
	void init_index(const int ntype,
        const int nat,
        std::vector<int> na,
        const int tot_inl,
        const LCAO_Orbitals &orb);
    // data structure that saves <psi|alpha>
    void allocate_nlm(const int nat);

//-------------------
// LCAO_dftu_new_psialpha.cpp
//-------------------

//wenfei 2022-1-5
//This file contains 2 subroutines:
//1. build_psialpha, which calculates the overlap
//between atomic basis and projector alpha : <psi_mu|alpha>
//which will be used in calculating pdm, gdmx, H_V_delta, F_delta;
//2. check_psialpha, which prints the results into .dat files
//for checking

public:

    //calculates <chi|alpha>
    void build_psialpha(const bool& cal_deri/**< [in] 0 for 2-center intergration, 1 for its derivation*/,
        const UnitCell &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver& GridD,
        const ORB_gen_tables &UOT);

    void check_psialpha(const bool& cal_deri/**< [in] 0 for 2-center intergration, 1 for its derivation*/,
        const UnitCell &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver& GridD,
        const ORB_gen_tables &UOT);

//-------------------
// LCAO_dftu_new_pdm.cpp
//-------------------

//This file contains subroutines for calculating pdm,
//which is defind as sum_mu,nu rho_mu,nu (<chi_mu|alpha><alpha|chi_nu>);
//as well as gdmx, which is the gradient of pdm, defined as
//sum_mu,nu rho_mu,nu d/dX(<chi_mu|alpha><alpha|chi_nu>)

//It also contains subroutines for printing pdm and gdmx
//for checking purpose

//There are 6 subroutines in this file:
//1. cal_projected_DM, which is used for calculating pdm for gamma point calculation
//2. cal_projected_DM_k, counterpart of 1, for multi-k
//3. check_projected_dm, which prints pdm to descriptor.dat

//4. cal_gdmx, calculating gdmx (and optionally gdm_epsl for stress) for gamma point
//5. cal_gdmx_k, counterpart of 3, for multi-k
//6. check_gdmx, which prints gdmx to a series of .dat files

public:

    ///calculate projected density matrix:
    ///pdm = sum_i,occ <phi_i|alpha1><alpha2|phi_k>
    void cal_projected_DM(const ModuleBase::matrix& dm/**< [in] density matrix*/,
        const UnitCell &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver& GridD);
    void cal_projected_DM_k(const std::vector<ModuleBase::ComplexMatrix>& dm,
        const UnitCell &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver& GridD,
        const int nks,
        const std::vector<ModuleBase::Vector3<double>> &kvec_d);
    void check_projected_dm(void);

//-------------------
// LCAO_dftu_new_vdelta.cpp
//-------------------

//This file contains subroutines related to V_delta, which is the deepks contribution to Hamiltonian
//defined as |alpha>V(D)<alpha|
//as well as subroutines for printing them for checking
//It also contains subroutine related to calculating e_delta_bands, which is basically
//tr (rho * V_delta)

//Four subroutines are contained in the file:
//1. add_v_delta : adds deepks contribution to hamiltonian, for gamma only
//2. add_v_delta_k : counterpart of 1, for multi-k
//3. check_v_delta : prints H_V_delta for checking
//4. check_v_delta_k : prints H_V_deltaR for checking
//5. cal_e_delta_band : calculates e_delta_bands for gamma only
//6. cal_e_delta_band_k : counterpart of 4, for multi-k

public:

    ///add dV to the Hamiltonian matrix
    void add_v_delta(const UnitCell &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver& GridD);
    void add_v_delta_k(const UnitCell &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver& GridD,
        const int nnr_in);

    void check_v_delta();
    void check_v_delta_k(const int nnr);

    ///calculate tr(\rho V_delta)
    void cal_e_delta_band(const std::vector<ModuleBase::matrix>& dm/**<[in] density matrix*/);
    void cal_e_delta_band_k(const std::vector<ModuleBase::ComplexMatrix>& dm/**<[in] density matrix*/,
        const int nks);

//-------------------
// LCAO_dftu_new_fdelta.cpp
//-------------------

//This file contains subroutines for calculating F_delta,
//which is defind as sum_mu,nu rho_mu,nu d/dX (<chi_mu|alpha>V(D)<alpha|chi_nu>)

//There are 3 subroutines in this file:
//1. cal_f_delta_gamma, which is used for gamma point calculation
//2. cal_f_delta_k, which is used for multi-k calculation
//3. check_f_delta, which prints F_delta into F_delta.dat for checking

    //for gamma only, pulay and HF terms of force are calculated together
    void cal_f_delta_gamma(const ModuleBase::matrix& dm/**< [in] density matrix*/,
        const UnitCell &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver& GridD,
        const bool isstress, ModuleBase::matrix& svnl_dalpha);

    //for multi-k, pulay and HF terms of force are calculated together
    void cal_f_delta_k(const std::vector<ModuleBase::ComplexMatrix>& dm/**<[in] density matrix*/,
        const UnitCell &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver& GridD,
        const int nks,
        const std::vector<ModuleBase::Vector3<double>> &kvec_d,
        const bool isstress, ModuleBase::matrix& svnl_dalpha);

    void check_f_delta(const int nat, ModuleBase::matrix& svnl_dalpha);

//-------------------
// LCAO_dftu_new_io.cpp
//-------------------

//This file contains subroutines that contains interface with libnpy
//since many arrays must be saved in numpy format
//It also contains subroutines for printing density matrices
//which is used in unit tests

//There are 2 subroutines for printing density matrices:
//1. print_dm : for gamma only
//2. print_dm_k : for multi-k

//And 6 which prints quantities in .npy format 
//3. save_npy_d : descriptor ->dm_eig.npy
//4. save_npy_gvx : gvx ->grad_vx.npy
//5. save_npy_e : energy
//6. save_npy_f : force
//7. save_npy_s : stress
//8. save_npy_o: orbital
//9. save_npy_orbital_precalc: orbital_precalc -> orbital_precalc.npy

public:
  
    ///print density matrices
    void print_dm(const ModuleBase::matrix &dm);
    void print_dm_k(const int nks, const std::vector<ModuleBase::ComplexMatrix>& dm);

//-------------------
// LCAO_dftu_new_mpi.cpp
//-------------------

//This file contains only one subroutine, allsum_deepks
//which is used to perform allsum on a two-level pointer
//It is used in a few places in the deepks code

#ifdef __MPI

public:
    //reduces a dim 2 array
    void allsum_deepks(
        int inlmax, //first dimension
        int ndim, //second dimension
        double** mat); //the array being reduced 
#endif

};

namespace GlobalC
{
    extern LCAO_Deepks ld;
}

#endif
