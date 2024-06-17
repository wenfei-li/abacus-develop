#ifndef ESOLVER_KS_H
#define ESOLVER_KS_H
#include <string.h>

#include <fstream>

#include "esolver_fp.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_elecstate/module_charge/charge_extra.h"
#include "module_elecstate/module_charge/charge_mixing.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_hsolver/hsolver.h"
#include "module_psi/psi.h"
#include "module_io/cal_test.h"
#include "module_io/output_potential.h"
#include "module_io/output_rho.h"
namespace ModuleESolver
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class ESolver_KS : public ESolver_FP
{
	public:

        //! Constructor
		ESolver_KS();

        //! Deconstructor
		virtual ~ESolver_KS();

		double scf_thr;   // scf threshold

		double drho;      // the difference between rho_in (before HSolver) and rho_out (After HSolver)

		int maxniter;     // maximum iter steps for scf

		int niter;        // iter steps actually used in scf

		bool conv_elec;   // If electron density is converged in scf.

		int out_freq_elec;// frequency for output

		virtual void before_all_runners(Input& inp, UnitCell& cell) override;

		virtual void init_after_vc(Input& inp, UnitCell& cell) override;    // liuyu add 2023-03-09

		virtual void runner(const int istep, UnitCell& cell) override;

		// calculate electron density from a specific Hamiltonian
		virtual void hamilt2density(const int istep, const int iter, const double ethr);

		// calculate electron states from a specific Hamiltonian
		virtual void hamilt2estates(const double ethr){};

		// get current step of Ionic simulation
		virtual int get_niter() override;

		// get maxniter used in current scf
		virtual int get_maxniter() override;

		// get conv_elec used in current scf
		virtual bool get_conv_elec() override;

	protected:
		//! Something to do before SCF iterations.
		virtual void before_scf(const int istep) {};

		//! Something to do before hamilt2density function in each iter loop.
		virtual void iter_init(const int istep, const int iter) {};

		//! Something to do after hamilt2density function in each iter loop.
		virtual void iter_finish(const int iter) {};

		//! Something to do after SCF iterations when SCF is converged or comes to the max iter step.
		virtual void after_scf(const int istep) {};

		//! <Temporary> It should be replaced by a function in Hamilt Class
		virtual void update_pot(const int istep, const int iter) {};

		//! choose strategy when charge density convergence achieved
		virtual bool do_after_converge(int& iter){return true;}

	protected:

		// Print the headline on the screen:
		// ITER   ETOT(eV)       EDIFF(eV)      DRHO    TIME(s) 
		void print_head();

		// Print inforamtion in each iter
		// G1    -3.435545e+03  0.000000e+00   3.607e-01  2.862e-01
		// for metaGGA
		// ITER   ETOT(eV)       EDIFF(eV)      DRHO       DKIN       TIME(s) 
		// G1    -3.435545e+03  0.000000e+00   3.607e-01  3.522e-01  2.862e-01
		void print_iter(
				const int iter, 
				const double drho, 
				const double dkin, 
				const double duration, 
				const double ethr);

		// Write the headline in the running_log file
		// "PW/LCAO" ALGORITHM --------------- ION=   1  ELEC=   1--------------------------------
		void write_head(
				std::ofstream& ofs_running, 
				const int istep, 
				const int iter);

		/// @brief create a new ModuleIO::Output_Rho object to output charge density
		ModuleIO::Output_Rho create_Output_Rho(int is, int iter, const std::string& prefix="None");

		/// @brief create a new ModuleIO::Output_Rho object to print kinetic energy density
		ModuleIO::Output_Rho create_Output_Kin(int is, int iter, const std::string& prefix = "None");

		/// @brief create a new ModuleIO::Output_Potential object to print potential
		ModuleIO::Output_Potential create_Output_Potential(int iter, const std::string& prefix = "None");

        //! Solve Hamitonian
		hsolver::HSolver<T, Device>* phsol = nullptr;

        //! Hamiltonian
		hamilt::Hamilt<T, Device>* p_hamilt = nullptr;

		ModulePW::PW_Basis_K* pw_wfc = nullptr;

		Charge_Mixing* p_chgmix = nullptr;

		wavefunc wf;

		Charge_Extra CE;

		// wavefunction coefficients
		psi::Psi<T>* psi = nullptr;

	protected:

		std::string basisname; //PW or LCAO

		void print_wfcfft(Input& inp, std::ofstream &ofs);
};	
} // end of namespace
#endif
