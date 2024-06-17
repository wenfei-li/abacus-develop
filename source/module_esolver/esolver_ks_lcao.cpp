#include "esolver_ks_lcao.h"

#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_io/dos_nao.h"
#include "module_io/nscf_band.h"
#include "module_io/output_dmk.h"
#include "module_io/output_log.h"
#include "module_io/output_mulliken.h"
#include "module_io/output_sk.h"
#include "module_io/to_qo.h"
#include "module_io/write_HS.h"
#include "module_io/write_Vxc.hpp"
#include "module_io/write_istate_info.h"
#include "module_io/write_proj_band_lcao.h"

//--------------temporary----------------------------
#include "module_base/global_function.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"

#include <memory>
#ifdef __EXX
#include "module_ri/RPA_LRI.h"
#endif

#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks_interface.h"
#endif
//-----force& stress-------------------
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"

//-----HSolver ElecState Hamilt--------
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hsolver/hsolver_lcao.h"
// function used by deepks
#include "module_elecstate/cal_dm.h"
//---------------------------------------------------

#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_io/write_wfc_nao.h"

namespace ModuleESolver
{

//------------------------------------------------------------------------------
//! the 1st function of ESolver_KS_LCAO: constructor
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
ESolver_KS_LCAO<TK, TR>::ESolver_KS_LCAO()
{
    this->classname = "ESolver_KS_LCAO";
    this->basisname = "LCAO";

// the following EXX code should be removed to other places, mohan 2024/05/11
#ifdef __EXX
    if (GlobalC::exx_info.info_ri.real_number)
    {
        this->exx_lri_double = std::make_shared<Exx_LRI<double>>(GlobalC::exx_info.info_ri);
        this->exd = std::make_shared<Exx_LRI_Interface<TK, double>>(this->exx_lri_double);
        this->LM.Hexxd = &this->exd->get_Hexxs();
    }
    else
    {
        this->exx_lri_complex = std::make_shared<Exx_LRI<std::complex<double>>>(GlobalC::exx_info.info_ri);
        this->exc = std::make_shared<Exx_LRI_Interface<TK, std::complex<double>>>(this->exx_lri_complex);
        this->LM.Hexxc = &this->exc->get_Hexxs();
    }
#endif
}

//------------------------------------------------------------------------------
//! the 2nd function of ESolver_KS_LCAO: deconstructor
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
ESolver_KS_LCAO<TK, TR>::~ESolver_KS_LCAO()
{
#ifndef USE_NEW_TWO_CENTER
    this->orb_con.clear_after_ions(*uot_, GlobalC::ORB, GlobalV::deepks_setorb, GlobalC::ucell.infoNL.nproj);
#endif
    delete uot_;
}

//------------------------------------------------------------------------------
//! the 3rd function of ESolver_KS_LCAO: init
//! mohan add 2024-05-11
//! 1) calculate overlap matrix S or initialize
//! 2) init ElecState
//! 3) init LCAO basis
//! 4) redundant ParaV and LM pointers
//! 5) initialize Density Matrix
//! 6) initialize Hamilt in LCAO
//! 7) initialize exx
//! 8) Quxin added for DFT+U
//! 9) ppcell
//! 10) init HSolver
//! 11) inititlize the charge density.
//! 12) initialize the potential.
//! 13) initialize deepks
//! 14) set occupations?
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_all_runners(Input& inp, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_all_runners");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_all_runners");

    // 1) calculate overlap matrix S
    if (GlobalV::CALCULATION == "get_S")
    {
        // 1.1) read pseudopotentials
        ucell.read_pseudo(GlobalV::ofs_running);

        // 1.2) symmetrize things
        if (ModuleSymmetry::Symmetry::symm_flag == 1)
        {
            ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
        }

        // 1.3) Setup k-points according to symmetry.
        this->kv
            .set(ucell.symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, ucell.G, ucell.latvec, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

        Print_Info::setup_parameters(ucell, this->kv);
    }
    else
    {
        // 1) else, call before_all_runners() in ESolver_KS
        ESolver_KS<TK>::before_all_runners(inp, ucell);
    } // end ifnot get_S

    // 2) init ElecState
    // autoset nbands in ElecState, it should before basis_init (for Psi 2d divid)
    if (this->pelec == nullptr)
    {
        // TK stands for double and complex<double>?
        this->pelec = new elecstate::ElecStateLCAO<TK>(&(this->chr), // use which parameter?
                                                       &(this->kv),
                                                       this->kv.get_nks(),
                                                       &(this->LOC),  // use which parameter?
                                                       &(this->GG),   // mohan add 2024-04-01
                                                       &(this->GK),   // mohan add 2024-04-01
                                                       &(this->LOWF), // use which parameter?
                                                       this->pw_rho,
                                                       this->pw_big);
    }

    // 3) init LCAO basis
    // reading the localized orbitals/projectors
    // construct the interpolation tables.
    this->init_basis_lcao(this->orb_con, inp, ucell);
    //------------------init Basis_lcao----------------------

    //! pass basis-pointer to EState and Psi
    /*
    Inform: on getting rid of ORB_control and Parallel_Orbitals

    Have to say it is all the stories start, the ORB_control instance pass its Parallel_Orbitals instance to
    Local_Orbital_Charge, Local_Orbital_Wfc and LCAO_Matrix, which is actually for getting information
    of 2D block-cyclic distribution.

    To remove LOC, LOWF and LM use in functions, one must make sure there is no more information imported
    to those classes. Then places where to get information from them can be substituted to orb_con

    Plan:
    1. Specifically for paraV, the thing to do first is to replace the use of ParaV to the oen of ORB_control.
    Then remove ORB_control and place paraV somewhere.
    */
    this->LOC.ParaV = &(this->orb_con.ParaV);
    this->LOWF.ParaV = &(this->orb_con.ParaV);
    this->LM.ParaV = &(this->orb_con.ParaV);

    // 5) initialize density matrix
    dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)
        ->init_DM(&this->kv, &(this->orb_con.ParaV), GlobalV::NSPIN);

    // this function should be removed outside of the function
    if (GlobalV::CALCULATION == "get_S")
    {
        ModuleBase::timer::tick("ESolver_KS_LCAO", "init");
        return;
    }

    // 6) initialize Hamilt in LCAO
    // * allocate H and S matrices according to computational resources
    // * set the 'trace' between local H/S and global H/S
    this->LM.divide_HS_in_frag(GlobalV::GAMMA_ONLY_LOCAL, orb_con.ParaV, this->kv.get_nks());

#ifdef __EXX
    // 7) initialize exx
    // PLEASE simplify the Exx_Global interface
    if (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax"
        || GlobalV::CALCULATION == "md")
    {
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            /* In the special "two-level" calculation case,
            first scf iteration only calculate the functional without exact exchange.
            but in "nscf" calculation, there is no need of "two-level" method. */
            if (ucell.atoms[0].ncpp.xc_func == "HF" || ucell.atoms[0].ncpp.xc_func == "PBE0"
                || ucell.atoms[0].ncpp.xc_func == "HSE")
            {
                XC_Functional::set_xc_type("pbe");
            }
            else if (ucell.atoms[0].ncpp.xc_func == "SCAN0")
            {
                XC_Functional::set_xc_type("scan");
            }

            // GlobalC::exx_lcao.init();
            if (GlobalC::exx_info.info_ri.real_number)
            {
                this->exx_lri_double->init(MPI_COMM_WORLD, this->kv);
            }
            else
            {
                this->exx_lri_complex->init(MPI_COMM_WORLD, this->kv);
            }
        }
    }
#endif

    // 8) initialize DFT+U
    if (GlobalV::dft_plus_u)
    {
        GlobalC::dftu.init(ucell, this->LM, this->kv.get_nks());
    }

    // 9) initialize ppcell
    GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc, this->pw_rho);

    // 10) initialize the HSolver
    if (this->phsol == nullptr)
    {
        this->phsol = new hsolver::HSolverLCAO<TK>(&(this->orb_con.ParaV));
        this->phsol->method = GlobalV::KS_SOLVER;
    }

    // 11) inititlize the charge density
    this->pelec->charge->allocate(GlobalV::NSPIN);
    this->pelec->omega = GlobalC::ucell.omega;

    // 12) initialize the potential
    if (this->pelec->pot == nullptr)
    {
        this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                    this->pw_rho,
                                                    &GlobalC::ucell,
                                                    &(GlobalC::ppcell.vloc),
                                                    &(this->sf),
                                                    &(this->pelec->f_en.etxc),
                                                    &(this->pelec->f_en.vtxc));
    }

#ifdef __DEEPKS
    // 13) initialize deepks
    if (GlobalV::deepks_scf)
    {
        // load the DeePKS model from deep neural network
        GlobalC::ld.load_model(INPUT.deepks_model);
    }
#endif

    // 14) set occupations
    if (GlobalV::ocp)
    {
        this->pelec->fixed_weights(GlobalV::ocp_kb, GlobalV::NBANDS, GlobalV::nelec);
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_all_runners");
    return;
}

//------------------------------------------------------------------------------
//! the 4th function of ESolver_KS_LCAO: init_after_vc
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::init_after_vc(Input& inp, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "init_after_vc");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "init_after_vc");

    ESolver_KS<TK>::init_after_vc(inp, ucell);
    /*
    Notes: on the removal of LOWF
    Following constructor of ElecStateLCAO requires LOWF. However, ElecState only need
    LOWF to do wavefunction 2dbcd (2D BlockCyclicDistribution) gathering. So, a free
    function is needed to replace the use of LOWF. The function indeed needs the information
    about 2dbcd, therefore another instance storing the information is needed instead.
    Then that instance will be the input of "the free function to gather".
    */
    if (GlobalV::md_prec_level == 2)
    {
        delete this->pelec;
        this->pelec = new elecstate::ElecStateLCAO<TK>(
            &(this->chr),
            &(this->kv),
            this->kv.get_nks(),
            &(this->LOC),
            &(this->GG),   // mohan add 2024-04-01
            &(this->GK),   // mohan add 2024-04-01
            &(this->LOWF), // should be replaced by a 2dbcd handle, if insist the "print_psi" must be in ElecState class
            this->pw_rho,
            this->pw_big);

        dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->init_DM(&this->kv, this->LM.ParaV, GlobalV::NSPIN);

        GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc, this->pw_rho);

        this->pelec->charge->allocate(GlobalV::NSPIN);
        this->pelec->omega = GlobalC::ucell.omega;

        // Initialize the potential.
        if (this->pelec->pot == nullptr)
        {
            this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                        this->pw_rho,
                                                        &GlobalC::ucell,
                                                        &(GlobalC::ppcell.vloc),
                                                        &(this->sf),
                                                        &(this->pelec->f_en.etxc),
                                                        &(this->pelec->f_en.vtxc));
        }
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "init_after_vc");
    return;
}

//------------------------------------------------------------------------------
//! the 5th function of ESolver_KS_LCAO: cal_energy
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
double ESolver_KS_LCAO<TK, TR>::cal_energy()
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "cal_energy");

    return this->pelec->f_en.etot;
}

//------------------------------------------------------------------------------
//! the 6th function of ESolver_KS_LCAO: cal_force
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::cal_force(ModuleBase::matrix& force)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "cal_force");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "cal_force");

	Force_Stress_LCAO<TK> fsl(this->RA, GlobalC::ucell.nat);

	fsl.getForceStress(
            GlobalV::CAL_FORCE,
			GlobalV::CAL_STRESS,
			GlobalV::TEST_FORCE,
			GlobalV::TEST_STRESS,
            this->orb_con.ParaV, 
			this->pelec,
			this->psi,
            this->LM,
            this->GG, // mohan add 2024-04-01
            this->GK, // mohan add 2024-04-01
            uot_,
			force,
			this->scs,
			this->sf,
			this->kv,
			this->pw_rho,
#ifdef __EXX
                       *this->exx_lri_double,
                       *this->exx_lri_complex,
#endif
                       &GlobalC::ucell.symm);

    // delete RA after cal_force

    this->RA.delete_grid();

    this->have_force = true;

    ModuleBase::timer::tick("ESolver_KS_LCAO", "cal_force");
}

//------------------------------------------------------------------------------
//! the 7th function of ESolver_KS_LCAO: cal_stress
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::cal_stress(ModuleBase::matrix& stress)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "cal_stress");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "cal_stress");

    if (!this->have_force)
    {
        ModuleBase::matrix fcs;
        this->cal_force(fcs);
    }
    stress = this->scs; // copy the stress
    this->have_force = false;

    ModuleBase::timer::tick("ESolver_KS_LCAO", "cal_stress");
}

//------------------------------------------------------------------------------
//! the 8th function of ESolver_KS_LCAO: after_all_runners
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::after_all_runners(void)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_all_runners");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_all_runners");

    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << this->pelec->f_en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

    if (INPUT.out_dos != 0 || INPUT.out_band[0] != 0 || INPUT.out_proj_band != 0)
    {
        GlobalV::ofs_running << "\n\n\n\n";
        GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        GlobalV::ofs_running << " |                                                                    |" << std::endl;
        GlobalV::ofs_running << " | Post-processing of data:                                           |" << std::endl;
        GlobalV::ofs_running << " | DOS (density of states) and bands will be output here.             |" << std::endl;
        GlobalV::ofs_running << " | If atomic orbitals are used, Mulliken charge analysis can be done. |" << std::endl;
        GlobalV::ofs_running << " | Also the .bxsf file containing fermi surface information can be    |" << std::endl;
        GlobalV::ofs_running << " | done here.                                                         |" << std::endl;
        GlobalV::ofs_running << " |                                                                    |" << std::endl;
        GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        GlobalV::ofs_running << "\n\n\n\n";
    }
    // qianrui modify 2020-10-18
    if (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "relax")
    {
        ModuleIO::write_istate_info(this->pelec->ekb, this->pelec->wg, this->kv, &(GlobalC::Pkpoints));
    }

    const int nspin0 = (GlobalV::NSPIN == 2) ? 2 : 1;

    if (INPUT.out_band[0])
    {
        for (int is = 0; is < nspin0; is++)
        {
            std::stringstream ss2;
            ss2 << GlobalV::global_out_dir << "BANDS_" << is + 1 << ".dat";
            GlobalV::ofs_running << "\n Output bands in file: " << ss2.str() << std::endl;
            ModuleIO::nscf_band(is,
                                ss2.str(),
                                GlobalV::NBANDS,
                                0.0,
                                INPUT.out_band[1],
                                this->pelec->ekb,
                                this->kv,
                                &(GlobalC::Pkpoints));
        }
    } // out_band

    if (INPUT.out_proj_band) // Projeced band structure added by jiyy-2022-4-20
    {
        ModuleIO::write_proj_band_lcao(this->psi, this->LM, this->pelec, this->kv, GlobalC::ucell, this->p_hamilt);
    }

    if (INPUT.out_dos)
    {
        ModuleIO::out_dos_nao(this->psi,
                              this->LM,
                              this->orb_con.ParaV,
                              this->pelec->ekb,
                              this->pelec->wg,
                              INPUT.dos_edelta_ev,
                              INPUT.dos_scale,
                              INPUT.dos_sigma,
                              *(this->pelec->klist),
                              GlobalC::Pkpoints,
                              GlobalC::ucell,
                              this->pelec->eferm,
                              GlobalV::NBANDS,
                              this->p_hamilt);
    }
    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_all_runners");
}

//------------------------------------------------------------------------------
//! the 9th function of ESolver_KS_LCAO: init_basis_lcao
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::init_basis_lcao(ORB_control& orb_con, Input& inp, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "init_basis_lcao");

    // autoset NB2D first
    if (GlobalV::NB2D == 0)
    {
        if (GlobalV::NLOCAL > 0)
        {
            GlobalV::NB2D = (GlobalV::NSPIN == 4) ? 2 : 1;
        }
        if (GlobalV::NLOCAL > 500)
        {
            GlobalV::NB2D = 32;
        }
        if (GlobalV::NLOCAL > 1000)
        {
            GlobalV::NB2D = 64;
        }
    }
    // Set the variables first
    this->orb_con.gamma_only = GlobalV::GAMMA_ONLY_LOCAL;
    this->orb_con.nlocal = GlobalV::NLOCAL;
    this->orb_con.nbands = GlobalV::NBANDS;
    this->orb_con.ParaV.nspin = GlobalV::NSPIN;
    this->orb_con.dsize = GlobalV::DSIZE;
    this->orb_con.nb2d = GlobalV::NB2D;
    this->orb_con.dcolor = GlobalV::DCOLOR;
    this->orb_con.drank = GlobalV::DRANK;
    this->orb_con.myrank = GlobalV::MY_RANK;
    this->orb_con.calculation = GlobalV::CALCULATION;
    this->orb_con.ks_solver = GlobalV::KS_SOLVER;
    this->orb_con.setup_2d = true;

    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.

    // NOTE: This following raw pointer serves as a temporary step in
    // LCAO refactoring. Eventually, it will be replaced by a shared_ptr,
    // which is the only owner of the ORB_gen_tables object. All other
    // usages will take a weak_ptr.
    uot_ = new ORB_gen_tables;
    auto& two_center_bundle = uot_->two_center_bundle;

    two_center_bundle.reset(new TwoCenterBundle);
    two_center_bundle->build_orb(ucell.ntype, ucell.orbital_fn);
    two_center_bundle->build_alpha(GlobalV::deepks_setorb, &ucell.descriptor_file);
    two_center_bundle->build_orb_onsite(ucell.ntype, GlobalV::onsite_radius);
    // currently deepks only use one descriptor file, so cast bool to int is fine

    // TODO Due to the omnipresence of GlobalC::ORB, we still have to rely
    // on the old interface for now.
    two_center_bundle->to_LCAO_Orbitals(GlobalC::ORB, inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax);

    ucell.infoNL.setupNonlocal(ucell.ntype, ucell.atoms, GlobalV::ofs_running, GlobalC::ORB);

    two_center_bundle->build_beta(ucell.ntype, ucell.infoNL.Beta);

    int Lmax = 0;
#ifdef __EXX
    Lmax = GlobalC::exx_info.info_ri.abfs_Lmax;
#endif

#ifndef USE_NEW_TWO_CENTER
    this->orb_con.set_orb_tables(GlobalV::ofs_running,
                                 *uot_,
                                 GlobalC::ORB,
                                 ucell.lat0,
                                 GlobalV::deepks_setorb,
                                 Lmax,
                                 ucell.infoNL.nprojmax,
                                 ucell.infoNL.nproj,
                                 ucell.infoNL.Beta);
#else
    two_center_bundle->tabulate();
#endif

    if (this->orb_con.setup_2d)
    {
        this->orb_con.setup_2d_division(GlobalV::ofs_running, GlobalV::ofs_warning);
        this->orb_con.ParaV.set_atomic_trace(GlobalC::ucell.get_iat2iwt(), GlobalC::ucell.nat, GlobalV::NLOCAL);
    }

    return;
}

//------------------------------------------------------------------------------
//! the 10th function of ESolver_KS_LCAO: iter_init
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::iter_init(const int istep, const int iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "iter_init");

    if (iter == 1)
    {
        this->p_chgmix->init_mixing(); // init mixing
        this->p_chgmix->mixing_restart_step = GlobalV::SCF_NMAX + 1;
        this->p_chgmix->mixing_restart_count = 0;
        // this output will be removed once the feeature is stable
        if (GlobalC::dftu.uramping > 0.01)
        {
            std::cout << " U-Ramping! Current U = ";
            for (int i = 0; i < GlobalC::dftu.U0.size(); i++)
            {
                std::cout << GlobalC::dftu.U[i] * ModuleBase::Ry_to_eV << " ";
            }
            std::cout << " eV " << std::endl;
        }
    }
    // for mixing restart
    if (iter == this->p_chgmix->mixing_restart_step && GlobalV::MIXING_RESTART > 0.0)
    {
        this->p_chgmix->init_mixing();
        this->p_chgmix->mixing_restart_count++;
        if (GlobalV::dft_plus_u)
        {
            GlobalC::dftu.uramping_update(); // update U by uramping if uramping > 0.01
            if (GlobalC::dftu.uramping > 0.01)
            {
                std::cout << " U-Ramping! Current U = ";
                for (int i = 0; i < GlobalC::dftu.U0.size(); i++)
                {
                    std::cout << GlobalC::dftu.U[i] * ModuleBase::Ry_to_eV << " ";
                }
                std::cout << " eV " << std::endl;
            }
            if (GlobalC::dftu.uramping > 0.01 && !GlobalC::dftu.u_converged())
            {
                this->p_chgmix->mixing_restart_step = GlobalV::SCF_NMAX + 1;
            }
        }
        if (GlobalV::MIXING_DMR) // for mixing_dmr
        {
            // allocate memory for dmr_mdata
            const elecstate::DensityMatrix<TK, double>* dm
                = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
            int nnr_tmp = dm->get_DMR_pointer(1)->get_nnr();
            this->p_chgmix->allocate_mixing_dmr(nnr_tmp);
        }
    }

    // mohan update 2012-06-05
    this->pelec->f_en.deband_harris = this->pelec->cal_delta_eband();

    // mohan move it outside 2011-01-13
    // first need to calculate the weight according to
    // electrons number.
    if (istep == 0 && this->wf.init_wfc == "file" // Note: on the removal of LOWF
        && this->LOWF.error == 0)                 // this means the wavefunction is read without any error.
    {                  // However the I/O of wavefunction are nonsence to be implmented in different places.
        if (iter == 1) // once the reading of wavefunction has any error, should exit immediately.
        {
            std::cout << " WAVEFUN -> CHARGE " << std::endl;

            // calculate the density matrix using read in wave functions
            // and the ncalculate the charge density on grid.

            this->pelec->skip_weights = true;
            this->pelec->psiToRho(*this->psi);
            this->pelec->skip_weights = false;

            // calculate the local potential(rho) again.
            // the grid integration will do in later grid integration.

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // a puzzle remains here.
            // if I don't renew potential,
            // The scf_thr is very small.
            // OneElectron, Hartree and
            // Exc energy are all correct
            // except the band energy.
            //
            // solved by mohan 2010-09-10
            // there are there rho here:
            // rho1: formed by read in orbitals.
            // rho2: atomic rho, used to construct H
            // rho3: generated by after diagonalize
            // here converged because rho3 and rho1
            // are very close.
            // so be careful here, make sure
            // rho1 and rho2 are the same rho.
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (GlobalV::NSPIN == 4)
            {
                GlobalC::ucell.cal_ux();
            }

            //! update the potentials by using new electron charge density
            this->pelec->pot->update_from_charge(this->pelec->charge, &GlobalC::ucell);

            //! compute the correction energy for metals
            this->pelec->f_en.descf = this->pelec->cal_delta_escf();
        }
    }

#ifdef __EXX
    // calculate exact-exchange
    if (GlobalC::exx_info.info_ri.real_number)
    {
        this->exd->exx_eachiterinit(*dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(), iter);
    }
    else
    {
        this->exc->exx_eachiterinit(*dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(), iter);
    }
#endif

    if (GlobalV::dft_plus_u)
    {
        if (istep != 0 || iter != 1)
        {
            GlobalC::dftu.set_dmr(dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM());
        }
        // Calculate U and J if Yukawa potential is used
        GlobalC::dftu.cal_slater_UJ(this->pelec->charge->rho, this->pw_rho->nrxx);
    }

#ifdef __DEEPKS
    // the density matrixes of DeePKS have been updated in each iter
    GlobalC::ld.set_hr_cal(true);

    // HR in HamiltLCAO should be recalculate
    if (GlobalV::deepks_scf)
    {
        this->p_hamilt->refresh();
    }
#endif

    if (GlobalV::VL_IN_H)
    {
        // update Gint_K
        if (!GlobalV::GAMMA_ONLY_LOCAL)
        {
            this->GK.renew();
        }
        // update real space Hamiltonian
        this->p_hamilt->refresh();
    }

    // run the inner lambda loop to contrain atomic moments with the DeltaSpin method
    if (GlobalV::sc_mag_switch && iter > GlobalV::sc_scf_nmin)
    {
        SpinConstrain<TK, base_device::DEVICE_CPU>& sc = SpinConstrain<TK, base_device::DEVICE_CPU>::getScInstance();
        sc.run_lambda_loop(iter - 1);
    }
}

//------------------------------------------------------------------------------
//! the 11th function of ESolver_KS_LCAO: hamilt2density
//! mohan add 2024-05-11
//! 1) save input rho
//! 2) save density matrix DMR for mixing
//! 3) solve the Hamiltonian and output band gap
//! 4) print bands for each k-point and each band
//! 5) EXX:
//! 6) DFT+U: compute local occupation number matrix and energy correction
//! 7) DeePKS: compute delta_e
//! 8) DeltaSpin:
//! 9) use new charge density to calculate energy
//! 10) symmetrize the charge density
//! 11) compute magnetization, only for spin==2
//! 12) calculate delta energy
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::hamilt2density(int istep, int iter, double ethr)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "hamilt2density");

    // 1) save input rho
    this->pelec->charge->save_rho_before_sum_band();

    // 2) save density matrix DMR for mixing
    if (GlobalV::MIXING_RESTART > 0 && GlobalV::MIXING_DMR && this->p_chgmix->mixing_restart_count > 0)
    {
        elecstate::DensityMatrix<TK, double>* dm = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        dm->save_DMR();
    }

    // 3) solve the Hamiltonian and output band gap
    if (this->phsol != nullptr)
    {
        // reset energy
        this->pelec->f_en.eband = 0.0;
        this->pelec->f_en.demet = 0.0;

        this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec, GlobalV::KS_SOLVER);

        if (GlobalV::out_bandgap)
        {
            if (!GlobalV::TWO_EFERMI)
            {
                this->pelec->cal_bandgap();
            }
            else
            {
                this->pelec->cal_bandgap_updw();
            }
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW", "HSolver has not been initialed!");
    }

    // 4) print bands for each k-point and each band
    for (int ik = 0; ik < this->kv.get_nks(); ++ik)
    {
        this->pelec->print_band(ik, INPUT.printe, iter);
    }

    // 5) what's the exd used for?
#ifdef __EXX
    if (GlobalC::exx_info.info_ri.real_number)
    {
        this->exd->exx_hamilt2density(*this->pelec, this->orb_con.ParaV, iter);
    }
    else
    {
        this->exc->exx_hamilt2density(*this->pelec, this->orb_con.ParaV, iter);
    }
#endif

    // 6) calculate the local occupation number matrix and energy correction in DFT+U
    if (GlobalV::dft_plus_u)
    {
        // only old DFT+U method should calculated energy correction in esolver,
        // new DFT+U method will calculate energy in calculating Hamiltonian
        if (GlobalV::dft_plus_u == 2)
        {
            if (GlobalC::dftu.omc != 2)
            {
                const std::vector<std::vector<TK>>& tmp_dm
                    = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMK_vector();
                this->dftu_cal_occup_m(iter, tmp_dm);
            }
            GlobalC::dftu.cal_energy_correction(istep);
        }
        GlobalC::dftu.output();
    }

    // (7) for deepks, calculate delta_e
#ifdef __DEEPKS
    if (GlobalV::deepks_scf)
    {
        const std::vector<std::vector<TK>>& dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMK_vector();

        this->dpks_cal_e_delta_band(dm);
    }
#endif

    // 8) for delta spin
    if (GlobalV::sc_mag_switch)
    {
        SpinConstrain<TK, base_device::DEVICE_CPU>& sc = SpinConstrain<TK, base_device::DEVICE_CPU>::getScInstance();
        sc.cal_MW(iter, &(this->LM));
    }

    // 9) use new charge density to calculate energy
    this->pelec->cal_energies(1);

    // 10) symmetrize the charge density
    Symmetry_rho srho;
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        srho.begin(is, *(this->pelec->charge), this->pw_rho, GlobalC::Pgrid, GlobalC::ucell.symm);
    }

    // 11) compute magnetization, only for spin==2
    GlobalC::ucell.magnet.compute_magnetization(this->pelec->charge->nrxx,
                                                this->pelec->charge->nxyz,
                                                this->pelec->charge->rho,
                                                this->pelec->nelec_spin.data());

    // 12) calculate delta energy
    this->pelec->f_en.deband = this->pelec->cal_delta_eband();
}

//------------------------------------------------------------------------------
//! the 12th function of ESolver_KS_LCAO: update_pot
//! mohan add 2024-05-11
//! 1) print Hamiltonian and Overlap matrix (why related to update_pot()?)
//! 2) print wavefunctions (why related to update_pot()?)
//! 3) print potential
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::update_pot(const int istep, const int iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "update_pot");

    // 1) print Hamiltonian and Overlap matrix
    if (this->conv_elec || iter == GlobalV::SCF_NMAX)
    {
        if (!GlobalV::GAMMA_ONLY_LOCAL && hsolver::HSolverLCAO<TK>::out_mat_hs[0])
        {
            this->GK.renew(true);
        }
        for (int ik = 0; ik < this->kv.get_nks(); ++ik)
        {
            if (hsolver::HSolverLCAO<TK>::out_mat_hs[0])
            {
                this->p_hamilt->updateHk(ik);
            }
            bool bit = false; // LiuXh, 2017-03-21
            // if set bit = true, there would be error in soc-multi-core calculation, noted by zhengdy-soc
            if (this->psi != nullptr && (istep % GlobalV::out_interval == 0))
            {
                hamilt::MatrixBlock<TK> h_mat, s_mat;
                this->p_hamilt->matrix(h_mat, s_mat);
                if (hsolver::HSolverLCAO<TK>::out_mat_hs[0])
                {
                    ModuleIO::save_mat(istep,
                                       h_mat.p,
                                       GlobalV::NLOCAL,
                                       bit,
                                       hsolver::HSolverLCAO<TK>::out_mat_hs[1],
                                       1,
                                       GlobalV::out_app_flag,
                                       "H",
                                       "data-" + std::to_string(ik),
                                       this->orb_con.ParaV,
                                       GlobalV::DRANK);
                    ModuleIO::save_mat(istep,
                                       s_mat.p,
                                       GlobalV::NLOCAL,
                                       bit,
                                       hsolver::HSolverLCAO<TK>::out_mat_hs[1],
                                       1,
                                       GlobalV::out_app_flag,
                                       "S",
                                       "data-" + std::to_string(ik),
                                       this->orb_con.ParaV,
                                       GlobalV::DRANK);
                }
            }
        }
    }

    // 2) print wavefunctions
    if (elecstate::ElecStateLCAO<TK>::out_wfc_lcao && (this->conv_elec || iter == GlobalV::SCF_NMAX)
        && (istep % GlobalV::out_interval == 0))
    {
        ModuleIO::write_wfc_nao(elecstate::ElecStateLCAO<TK>::out_wfc_lcao,
                                this->psi[0],
                                this->pelec->ekb,
                                this->pelec->wg,
                                this->pelec->klist->kvec_c,
                                this->orb_con.ParaV,
                                istep);
    }

    // 3) print potential
    if (this->conv_elec || iter == GlobalV::SCF_NMAX)
    {
        if (GlobalV::out_pot < 0) // mohan add 2011-10-10
        {
            GlobalV::out_pot = -2;
        }
    }

    if (!this->conv_elec)
    {
        if (GlobalV::NSPIN == 4)
        {
            GlobalC::ucell.cal_ux();
        }
        this->pelec->pot->update_from_charge(this->pelec->charge, &GlobalC::ucell);
        this->pelec->f_en.descf = this->pelec->cal_delta_escf();
    }
    else
    {
        this->pelec->cal_converged();
    }
}

//------------------------------------------------------------------------------
//! the 13th function of ESolver_KS_LCAO: iter_finish
//! mohan add 2024-05-11
//! 1) mix density matrix
//! 2) output charge density
//! 3) output exx matrix
//! 4) output charge density and density matrix
//! 5) cal_MW? (why put it here?)
//! 6) calculate the total energy?
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::iter_finish(int iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "iter_finish");

    // 1) mix density matrix if mixing_restart + mixing_dmr + not first mixing_restart at every iter
    if (GlobalV::MIXING_RESTART > 0 && this->p_chgmix->mixing_restart_count > 0 && GlobalV::MIXING_DMR)
    {
        elecstate::DensityMatrix<TK, double>* dm = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        this->p_chgmix->mix_dmr(dm);
    }

    // 2) save charge density
    // Peize Lin add 2020.04.04
    if (GlobalC::restart.info_save.save_charge)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            GlobalC::restart.save_disk("charge", is, this->pelec->charge->nrxx, this->pelec->charge->rho[is]);
        }
    }

#ifdef __EXX
    // 3) save exx matrix
    int two_level_step = GlobalC::exx_info.info_ri.real_number ? this->exd->two_level_step : this->exc->two_level_step;

    if (GlobalC::restart.info_save.save_H && two_level_step > 0
        && (!GlobalC::exx_info.info_global.separate_loop || iter == 1)) // to avoid saving the same value repeatedly
    {
        std::vector<TK> Hexxk_save(this->orb_con.ParaV.get_local_size());
        for (int ik = 0; ik < this->kv.get_nks(); ++ik)
        {
            ModuleBase::GlobalFunc::ZEROS(Hexxk_save.data(), Hexxk_save.size());

            hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>> opexx_save(&this->LM, nullptr, &Hexxk_save, this->kv);

            opexx_save.contributeHk(ik);

            GlobalC::restart.save_disk("Hexx", ik, this->orb_con.ParaV.get_local_size(), Hexxk_save.data());
        }
        if (GlobalV::MY_RANK == 0)
        {
            GlobalC::restart.save_disk("Eexx", 0, 1, &this->pelec->f_en.exx);
        }
    }
#endif

    // 4) output charge density and density matrix
    bool print = false;
    if (this->out_freq_elec && iter % this->out_freq_elec == 0)
    {
        print = true;
    }

    if (print)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            this->create_Output_Rho(is, iter, "tmp_").write();
            this->create_Output_DM(is, iter).write();
            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                this->create_Output_Kin(is, iter, "tmp_").write();
            }
        }
    }

    // 5) cal_MW?
    // escon: energy of spin constraint depends on Mi, so cal_energies should be called after cal_MW
    if (GlobalV::sc_mag_switch)
    {
        SpinConstrain<TK, base_device::DEVICE_CPU>& sc = SpinConstrain<TK, base_device::DEVICE_CPU>::getScInstance();
        sc.cal_MW(iter, &(this->LM));
    }

    // 6) calculate the total energy.
    this->pelec->cal_energies(2);
}

//------------------------------------------------------------------------------
//! the 14th function of ESolver_KS_LCAO: after_scf
//! mohan add 2024-05-11
//! 1) write charge difference into files for charge extrapolation
//! 2) write density matrix for sparse matrix
//! 3) write charge density
//! 4) write density matrix
//! 5) write Vxc
//! 6) write Exx matrix
//! 7) write potential
//! 8) write convergence
//! 9) write fermi energy
//! 10) write eigenvalues
//! 11) write deepks information
//! 12) write rpa information
//! 13) write HR in npz format
//! 14) write dm in npz format
//! 15) write md related
//! 16) write spin constrian MW?
//! 17) delete grid
//! 18) write quasi-orbitals
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::after_scf(const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_scf");

    // 1) write charge difference into files for charge extrapolation
    if (GlobalV::CALCULATION != "scf")
    {
        this->CE.save_files(istep,
                            GlobalC::ucell,
#ifdef __MPI
                            this->pw_big,
#endif
                            this->pelec->charge,
                            &this->sf);
    }

    // 2) write density matrix for sparse matrix
    if (this->LOC.out_dm1 == 1)
    {
        this->create_Output_DM1(istep).write();
    }

    // 3) write charge density
    if (GlobalV::out_chg)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            this->create_Output_Rho(is, istep).write();
            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                this->create_Output_Kin(is, istep).write();
            }
        }
    }

    // 4) write density matrix
    if (this->LOC.out_dm)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            this->create_Output_DM(is, istep).write();
        }
    }

    // 5) write Vxc
    bool out_exc = true; // tmp, add parameter!
    if (GlobalV::out_mat_xc)
    {
        ModuleIO::write_Vxc<TK, TR>(GlobalV::NSPIN,
                                    GlobalV::NLOCAL,
                                    GlobalV::DRANK,
                                    *this->psi,
                                    GlobalC::ucell,
                                    this->sf,
                                    *this->pw_rho,
                                    *this->pw_rhod,
                                    GlobalC::ppcell.vloc,
                                    *this->pelec->charge,
                                    this->GG,
                                    this->GK,
                                    this->LM,
                                    this->LOC,
                                    this->kv,
                                    this->pelec->wg,
                                    GlobalC::GridD);
    }

#ifdef __EXX
    // 6) write Exx matrix
    if (GlobalC::exx_info.info_global.cal_exx) // Peize Lin add if 2022.11.14
    {
        const std::string file_name_exx = GlobalV::global_out_dir + "HexxR" + std::to_string(GlobalV::MY_RANK);
        if (GlobalC::exx_info.info_ri.real_number)
        {
            this->exd->write_Hexxs_csr(file_name_exx, GlobalC::ucell);
        }
        else
        {
            this->exc->write_Hexxs_csr(file_name_exx, GlobalC::ucell);
        }
    }
#endif

    // 7) write potential
    this->create_Output_Potential(istep).write();

    // 8) write convergence
    ModuleIO::output_convergence_after_scf(this->conv_elec, this->pelec->f_en.etot);

    // 9) write fermi energy
    ModuleIO::output_efermi(this->conv_elec, this->pelec->eferm.ef);

    // 10) write eigenvalues
    if (GlobalV::OUT_LEVEL != "m")
    {
        this->pelec->print_eigenvalue(GlobalV::ofs_running);
    }

    // 11) write deepks information
#ifdef __DEEPKS
    std::shared_ptr<LCAO_Deepks> ld_shared_ptr(&GlobalC::ld, [](LCAO_Deepks*) {});
    LCAO_Deepks_Interface LDI = LCAO_Deepks_Interface(ld_shared_ptr);
    ModuleBase::timer::tick("ESolver_KS_LCAO", "out_deepks_labels");
    LDI.out_deepks_labels(this->pelec->f_en.etot,
                          this->pelec->klist->get_nks(),
                          GlobalC::ucell.nat,
                          this->pelec->ekb,
                          this->pelec->klist->kvec_d,
                          GlobalC::ucell,
                          GlobalC::ORB,
                          GlobalC::GridD,
                          &(this->orb_con.ParaV),
                          *(this->psi),
                          dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM());

    ModuleBase::timer::tick("ESolver_KS_LCAO", "out_deepks_labels");
#endif

#ifdef __EXX
    // 12) write rpa information
    if (INPUT.rpa)
    {
        // ModuleRPA::DFT_RPA_interface rpa_interface(GlobalC::exx_info.info_global);
        // rpa_interface.rpa_exx_lcao().info.files_abfs = GlobalV::rpa_orbitals;
        // rpa_interface.out_for_RPA(*(this->LOWF.ParaV), *(this->psi), this->LOC, this->pelec);
        RPA_LRI<TK, double> rpa_lri_double(GlobalC::exx_info.info_ri);
        rpa_lri_double.cal_postSCF_exx(*dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                                       MPI_COMM_WORLD,
                                       this->kv);
        rpa_lri_double.init(MPI_COMM_WORLD, this->kv);
        rpa_lri_double.out_for_RPA(this->orb_con.ParaV, *(this->psi), this->pelec);
    }
#endif

    // 13) write HR in npz format
    if (GlobalV::out_hr_npz)
    {
        this->p_hamilt->updateHk(0); // first k point, up spin
        hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao
            = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(this->p_hamilt);
        std::string zipname = "output_HR0.npz";
        this->output_mat_npz(zipname, *(p_ham_lcao->getHR()));

        if (GlobalV::NSPIN == 2)
        {
            this->p_hamilt->updateHk(this->kv.get_nks() / 2); // the other half of k points, down spin
            hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao
                = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(this->p_hamilt);
            zipname = "output_HR1.npz";
            this->output_mat_npz(zipname, *(p_ham_lcao->getHR()));
        }
    }

    // 14) write dm in npz format
    if (GlobalV::out_dm_npz)
    {
        const elecstate::DensityMatrix<TK, double>* dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        std::string zipname = "output_DM0.npz";
        this->output_mat_npz(zipname, *(dm->get_DMR_pointer(1)));

        if (GlobalV::NSPIN == 2)
        {
            zipname = "output_DM1.npz";
            this->output_mat_npz(zipname, *(dm->get_DMR_pointer(2)));
        }
    }

    // 15) write md related
    if (!md_skip_out(GlobalV::CALCULATION, istep, GlobalV::out_interval))
    {
        this->create_Output_Mat_Sparse(istep).write();
        // mulliken charge analysis
        if (GlobalV::out_mul)
        {
            this->cal_mag(istep, true);
        }
    }

    // 16) write spin constrian MW?
    // spin constrain calculations, added by Tianqi Zhao.
    if (GlobalV::sc_mag_switch)
    {
        SpinConstrain<TK, base_device::DEVICE_CPU>& sc = SpinConstrain<TK, base_device::DEVICE_CPU>::getScInstance();
        sc.cal_MW(istep, &(this->LM), true);
        sc.print_Mag_Force();
    }

    // 17) delete grid
    if (!GlobalV::CAL_FORCE && !GlobalV::CAL_STRESS)
    {
        RA.delete_grid();
    }

    // 18) write quasi-orbitals, added by Yike Huang.
    if (GlobalV::qo_switch)
    {
        toQO tqo(GlobalV::qo_basis, GlobalV::qo_strategy, GlobalV::qo_thr, GlobalV::qo_screening_coeff);
        tqo.initialize(GlobalV::global_out_dir,
                       GlobalV::global_pseudo_dir,
                       GlobalV::global_orbital_dir,
                       &GlobalC::ucell,
                       this->kv.kvec_d,
                       GlobalV::ofs_running,
                       GlobalV::MY_RANK,
                       GlobalV::NPROC);
        tqo.calculate();
    }
}

//------------------------------------------------------------------------------
//! the 15th function of ESolver_KS_LCAO: do_after_converge
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
bool ESolver_KS_LCAO<TK, TR>::do_after_converge(int& iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "do_after_converge");

#ifdef __EXX
    if (GlobalC::exx_info.info_ri.real_number)
    {
        return this->exd->exx_after_converge(*this->p_hamilt,
                                             this->LM,
                                             *dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                                             this->kv,
                                             iter);
    }
    else
    {
        return this->exc->exx_after_converge(*this->p_hamilt,
                                             this->LM,
                                             *dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                                             this->kv,
                                             iter);
    }
#endif // __EXX

    if (GlobalV::dft_plus_u)
    {
        // use the converged occupation matrix for next MD/Relax SCF calculation
        GlobalC::dftu.initialed_locale = true;
    }

    return true;
}

//------------------------------------------------------------------------------
//! the 16th function of ESolver_KS_LCAO: create_Output_DM
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
ModuleIO::Output_DM ESolver_KS_LCAO<TK, TR>::create_Output_DM(int is, int iter)
{
    const int precision = 3;

    return ModuleIO::Output_DM(this->GridT,
                               is,
                               iter,
                               precision,
                               this->LOC.out_dm,
                               this->LOC.DM,
                               this->pelec->eferm.get_efval(is),
                               &(GlobalC::ucell),
                               GlobalV::global_out_dir,
                               GlobalV::GAMMA_ONLY_LOCAL);
}

//------------------------------------------------------------------------------
//! the 17th function of ESolver_KS_LCAO: create_Output_DM1
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
ModuleIO::Output_DM1 ESolver_KS_LCAO<TK, TR>::create_Output_DM1(int istep)
{
    const elecstate::DensityMatrix<complex<double>, double>* DM
        = dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)->get_DM();
    return ModuleIO::Output_DM1(GlobalV::NSPIN, istep, this->LOC, this->RA, this->kv, DM);
}

//------------------------------------------------------------------------------
//! the 18th function of ESolver_KS_LCAO: create_Output_Mat_Sparse
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
ModuleIO::Output_Mat_Sparse<TK> ESolver_KS_LCAO<TK, TR>::create_Output_Mat_Sparse(int istep)
{
	return ModuleIO::Output_Mat_Sparse<TK>(
            hsolver::HSolverLCAO<TK>::out_mat_hsR,
			hsolver::HSolverLCAO<TK>::out_mat_dh,
			hsolver::HSolverLCAO<TK>::out_mat_t,
			INPUT.out_mat_r,
			istep,
			this->pelec->pot->get_effective_v(),
			this->orb_con.ParaV,
            this->GK, // mohan add 2024-04-01
            uot_,
			this->LM,
            GlobalC::GridD, // mohan add 2024-04-06
			this->kv,
			this->p_hamilt);
}

//------------------------------------------------------------------------------
//! the 19th function of ESolver_KS_LCAO: md_skip_out
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
bool ESolver_KS_LCAO<TK, TR>::md_skip_out(std::string calculation, int istep, int interval)
{
    if (calculation == "md")
    {
        if (istep % interval != 0)
        {
            return true;
        }
    }
    return false;
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::cal_mag(const int istep, const bool print)
{
    auto cell_index = CellIndex(GlobalC::ucell.get_atomLabels(),
                                GlobalC::ucell.get_atomCounts(),
                                GlobalC::ucell.get_lnchiCounts(),
                                GlobalV::NSPIN);
    auto out_sk = ModuleIO::Output_Sk<TK>(&(this->LM),
                                          this->p_hamilt,
                                          &(this->orb_con.ParaV),
                                          GlobalV::NSPIN,
                                          this->kv.get_nks());
    auto out_dmk = ModuleIO::Output_DMK<TK>(dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                                            &(this->orb_con.ParaV),
                                            GlobalV::NSPIN,
                                            this->kv.get_nks());
    auto mulp = ModuleIO::Output_Mulliken<TK>(&(out_sk),
                                              &(out_dmk),
                                              &(this->orb_con.ParaV),
                                              &cell_index,
                                              this->kv.isk,
                                              GlobalV::NSPIN);
    auto atom_chg = mulp.get_atom_chg();
    /// used in updating mag info in STRU file
    GlobalC::ucell.atom_mulliken = mulp.get_atom_mulliken(atom_chg);
    if (print && GlobalV::MY_RANK == 0)
    {
        /// write the Orbital file
        cell_index.write_orb_info(GlobalV::global_out_dir);
        /// write mulliken.txt
        mulp.write(istep, GlobalV::global_out_dir);
        /// write atomic mag info in running log file
        mulp.print_atom_mag(atom_chg, GlobalV::ofs_running);
    }
}

//------------------------------------------------------------------------------
//! the 20th,21th,22th functions of ESolver_KS_LCAO
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
