#include "driver.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_esolver/esolver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/cal_test.h"
#include "module_io/input.h"
#include "module_io/input_conv.h"
#include "module_io/para_json.h"
#include "module_io/print_info.h"
#include "module_io/winput.h"
Driver::Driver()
{
}

Driver::~Driver()
{
    // Release the device memory within singleton object GlobalC::ppcell
    // before the main function exits.
    GlobalC::ppcell.release_memory();
}

void Driver::init()
{
    ModuleBase::TITLE("Driver", "init");

    time_t time_start = std::time(nullptr);
    ModuleBase::timer::start();

    // (1) read the input parameters.
    // INPUT should be initalized here and then pass to atomic world, mohan 2024-05-12
    // INPUT should not be GlobalC, mohan 2024-05-12
    Driver::reading();

    // (2) welcome to the atomic world!
    this->atomic_world();

    // (3) output information
    time_t time_finish = std::time(nullptr);
    Print_Info::print_time(time_start, time_finish);

    // (4) close all of the running logs
    INPUT.close_log();

    // (5) output the json file
    // Json::create_Json(&GlobalC::ucell.symm,GlobalC::ucell.atoms,&INPUT);
    Json::create_Json(&GlobalC::ucell, &INPUT);
}

void Driver::reading()
{
    ModuleBase::timer::tick("Driver", "reading");

    // (1) read INPUT
    INPUT.Init(GlobalV::global_in_card);

    // (2) copy the variables from INPUT to each class
    Input_Conv::Convert();

    // (3) define the 'DIAGONALIZATION' world in MPI
    Parallel_Global::split_diag_world(GlobalV::DIAGO_PROC);
    Parallel_Global::split_grid_world(GlobalV::DIAGO_PROC);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "DRANK", GlobalV::DRANK + 1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "DSIZE", GlobalV::DSIZE);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "DCOLOR", GlobalV::DCOLOR + 1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "GRANK", GlobalV::GRANK + 1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "GSIZE", GlobalV::GSIZE);

#ifdef __MPI
    // (4)  divide the GlobalV::NPROC processors into GlobalV::KPAR for k-points parallelization.
    Parallel_Global::init_pools();
#endif

    // (5) Read in parameters about wannier functions.
    winput::Init(GlobalV::global_wannier_card);

    // (6) Print the parameters into INPUT file.
    std::stringstream ss1;
    ss1 << GlobalV::global_out_dir << GlobalV::global_in_card;
    INPUT.Print(ss1.str());
    // ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"READING CARDS");

    ModuleBase::timer::tick("Driver", "reading");
}

void Driver::atomic_world()
{
    ModuleBase::TITLE("Driver", "atomic_world");
    //--------------------------------------------------
    // choose basis sets:
    // pw: plane wave basis set
    // lcao_in_pw: LCAO expaned by plane wave basis set
    // lcao: linear combination of atomic orbitals
    //--------------------------------------------------

    // where the actual stuff is done
    this->driver_run();

    ModuleBase::timer::finish(GlobalV::ofs_running);
    ModuleBase::Memory::print_all(GlobalV::ofs_running);
}
