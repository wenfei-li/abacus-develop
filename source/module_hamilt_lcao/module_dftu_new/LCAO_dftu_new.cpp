//wenfei 2022-1-5
//This file contains constructor and destructor of the class LCAO_dftu_new, 
//as well as subroutines for initializing and releasing relevant data structures 

//Other than the constructor and the destructor, it contains 3 types of subroutines:
//1. subroutines that are related to calculating onsite dm:
//  - init : allocates some arrays
//  - init_index : records the index (inl)
//  - allocate_nlm : allocates data structures (nlm_save) which is used to store <chi|alpha>
//4. subroutines that are related to V_delta:
//  - allocate_V_delta : allocates H_V_delta; if calculating force, it also calls
//      init_gdmx, as well as allocating F_delta
//  - allocate_V_deltaR : allcoates H_V_deltaR, for multi-k calculations

#include "LCAO_dftu_new.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/parallel_common.h"

namespace GlobalC
{
    LCAO_DftU_New lcao_dftu_new;
}

//Constructor of the class
LCAO_DftU_New::LCAO_DftU_New()
{
    inl_index = new ModuleBase::IntArray[1];
    inl_l = nullptr;

    allocated_H = false;
    allocated_HR = false;
}

//Desctructor of the class
LCAO_DftU_New::~LCAO_DftU_New()
{
    delete[] inl_index;
    delete[] inl_l;
    if(allocated_H)
    {
        if(GlobalV::GAMMA_ONLY_LOCAL)
        {
            for(int is = 0; is < GlobalV::NSPIN; is ++)
            {
                delete[] H_V_delta[is];
            }
            delete[] H_V_delta;
        }
        else
        {
            for(int ik=0;ik<nks;ik++)
            {
                delete[] this->H_V_delta_k[ik];
            }
            delete[] H_V_delta_k;
        }
    }

    if(allocated_HR)
    {
        for(int is = 0; is < GlobalV::NSPIN; is ++)
        {
            delete[] H_V_deltaR[is];
        }
        delete[] H_V_deltaR;
    }

    //=======1. to use deepks, pdm is required==========
    //delete pdm**
    if(allocated_pdm)
    {
        for(int is = 0; is < GlobalV::NSPIN; is ++)
        {
            for (int inl = 0;inl < this->inlmax;inl++)
            {
                delete[]  pdm[is][inl];
                delete[] gedm[is][inl];
            }
            delete[]  pdm[is];
            delete[] gedm[is];
        }
        delete[]  pdm;
        delete[] gedm;
    }
}

void LCAO_DftU_New::init(
    const LCAO_Orbitals& orb,
    const int nat,
    const int ntype,
    const Parallel_Orbitals& pv_in,
    std::vector<int> na,
    std::vector<std::vector<double>> & uvalue_in)
{
    ModuleBase::TITLE("LCAO_DftU_New", "init");

    // I am not so fond of the way this part is done, but it suffices for now
    
    GlobalV::ofs_running << " Initialize the projector index for DeePKS (lcao line)" << std::endl;

    const int lm = orb.get_lmax_d();
    const int nm = orb.get_nchimax_d();
    const int tot_inl_per_atom = orb.Alpha[0].getTotal_nchi();

    assert(lm >= 0);
    assert(nm >= 0);
    assert(tot_inl_per_atom >= 0);
    
    const int tot_inl = tot_inl_per_atom * nat;

    this->lmaxd = lm;
    this->nmaxd = nm;
    this->inlmax = tot_inl;
    GlobalV::ofs_running << " lmax of projector = " << this->lmaxd << std::endl;
    GlobalV::ofs_running << " nmax of projector = " << nmaxd << std::endl;
    
    //init pdm ***
    if(!allocated_pdm)
    {
        const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
        this->pdm = new double** [GlobalV::NSPIN];
        for(int is = 0; is < GlobalV::NSPIN; is ++)
        {
            this->pdm[is] = new double* [this->inlmax];
            for (int inl = 0;inl < this->inlmax;inl++)
            {
                this->pdm[is][inl] = new double[pdm_size];
                ModuleBase::GlobalFunc::ZEROS(this->pdm[is][inl], pdm_size);
            }
        }
        allocated_pdm = true;
    }

    this->init_index(ntype, nat, na, tot_inl, orb);
    this->allocate_nlm(nat);

    this->pv = &pv_in;

    if_has_u.resize(nat);
    uvalue.resize(nat);
    for(int iat = 0; iat < nat ; iat ++)
    {
        uvalue[iat].resize(lm+1);
        ModuleBase::GlobalFunc::ZEROS(uvalue[iat].data(), lm+1);
    }

    int iat = 0;
    for(int it = 0; it < ntype; it ++)
    {
        for(int ia = 0; ia < na[it]; ia ++)
        {
            bool u_applied = false;
            for(int il = 0; il < uvalue_in[it].size(); il ++)
            {
                if(std::abs(uvalue_in[it][il]) > 1e-8)
                {
                    u_applied = true;
                    uvalue[iat][il] = uvalue_in[it][il];
                }
            }
            if_has_u[iat] = u_applied;
            iat ++;
        }
    }

    return;
}

void LCAO_DftU_New::init_index(const int ntype, const int nat, std::vector<int> na, const int Total_nchi, const LCAO_Orbitals &orb)
{
    delete[] this->inl_index;
    this->inl_index = new ModuleBase::IntArray[ntype];
    delete[] this->inl_l;
    this->inl_l = new int[this->inlmax];
    ModuleBase::GlobalFunc::ZEROS(this->inl_l, this->inlmax);

    int inl = 0;
    int alpha = 0;
    for (int it = 0; it < ntype; it++)
    {
        this->inl_index[it].create(
            na[it],
            this->lmaxd + 1,
            this->nmaxd); 

        GlobalV::ofs_running << " Type " << it + 1
                    << " number_of_atoms " << na[it] << std::endl;

        for (int ia = 0; ia < na[it]; ia++)
        {
            for (int l = 0; l < this->lmaxd + 1; l++)
            {
                for (int n = 0; n < orb.Alpha[0].getNchi(l); n++)
                {
                    this->inl_index[it](ia, l, n) = inl;
                    this->inl_l[inl] = l;
                    inl++;
                }
            }
        }//end ia
    }//end it
    assert(Total_nchi == inl);

	return;
}

void LCAO_DftU_New::allocate_nlm(const int nat)
{
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        this->nlm_save.resize(nat);
    }
    else
    {
        this->nlm_save_k.resize(nat);
    }
}

void LCAO_DftU_New::allocate_V_delta(const int nat, const int nks_in)
{
    ModuleBase::TITLE("LCAO_DftU_New", "allocate_V_delta");

    nks = nks_in;

    //initialize the H matrix H_V_delta
    if(!allocated_H)
    {
        if(GlobalV::GAMMA_ONLY_LOCAL)
        {
            this->H_V_delta = new double* [GlobalV::NSPIN];
            for(int is = 0; is < GlobalV::NSPIN; is ++)
            {
                H_V_delta[is] = new double[pv->nloc];
                ModuleBase::GlobalFunc::ZEROS(this->H_V_delta[is], pv->nloc);
            }
        }
        else
        {
            H_V_delta_k = new std::complex<double>* [nks];
            for(int ik=0;ik<nks;ik++)
            {
                this->H_V_delta_k[ik] = new std::complex<double>[pv->nloc];
                ModuleBase::GlobalFunc::ZEROS(this->H_V_delta_k[ik], pv->nloc);
            }
        }
        allocated_H = true;
    }

    //init gedm**
    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    this->gedm = new double** [GlobalV::NSPIN];
    for(int is = 0; is < GlobalV::NSPIN; is ++)
    {
        this->gedm[is] = new double* [this->inlmax];
        for (int inl = 0;inl < this->inlmax;inl++)
        {
            this->gedm[is][inl] = new double[pdm_size];
            ModuleBase::GlobalFunc::ZEROS(this->gedm[is][inl], pdm_size);
        }
    }
    if (GlobalV::CAL_FORCE)
    {
        //init F_delta
        F_delta.create(nat, 3);
    }

    return;
}

void LCAO_DftU_New::allocate_V_deltaR(const int nnr)
{
    ModuleBase::TITLE("LCAO_DftU_New", "allocate_V_deltaR");
    if(!allocated_HR)
    {
        H_V_deltaR = new double* [GlobalV::NSPIN];
        for(int is = 0; is < GlobalV::NSPIN; is ++)
        {
            H_V_deltaR[is] = new double[nnr];
            ModuleBase::GlobalFunc::ZEROS(H_V_deltaR[is], nnr);
        }
        allocated_HR = true;
    }
}

void LCAO_DftU_New::read_info(std::vector<std::vector<double>> & uvalue_in, const int ntype)
{
    ModuleBase::TITLE("LCAO_DftU_New","read_info");

    int lmax;
    if(GlobalV::MY_RANK == 0)
    {
        std::ifstream ifs("INPUT");
        if(ifs.fail())
        {
            ModuleBase::WARNING_QUIT("lcao_dftu_new.cpp","input file not found!");
        }

        std::string line;
        line = this->scan_file(ifs,"dftu_lmax");
        
        std::string key;
        std::stringstream tmp;

        tmp << line;
        tmp >> key >> lmax;

        uvalue_in.resize(ntype);
        for(int it = 0; it < ntype; it ++)
        {
            uvalue_in[it].resize(lmax+1);
        }

        ifs.clear();
        ifs.seekg (0, std::ios::beg);    

        line = this->scan_file(ifs,"dftu_uvalue");
        tmp.clear();
        tmp << line;
        tmp >> key;

        for(int it = 0; it < ntype; it ++)
        {
            for(int il = 0; il < lmax+1; il ++)
            {
                tmp >> uvalue_in[it][il];
            }
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_int(lmax);
    if(GlobalV::MY_RANK != 0)
    {
        uvalue_in.resize(ntype);
        for(int it = 0; it < ntype; it ++)
        {
            uvalue_in[it].resize(lmax+1);
        }
    }

    for(int it = 0; it < ntype; it ++)
    {
        for(int il = 0; il < lmax+1; il ++)
        {
            Parallel_Common::bcast_double(uvalue_in[it][il]);
        }
    }
#endif

    for(int it = 0; it < ntype; it ++)
    {
        for(int il = 0; il < lmax+1; il ++)
        {
            GlobalV::ofs_running << uvalue_in[it][il];
        }
    }    

}

std::string LCAO_DftU_New::scan_file(std::ifstream &ifs, std::string pattern)
{
    std::string line;

    while (!ifs.eof())
    {
        getline(ifs,line);
        if (line.find(pattern) != std::string::npos)
        {
            //replace all quotation marks by space
            //to make it easier for later operations
            std::replace( line.begin(), line.end(), '"', ' ');

            return line;
        }
    }

    ModuleBase::WARNING_QUIT("Paw_Element::scan_file","pattern not found in xml file!");
    return 0;    
}