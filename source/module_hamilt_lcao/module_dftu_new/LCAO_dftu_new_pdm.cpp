//wenfei 2022-1-11
//This file contains subroutines for calculating pdm,
//which is defind as sum_mu,nu rho_mu,nu (<chi_mu|alpha><alpha|chi_nu>);

//It also contains subroutines for printing pdm for checking purpose

//There are 6 subroutines in this file:
//1. cal_projected_DM, which is used for calculating pdm for gamma point calculation
//2. cal_projected_DM_k, counterpart of 1, for multi-k
//3. check_projected_dm, which prints pdm to pdm.dat

#include "LCAO_dftu_new.h"
#include "module_base/vector3.h"
#include "module_base/timer.h"
#include "module_base/constants.h"

//this subroutine performs the calculation of projected density matrices
//pdm_m,m'=\sum_{mu,nu} rho_{mu,nu} <chi_mu|alpha_m><alpha_m'|chi_nu>
void LCAO_DftU_New::cal_projected_DM(const std::vector<ModuleBase::matrix> &dm, 
    const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver& GridD)
{
    ModuleBase::TITLE("LCAO_DftU_New", "cal_projected_DM");
    ModuleBase::timer::tick("LCAO_DftU_New","cal_projected_DM");

    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    if (GlobalV::init_chg == "file" && !this->init_pdm) //for DeePKS NSCF calculation 
    {
        std::ifstream ifs("pdm.dat");
        if (!ifs)
        {
            ModuleBase::WARNING_QUIT("LCAO_DftU_New::cal_projected_DM", "Can not find the file pdm.dat . Please do DeePKS SCF calculation first.");
        }
        for(int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for(int inl=0;inl<this->inlmax;inl++)
            {
                for(int ind=0;ind<pdm_size;ind++)
                {
                    double c;
                    ifs >> c;
                    pdm[is][inl][ind] = c;
                }
            }
        }
        this->init_pdm = true;
        return;
    }

    if(dm[0].nr == 0 && dm[0].nc ==0)
    {
        return;
    }

    for(int is = 0; is < GlobalV::NSPIN; is ++)
    {
        for(int inl=0;inl<inlmax;inl++)
        {
            ModuleBase::GlobalFunc::ZEROS(pdm[is][inl],pdm_size);
        }
    }

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = orb.Phi[T1].getRcut(); 

				for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GridD.getType(ad2);
					const int I2 = GridD.getNatom(ad2);
					const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
					const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
					
					const double Rcut_AO2 = orb.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * ucell.lat0;

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

					for (int iw1=0; iw1<nw1_tot; ++iw1)
					{
						const int iw1_all = start1 + iw1;
                        const int iw1_local = pv->global2local_col(iw1_all);
						if(iw1_local < 0)continue;
						const int iw1_0 = iw1/GlobalV::NPOL;

						for (int iw2=0; iw2<nw2_tot; ++iw2)
						{
							const int iw2_all = start2 + iw2;
                            const int iw2_local = pv->global2local_row(iw2_all);
							if(iw2_local < 0)continue;
							const int iw2_0 = iw2/GlobalV::NPOL;

                            std::vector<double> nlm1 = this->nlm_save[iat][ad1][iw1_all][0];
                            std::vector<double> nlm2 = this->nlm_save[iat][ad2][iw2_all][0];
                            assert(nlm1.size()==nlm2.size());

                            int ib=0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < 2 * L0 + 1;++m1)
                                    {
                                        for (int m2 = 0; m2 < 2 * L0 + 1; ++m2)
                                        {
                                            int ind = m1*nm + m2;
                                            for(int is = 0; is < dm.size(); ++is)
                                            {
                                                pdm[is][inl][ind] += dm[is](iw1_local, iw2_local)*nlm1[ib+m1]*nlm2[ib+m2];
                                            }
                                        }
                                    }
                                    ib+=nm;
                                }
                            }
                            assert(ib==nlm1.size());               
						}//iw2
					}//iw1
				}//ad2
			}//ad1
        }//I0
    }//T0

#ifdef __MPI
    for(int is = 0; is < GlobalV::NSPIN; is ++)
    {
        allsum_deepks(this->inlmax,pdm_size,this->pdm[is]);
    }
#endif
    ModuleBase::timer::tick("LCAO_DftU_New","cal_projected_DM"); 
    return;
}

void LCAO_DftU_New::cal_projected_DM_k(const std::vector<ModuleBase::ComplexMatrix>& dm,
    const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver& GridD,
    const int nks,
    const int * isk,
    const std::vector<ModuleBase::Vector3<double>> &kvec_d)
{

    ModuleBase::timer::tick("LCAO_DftU_New","cal_projected_DM_k");
    ModuleBase::TITLE("LCAO_DftU_New", "cal_projected_DM_k");

    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);

    if (GlobalV::init_chg == "file" && !this->init_pdm) //for DeePKS NSCF calculation 
    {
        std::ifstream ifs("pdm.dat");
        if (!ifs)
        {
            ModuleBase::WARNING_QUIT("LCAO_DftU_New::cal_projected_DM_k","Can not find the file pdm.dat . Please do DeePKS SCF calculation first.");
        }
        for(int is = 0; is < GlobalV::NSPIN; is ++)
        {
            for(int inl=0;inl<this->inlmax;inl++)
            {
                for(int ind=0;ind<pdm_size;ind++)
                {
                    double c;
                    ifs >> c;
                    pdm[is][inl][ind] = c;
                }
            }
            this->init_pdm = true;
        }
        return;
    }

    //check for skipping
    if(dm[0].nr == 0 && dm[0].nc ==0)
    {
        ModuleBase::timer::tick("LCAO_DftU_New","cal_projected_DM_k");
        return;
    }

    for(int is = 0; is < GlobalV::NSPIN; is ++)
    {
        for(int inl=0;inl<inlmax;inl++)
        {
            ModuleBase::GlobalFunc::ZEROS(pdm[is][inl],pdm_size);
        }
    }

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1,I1);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = orb.Phi[T1].getRcut();

                ModuleBase::Vector3<double> dR1(GridD.getBox(ad1).x, GridD.getBox(ad1).y, GridD.getBox(ad1).z); 

				for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GridD.getType(ad2);
					const int I2 = GridD.getNatom(ad2);
                    const int ibt2 = ucell.itia2iat(T2,I2);
					const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
					const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
                    ModuleBase::Vector3<double> dR2(GridD.getBox(ad2).x, GridD.getBox(ad2).y, GridD.getBox(ad2).z);
					
					const double Rcut_AO2 = orb.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * ucell.lat0;

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

					for (int iw1=0; iw1<nw1_tot; ++iw1)
					{
						const int iw1_all = start1 + iw1;
                        const int iw1_local = pv->global2local_col(iw1_all);
						if(iw1_local < 0)continue;
						const int iw1_0 = iw1/GlobalV::NPOL;
						for (int iw2=0; iw2<nw2_tot; ++iw2)
						{
							const int iw2_all = start2 + iw2;
                            const int iw2_local = pv->global2local_row(iw2_all);
							if(iw2_local < 0)continue;
							const int iw2_0 = iw2/GlobalV::NPOL;
 
                            std::vector<double> dm_current;
                            dm_current.resize(GlobalV::NSPIN);
                            std::vector<std::complex<double>> tmp;
                            tmp.resize(GlobalV::NSPIN);
                            for(int is = 0; is < GlobalV::NSPIN; is ++)
                            {
                                tmp[is]=0.0;
                            }
                            for(int ik=0;ik<nks;ik++)
                            {
                                const double arg = ( kvec_d[ik] * (dR1-dR2) ) * ModuleBase::TWO_PI;
                                const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );
                                const int is = isk[ik];
                                tmp[is] += dm[ik](iw1_local,iw2_local)*kphase;
                            }
                            for(int is = 0; is < GlobalV::NSPIN; is ++)
                            {
                                dm_current[is]=tmp[is].real();
                            }

                            key_tuple key_1(ibt1,dR1.x,dR1.y,dR1.z);
                            key_tuple key_2(ibt2,dR2.x,dR2.y,dR2.z);
                            std::vector<double> nlm1 = this->nlm_save_k[iat][key_1][iw1_all][0];
                            std::vector<double> nlm2 = this->nlm_save_k[iat][key_2][iw2_all][0];
                            assert(nlm1.size()==nlm2.size());

                            int ib=0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < 2 * L0 + 1;++m1)
                                    {
                                        for (int m2 = 0; m2 < 2 * L0 + 1; ++m2)
                                        {
                                            int ind = m1*nm + m2;
                                            for(int is = 0; is < GlobalV::NSPIN; is ++)
                                            {
                                                pdm[is][inl][ind] += dm_current[is]*nlm1[ib+m1]*nlm2[ib+m2];
                                            }
                                        }
                                    }
                                    ib+=nm;
                                }
                            }
                            assert(ib==nlm1.size());               
						}//iw2
					}//iw1
				}//ad2
			}//ad1
        }//I0
    }//T0

#ifdef __MPI
    for(int is = 0; is < GlobalV::NSPIN; is ++)
    {
        allsum_deepks(this->inlmax,pdm_size,this->pdm[is]);
    }
#endif
    ModuleBase::timer::tick("LCAO_DftU_New","cal_projected_DM_k");
    return;
    
}

void LCAO_DftU_New::check_projected_dm(void)
{
    std::ofstream ofs("pdm.dat");
    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    ofs<<std::setprecision(10);
    for(int is = 0; is < GlobalV::NSPIN; is ++)
    {
        for(int inl=0;inl<inlmax;inl++)
        {
            for(int ind=0;ind<pdm_size;ind++)
            {
                ofs << pdm[is][inl][ind] << " ";
            }
            ofs << std::endl;
        }
    }
}