#include "overlap_new.h"

#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

template <typename TK, typename TR>
hamilt::OverlapNew<hamilt::OperatorLCAO<TK, TR>>::OverlapNew(LCAO_Matrix* LM_in,
                                                             const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                                             hamilt::HContainer<TR>* hR_in,
                                                             std::vector<TK>* hK_in,
                                                             hamilt::HContainer<TR>* SR_in,
                                                             std::vector<TK>* SK_pointer_in,
                                                             const UnitCell* ucell_in,
                                                             Grid_Driver* GridD_in,
                                                             const ORB_gen_tables* uot,
                                                             const Parallel_Orbitals* paraV)
    : hamilt::OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in), uot_(uot)
{
    this->cal_type = lcao_overlap;
    this->ucell = ucell_in;
    this->SR = SR_in;
    this->SK_pointer = SK_pointer_in;
#ifdef __DEBUG
    assert(this->ucell != nullptr);
    assert(this->SR != nullptr);
    assert(this->SK_pointer != nullptr);
#endif
    // initialize SR to allocate sparse overlap matrix memory
    this->initialize_SR(GridD_in, paraV);
}

// initialize_SR()
template <typename TK, typename TR>
void hamilt::OverlapNew<hamilt::OperatorLCAO<TK, TR>>::initialize_SR(Grid_Driver* GridD, const Parallel_Orbitals* paraV)
{
    ModuleBase::TITLE("OverlapNew", "initialize_SR");
    ModuleBase::timer::tick("OverlapNew", "initialize_SR");
    for (int iat1 = 0; iat1 < ucell->nat; iat1++)
    {
        auto tau1 = ucell->get_tau(iat1);
        int T1, I1;
        ucell->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD->Find_atom(*ucell, tau1, T1, I1, &adjs);
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            int iat2 = ucell->itia2iat(T2, I2);
            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            // choose the real adjacent atoms
            const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
            // Note: the distance of atoms should less than the cutoff radius, 
            // When equal, the theoretical value of matrix element is zero, 
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (this->ucell->cal_dtau(iat1, iat2, R_index).norm() * this->ucell->lat0
                >= orb.Phi[T1].getRcut() + orb.Phi[T2].getRcut())
            {
                continue;
            }
            hamilt::AtomPair<TR> tmp(iat1, iat2, R_index, paraV);
            SR->insert_pair(tmp);
        }
    }
    // allocate the memory of BaseMatrix in SR, and set the new values to zero
    SR->allocate(nullptr, true);
    ModuleBase::timer::tick("OverlapNew", "initialize_SR");
}

template <typename TK, typename TR>
void hamilt::OverlapNew<hamilt::OperatorLCAO<TK, TR>>::calculate_SR()
{
    ModuleBase::TITLE("OverlapNew", "calculate_SR");
    ModuleBase::timer::tick("OverlapNew", "calculate_SR");
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iap = 0; iap < this->SR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<TR>& tmp = this->SR->get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        const Parallel_Orbitals* paraV = tmp.get_paraV();

        for (int iR = 0; iR < tmp.get_R_size(); ++iR)
        {
            const ModuleBase::Vector3<int> R_index = tmp.get_R_index(iR);
            auto dtau = ucell->cal_dtau(iat1, iat2, R_index);
            TR* data_pointer = tmp.get_pointer(iR);
            this->cal_SR_IJR(iat1, iat2, paraV, dtau, data_pointer);
        }
    }
    // if TK == double, then SR should be fixed to gamma case
    // judge type of TK equal to double
    if (std::is_same<TK, double>::value)
    {
        this->SR->fix_gamma();
    }
    ModuleBase::timer::tick("OverlapNew", "calculate_SR");
}

// cal_SR_IJR()
template <typename TK, typename TR>
void hamilt::OverlapNew<hamilt::OperatorLCAO<TK, TR>>::cal_SR_IJR(const int& iat1,
                                                                  const int& iat2,
                                                                  const Parallel_Orbitals* paraV,
                                                                  const ModuleBase::Vector3<double>& dtau,
                                                                  TR* data_pointer)
{
    const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
    // ---------------------------------------------
    // get info of orbitals of atom1 and atom2 from ucell
    // ---------------------------------------------
    int T1, I1;
    this->ucell->iat2iait(iat1, &I1, &T1);
    int T2, I2;
    this->ucell->iat2iait(iat2, &I2, &T2);
    Atom& atom1 = this->ucell->atoms[T1];
    Atom& atom2 = this->ucell->atoms[T2];

    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();

    const int* iw2l1 = atom1.iw2l;
    const int* iw2n1 = atom1.iw2n;
    const int* iw2m1 = atom1.iw2m;
    const int* iw2l2 = atom2.iw2l;
    const int* iw2n2 = atom2.iw2n;
    const int* iw2m2 = atom2.iw2m;
#ifndef USE_NEW_TWO_CENTER
    // ---------------------------------------------
    // get tau1 (in cell <0,0,0>) and tau2 (in cell R)
    // in principle, only dtau is needed in this function
    // snap_psipsi should be refactored to use dtau directly
    // ---------------------------------------------
    const ModuleBase::Vector3<double>& tau1 = this->ucell->get_tau(iat1);
    const ModuleBase::Vector3<double> tau2 = tau1 + dtau;
#endif
    // ---------------------------------------------
    // calculate the overlap matrix for each pair of orbitals
    // ---------------------------------------------
    double olm[3] = {0, 0, 0};
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    const int step_trace = col_indexes.size() + 1;
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const int iw1 = row_indexes[iw1l] / npol;
        const int L1 = iw2l1[iw1];
        const int N1 = iw2n1[iw1];
        const int m1 = iw2m1[iw1];
#ifdef USE_NEW_TWO_CENTER
        int M1 = (m1 % 2 == 0) ? -m1/2 : (m1+1)/2;
#endif
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const int iw2 = col_indexes[iw2l] / npol;
            const int L2 = iw2l2[iw2];
            const int N2 = iw2n2[iw2];
            const int m2 = iw2m2[iw2];
#ifdef USE_NEW_TWO_CENTER
            //=================================================================
            //          new two-center integral (temporary)
            //=================================================================
            // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
            int M2 = (m2 % 2 == 0) ? -m2/2 : (m2+1)/2;
            uot_->two_center_bundle->overlap_orb->calculate(T1, L1, N1, M1,
                    T2, L2, N2, M2, dtau * this->ucell->lat0, olm);
#else
            uot_->snap_psipsi(orb, // orbitals
                            olm,
                            0,
                            'S', // olm, job of derivation, dtype of Operator
                            tau1,
                            T1,
                            L1,
                            m1,
                            N1, // info of atom1
                            tau2,
                            T2,
                            L2,
                            m2,
                            N2 // info of atom2
            );
#endif
            for (int ipol = 0; ipol < npol; ipol++)
            {
                data_pointer[ipol * step_trace] += olm[0];
            }
            data_pointer += npol;
        }
        data_pointer += (npol - 1) * col_indexes.size();
    }
}

// contributeHR()
template <typename TK, typename TR>
void hamilt::OverlapNew<hamilt::OperatorLCAO<TK, TR>>::contributeHR()
{
    if (this->SR_fixed_done)
    {
        return;
    }
    this->calculate_SR();
    this->SR_fixed_done = true;
}

// contributeHk()
template <typename TK, typename TR>
void hamilt::OverlapNew<hamilt::OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    // if k vector is not changed, then do nothing and return, only for gamma_only case
    if (this->kvec_d[ik] == this->kvec_d_old && std::is_same<TK, double>::value)
    {
        return;
    }
    ModuleBase::TITLE("OverlapNew", "contributeHk");
    ModuleBase::timer::tick("OverlapNew", "contributeHk");
    // set SK to zero and then calculate SK for each k vector
    ModuleBase::GlobalFunc::ZEROS(this->SK_pointer->data(), this->SK_pointer->size());
    if(ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
    {
        const int nrow = this->SR->get_atom_pair(0).get_paraV()->get_row_size();
        hamilt::folding_HR(*this->SR, this->SK_pointer->data(), this->kvec_d[ik], nrow, 1);
    }
    else
    {
        const int ncol = this->SR->get_atom_pair(0).get_paraV()->get_col_size();
        hamilt::folding_HR(*this->SR, this->SK_pointer->data(), this->kvec_d[ik], ncol, 0);
    }
    // update kvec_d_old
    this->kvec_d_old = this->kvec_d[ik];

    ModuleBase::timer::tick("OverlapNew", "contributeHk");
}

template<typename TK, typename TR>
TK* hamilt::OverlapNew<hamilt::OperatorLCAO<TK, TR>>::getSk()
{
    if(this->SK_pointer != nullptr)
    {
        return this->SK_pointer->data();
    }
    return nullptr;
}

template class hamilt::OverlapNew<hamilt::OperatorLCAO<double, double>>;
template class hamilt::OverlapNew<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::OverlapNew<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;
