#ifdef __CUSOLVERMP

#include "diago_cusolvermp.h"

#include "module_base/timer.h"

using complex = std::complex<double>;

namespace hsolver
{
template <typename T>
void DiagoCusolverMP<T>::diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoCusolverMP", "diag");
    hamilt::MatrixBlock<T> h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);

    std::vector<Real> eigen(GlobalV::NLOCAL, 0.0);

    MPI_Comm COMM_DIAG = MPI_COMM_WORLD; // use all processes
    {
        Diag_CusolverMP_gvd<T> es(COMM_DIAG, (const int)h_mat.row, (const int)h_mat.col, (const int*)h_mat.desc);

        ModuleBase::timer::tick("DiagoCusolverMP", "Diag_CusolverMP_gvd");
        es.generalized_eigenvector(h_mat.p, s_mat.p, eigen.data(), psi.get_pointer());
        ModuleBase::timer::tick("DiagoCusolverMP", "Diag_CusolverMP_gvd");
    }
    const int inc = 1;
    BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
}
template class DiagoCusolverMP<double>;
template class DiagoCusolverMP<complex>;

} // namespace hsolver

#endif