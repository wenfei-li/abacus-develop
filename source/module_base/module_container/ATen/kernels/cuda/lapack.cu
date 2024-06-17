#include <ATen/kernels/lapack.h>
#include <base/third_party/lapack.h>

#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuda_runtime.h>
#include <thrust/complex.h>

namespace container {
namespace kernels {


static cusolverDnHandle_t cusolver_handle = nullptr;

void createGpuSolverHandle() {
    if (cusolver_handle == nullptr) {
        cusolverErrcheck(cusolverDnCreate(&cusolver_handle));
    }
}

void destroyGpuSolverHandle() {
    if (cusolver_handle != nullptr) {
        cusolverErrcheck(cusolverDnDestroy(cusolver_handle));
        cusolver_handle = nullptr;
    }
}

template <typename T>
__global__ void set_matrix_kernel(
    const char uplo,
    T* A,
    const int dim)
{
    int bid = blockIdx.x;
    int tid = threadIdx.x;

    for (int ii = tid; ii < bid + 1; ii += THREADS_PER_BLOCK) {
        if (uplo == 'L') {
            A[ii * dim + bid + 1] = static_cast<T>(0);
        }
        else {
            A[(bid + 1) * dim + ii] = static_cast<T>(0);
        }
    }
}

template <typename T>
struct set_matrix<T, DEVICE_GPU> {
    using Type = typename GetTypeThrust<T>::type;
    void operator() (
        const char& uplo,
        T* A,
        const int& dim)
    {
        set_matrix_kernel<Type><<<dim - 1, THREADS_PER_BLOCK>>>(
            uplo, reinterpret_cast<Type*>(A), dim);

        cudaCheckOnDebug();
    }
};

template <typename T>
struct lapack_trtri<T, DEVICE_GPU> {
    void operator()(
        const char& uplo,
        const char& diag,
        const int& dim,
        T* Mat,
        const int& lda)
    {
        // TODO: trtri is not implemented in this method yet
        // Cause the trtri in cuSolver is not stable for ABACUS!
        //cuSolverConnector::trtri(cusolver_handle, uplo, diag, dim, Mat, lda);
        cuSolverConnector::potri(cusolver_handle, uplo, diag, dim, Mat, lda);
    }
};

template <typename T>
struct lapack_potrf<T, DEVICE_GPU> {
    void operator()(
        const char& uplo,
        const int& dim,
        T* Mat,
        const int& lda)
    {
        cuSolverConnector::potrf(cusolver_handle, uplo, dim, Mat, dim);
    }
};

template <typename T>
struct lapack_dnevd<T, DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const char& jobz,
        const char& uplo,
        T* Mat,
        const int& dim,
        Real* eigen_val)
    {
        cuSolverConnector::dnevd(cusolver_handle, jobz, uplo, dim, Mat, dim, eigen_val);
    }
};

template <typename T>
struct lapack_dngvd<T, DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int& itype,
        const char& jobz,
        const char& uplo,
        T* Mat_A,
        T* Mat_B,
        const int& dim,
        Real* eigen_val)
    {
        cuSolverConnector::dngvd(cusolver_handle, itype, jobz, uplo, dim, Mat_A, dim, Mat_B, dim, eigen_val);
    }
};

template struct set_matrix<float,  DEVICE_GPU>;
template struct set_matrix<double, DEVICE_GPU>;
template struct set_matrix<std::complex<float>,  DEVICE_GPU>;
template struct set_matrix<std::complex<double>, DEVICE_GPU>;

template struct lapack_trtri<float,  DEVICE_GPU>;
template struct lapack_trtri<double, DEVICE_GPU>;
template struct lapack_trtri<std::complex<float>,  DEVICE_GPU>;
template struct lapack_trtri<std::complex<double>, DEVICE_GPU>;

template struct lapack_potrf<float,  DEVICE_GPU>;
template struct lapack_potrf<double, DEVICE_GPU>;
template struct lapack_potrf<std::complex<float>,  DEVICE_GPU>;
template struct lapack_potrf<std::complex<double>, DEVICE_GPU>;

template struct lapack_dnevd<float,  DEVICE_GPU>;
template struct lapack_dnevd<double, DEVICE_GPU>;
template struct lapack_dnevd<std::complex<float>,  DEVICE_GPU>;
template struct lapack_dnevd<std::complex<double>, DEVICE_GPU>;

template struct lapack_dngvd<float,  DEVICE_GPU>;
template struct lapack_dngvd<double, DEVICE_GPU>;
template struct lapack_dngvd<std::complex<float>,  DEVICE_GPU>;
template struct lapack_dngvd<std::complex<double>, DEVICE_GPU>;

} // namespace kernels
} // namespace container