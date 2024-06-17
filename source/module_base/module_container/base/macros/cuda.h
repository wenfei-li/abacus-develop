#ifndef BASE_MACROS_CUDA_H_
#define BASE_MACROS_CUDA_H_

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <thrust/complex.h>

#define THREADS_PER_BLOCK 256

template <typename T>
struct GetTypeThrust
{
    using type = T;
};

template <>
struct GetTypeThrust<std::complex<float>>
{
    using type = thrust::complex<float>; /**< The return type specialization for std::complex<float>. */
};

template <>
struct GetTypeThrust<std::complex<double>>
{
    using type = thrust::complex<double>; /**< The return type specialization for std::complex<float>. */
};

static inline cublasOperation_t GetCublasOperation(const char& trans)
{
    cublasOperation_t cutrans = {};
    if (trans == 'N')
    {
        cutrans = CUBLAS_OP_N;
    }
    else if (trans == 'T')
    {
        cutrans = CUBLAS_OP_T;
    }
    else if (trans == 'C')
    {
        cutrans = CUBLAS_OP_C;
    }
    return cutrans;
}

template <typename T>
struct GetTypeCuda
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_R_32F;
};
// Specializations of DataTypeToEnum for supported types.
template <>
struct GetTypeCuda<int>
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_R_32I;
};
template <>
struct GetTypeCuda<float>
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_R_32F;
};
template <>
struct GetTypeCuda<double>
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_R_64F;
};
template <>
struct GetTypeCuda<int64_t>
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_R_64I;
};
template <>
struct GetTypeCuda<std::complex<float>>
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_C_32F;
};
template <>
struct GetTypeCuda<std::complex<double>>
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_C_64F;
};

static inline cublasFillMode_t cublas_fill_mode(const char& uplo)
{
    if (uplo == 'U' || uplo == 'u')
        return CUBLAS_FILL_MODE_UPPER;
    else if (uplo == 'L' || uplo == 'l')
        return CUBLAS_FILL_MODE_LOWER;
    else
        throw std::runtime_error("cublas_fill_mode: unknown uplo");
}

static inline cublasDiagType_t cublas_diag_type(const char& diag)
{
    if (diag == 'U' || diag == 'u')
        return CUBLAS_DIAG_UNIT;
    else if (diag == 'N' || diag == 'n')
        return CUBLAS_DIAG_NON_UNIT;
    else
        throw std::runtime_error("cublas_diag_type: unknown diag");
}

static inline cusolverEigMode_t cublas_eig_mode(const char& jobz)
{
    if (jobz == 'N' || jobz == 'n')
        return CUSOLVER_EIG_MODE_NOVECTOR;
    else if (jobz == 'V' || jobz == 'v')
        return CUSOLVER_EIG_MODE_VECTOR;
    else
        throw std::runtime_error("cublas_eig_mode: unknown diag");
}

static inline cusolverEigType_t cublas_eig_type(const int& itype)
{
    if (itype == 1)
        return CUSOLVER_EIG_TYPE_1;
    else if (itype == 2)
        return CUSOLVER_EIG_TYPE_2;
    else
        throw std::runtime_error("cublas_eig_mode: unknown diag");
}

// cuSOLVER API errors
static const char* cusolverGetErrorEnum(cusolverStatus_t error)
{
    switch (error)
    {
    case CUSOLVER_STATUS_SUCCESS:
        return "CUSOLVER_STATUS_SUCCESS";
    case CUSOLVER_STATUS_NOT_INITIALIZED:
        return "CUSOLVER_STATUS_NOT_INITIALIZED";
    case CUSOLVER_STATUS_ALLOC_FAILED:
        return "CUSOLVER_STATUS_ALLOC_FAILED";
    case CUSOLVER_STATUS_INVALID_VALUE:
        return "CUSOLVER_STATUS_INVALID_VALUE";
    case CUSOLVER_STATUS_ARCH_MISMATCH:
        return "CUSOLVER_STATUS_ARCH_MISMATCH";
    case CUSOLVER_STATUS_MAPPING_ERROR:
        return "CUSOLVER_STATUS_MAPPING_ERROR";
    case CUSOLVER_STATUS_EXECUTION_FAILED:
        return "CUSOLVER_STATUS_EXECUTION_FAILED";
    case CUSOLVER_STATUS_INTERNAL_ERROR:
        return "CUSOLVER_STATUS_INTERNAL_ERROR";
    case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
        return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
    case CUSOLVER_STATUS_NOT_SUPPORTED:
        return "CUSOLVER_STATUS_NOT_SUPPORTED ";
    case CUSOLVER_STATUS_ZERO_PIVOT:
        return "CUSOLVER_STATUS_ZERO_PIVOT";
    case CUSOLVER_STATUS_INVALID_LICENSE:
        return "CUSOLVER_STATUS_INVALID_LICENSE";
    default:
        return "Unknown cusolverStatus_t message";
    }
}

inline void cusolverAssert(cusolverStatus_t code, const char* file, int line)
{
    if (code != CUSOLVER_STATUS_SUCCESS)
    {
        fprintf(stderr, " Unexpected cuSOLVER Error: %s %s %d\n", cusolverGetErrorEnum(code), file, line);
        exit(code);
    }
}

// cuSOLVER API errors
static const char* cublasGetErrorEnum(cublasStatus_t error)
{
    switch (error)
    {
    case CUBLAS_STATUS_SUCCESS:
        return "CUBLAS_STATUS_SUCCESS";
    case CUBLAS_STATUS_NOT_INITIALIZED:
        return "CUBLAS_STATUS_NOT_INITIALIZED";
    case CUBLAS_STATUS_ALLOC_FAILED:
        return "CUBLAS_STATUS_ALLOC_FAILED";
    case CUBLAS_STATUS_INVALID_VALUE:
        return "CUBLAS_STATUS_INVALID_VALUE";
    case CUBLAS_STATUS_ARCH_MISMATCH:
        return "CUBLAS_STATUS_ARCH_MISMATCH";
    case CUBLAS_STATUS_MAPPING_ERROR:
        return "CUBLAS_STATUS_MAPPING_ERROR";
    case CUBLAS_STATUS_EXECUTION_FAILED:
        return "CUBLAS_STATUS_EXECUTION_FAILED";
    case CUBLAS_STATUS_INTERNAL_ERROR:
        return "CUBLAS_STATUS_INTERNAL_ERROR";
    default:
        return "Unknown";
    }
}

inline void cublasAssert(cublasStatus_t res, const char* file, int line)
{
    if (res != CUBLAS_STATUS_SUCCESS)
    {
        fprintf(stderr, " Unexpected cuBLAS Error: %s %s %d\n", cublasGetErrorEnum(res), file, line);
        exit(res);
    }
}

#define cusolverErrcheck(res)                                                                                          \
    {                                                                                                                  \
        cusolverAssert((res), __FILE__, __LINE__);                                                                     \
    }

#define cublasErrcheck(res)                                                                                            \
    {                                                                                                                  \
        cublasAssert((res), __FILE__, __LINE__);                                                                       \
    }

// CUDA API errors
#define cudaErrcheck(res)                                                                                              \
    {                                                                                                                  \
        if (res != cudaSuccess)                                                                                        \
        {                                                                                                              \
            fprintf(stderr, " Unexpected Device Error %s:%d: %s, %s\n", __FILE__, __LINE__, cudaGetErrorName(res),     \
                    cudaGetErrorString(res));                                                                          \
            exit(res);                                                                                                 \
        }                                                                                                              \
    }

#ifdef __DEBUG
#define cudaCheckOnDebug() cudaErrcheck(cudaDeviceSynchronize())
#else
#define cudaCheckOnDebug()
#endif

#endif // BASE_MACROS_CUDA_H_