#include "parallel_reduce.h"

#include "module_base/global_function.h"
#include "parallel_global.h"

template <>
void Parallel_Reduce::reduce_all<int>(int& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

void Parallel_Reduce::reduce_int_diag(int& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_INT, MPI_SUM, DIAG_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<double>(double& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<float>(float& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<int>(int* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

void Parallel_Reduce::reduce_int_grid(int* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, GRID_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<double>(double* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

void Parallel_Reduce::reduce_double_grid(double* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, GRID_WORLD);
#endif
    return;
}

void Parallel_Reduce::reduce_double_diag(double* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, DIAG_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<float>(float& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_FLOAT, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<double>(double& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<double>(double* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, POOL_WORLD);
#endif
    return;
}

// (1) the value is same in each pool.
// (2) we need to reduce the value from different pool.
void Parallel_Reduce::reduce_double_allpool(const int& kpar, const int& nproc_in_pool, double& object)
{
    if (kpar == 1)
        return;
#ifdef __MPI
    double swap = object / nproc_in_pool;
    MPI_Allreduce(&swap, &object, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}

// (1) the value is same in each pool.
// (2) we need to reduce the value from different pool.
void Parallel_Reduce::reduce_double_allpool(const int& kpar, const int& nproc_in_pool, double* object, const int n)
{
    if (kpar == 1)
        return;
#ifdef __MPI
    double* swap = new double[n];
    for (int i = 0; i < n; i++)
    {
        swap[i] = object[i] / nproc_in_pool;
    }
    MPI_Allreduce(swap, object, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    delete[] swap;
#endif
}

template <>
void Parallel_Reduce::reduce_all<std::complex<double>>(std::complex<double>& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

// LiuXh add 2019-07-16
template <>
void Parallel_Reduce::reduce_all<std::complex<double>>(std::complex<double>* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<std::complex<double>>(std::complex<double>& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<std::complex<float>>(std::complex<float>* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_C_FLOAT_COMPLEX, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<std::complex<double>>(std::complex<double>* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#endif
    return;
}

void Parallel_Reduce::gather_int_all(int& v, int* all)
{
#ifdef __MPI
    assert(all != nullptr);
    MPI_Allgather(&v, 1, MPI_INT, all, 1, MPI_INT, MPI_COMM_WORLD);
#endif
    return;
}

void Parallel_Reduce::gather_min_int_all(const int& nproc, int& v)
{
#ifdef __MPI
    int* all = new int[nproc];
    MPI_Allgather(&v, 1, MPI_INT, all, 1, MPI_INT, MPI_COMM_WORLD);
    for (int i = 0; i < nproc; i++)
    {
        if (v > all[i])
        {
            v = all[i];
        }
    }
    delete[] all;
#endif
}

void Parallel_Reduce::gather_max_double_all(const int& nproc, double& v)
{
#ifdef __MPI
    double* value = new double[nproc];
    ModuleBase::GlobalFunc::ZEROS(value, nproc);
    MPI_Allgather(&v, 1, MPI_DOUBLE, value, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    for (int i = 0; i < nproc; i++)
    {
        if (v < value[i])
        {
            v = value[i];
        }
    }
    delete[] value;
#endif
}

void Parallel_Reduce::gather_max_double_pool(const int& nproc_in_pool, double& v)
{
#ifdef __MPI
    if (nproc_in_pool == 1)
        return;
    double* value = new double[nproc_in_pool];
    ModuleBase::GlobalFunc::ZEROS(value, nproc_in_pool);
    MPI_Allgather(&v, 1, MPI_DOUBLE, value, 1, MPI_DOUBLE, POOL_WORLD);
    for (int i = 0; i < nproc_in_pool; i++)
    {
        if (v < value[i])
        {
            v = value[i];
        }
    }
    delete[] value;
#endif
}

void Parallel_Reduce::gather_min_double_pool(const int& nproc_in_pool, double& v)
{
#ifdef __MPI
    if (nproc_in_pool == 1)
        return;
    double* value = new double[nproc_in_pool];
    ModuleBase::GlobalFunc::ZEROS(value, nproc_in_pool);
    MPI_Allgather(&v, 1, MPI_DOUBLE, value, 1, MPI_DOUBLE, POOL_WORLD);
    for (int i = 0; i < nproc_in_pool; i++)
    {
        if (v > value[i])
        {
            v = value[i];
        }
    }
    delete[] value;
#endif
}

void Parallel_Reduce::gather_min_double_all(const int& nproc, double& v)
{
#ifdef __MPI
    double* value = new double[nproc];
    ModuleBase::GlobalFunc::ZEROS(value, nproc);
    MPI_Allgather(&v, 1, MPI_DOUBLE, value, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    for (int i = 0; i < nproc; i++)
    {
        if (v > value[i])
        {
            v = value[i];
        }
    }
    delete[] value;
#endif
}