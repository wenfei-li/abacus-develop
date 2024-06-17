#ifndef CUDA_TOOLS_CUH
#define CUDA_TOOLS_CUH
#include <assert.h> // for assert
#include <cublas_v2.h>
#include <cuda.h> // for CUDA_VERSION
#include <cuda_runtime.h>

#include <fstream>
#include <iostream>
#include <sstream>

#define checkCuda(val) check(val, #val, __FILE__, __LINE__)
#define checkCudaLastError() __checkCudaLastError(__FILE__, __LINE__)

cudaError_t check(cudaError_t result, const char *const func, const char *const file, const int line);
cudaError_t __checkCudaLastError(const char *file, const int line);

void dump_cuda_array_to_file(double* cuda_array,
                             int width,
                             int hight,
                             const std::string& filename);

static inline int ceil_div(int a, int b)
{
    return (a + b - 1) / b;
}

/*
 * @brief: A simple wrapper for cudaMalloc and cudaFree, sync and async CUDA
 * memory copy
 * @param: T: the type of the data
 *
 * @note:
 * Manual management of CUDA memory is a very delicate task; complex pointers
 * and malloc/free operations make it easy for us to encounter memory bugs. The
 * severity of the issues increases significantly when introducing multi-node,
 * multi-GPU, and multi-stream parallelism.
 * Debugging after encountering bugs is also very difficult, finding the leaking
 * pointer from dozens of variables can be quite a headache.
 * Therefore, considering that our use and management of memory have some
 * homogeneity, we have abstracted these needs into the following encapsulations
 * to reduce the cost of maintenance and development. The memory is allocated in
 * the constructor and freed in the destructor.
 *
 * The following interface is primarily designed for the following requirements:
 * 1. We need to split a large task into multiple subtasks to run on multiple
 *    streams across multiple GPUs on multiple nodes.
 * 2. It is necessary to allocate memory of the same shape on both host and
 * device.
 * 3. Data copying between host and device sync or async is required.
 */

template <typename T>
class Cuda_Mem_Wrapper
{
  public:

    Cuda_Mem_Wrapper();
    Cuda_Mem_Wrapper(int one_stream_size,
                     int one_stream_size_aligned,
                     int stream_number = 1,
                     bool malloc_host = true);
    Cuda_Mem_Wrapper(int one_stream_size,
                     int stream_number = 1,
                     bool malloc_host = true);

    Cuda_Mem_Wrapper(const Cuda_Mem_Wrapper& other) = delete;
    Cuda_Mem_Wrapper& operator=(const Cuda_Mem_Wrapper& other) = delete;
    Cuda_Mem_Wrapper(Cuda_Mem_Wrapper&& other) noexcept;
    Cuda_Mem_Wrapper& operator=(Cuda_Mem_Wrapper&& other) noexcept;
    
    ~Cuda_Mem_Wrapper();
    void copy_host_to_device_sync(const int stream_id = 0);
    void copy_host_to_device_async(const cudaStream_t stream, const int stream_id);
    void copy_host_to_device_async(const cudaStream_t stream, const int stream_id, const int size);
    void copy_device_to_host_sync(const int stream_id = 0);
    void copy_device_to_host_async(const cudaStream_t stream, const int stream_id);
    void memset_device_sync(const int stream_id = 0, const int value = 0);
    void memset_device_async(const cudaStream_t stream, 
                             const int stream_id = 0,
                             const int value = 0);
    void memset_host(const int stream_id = 0, const int value = 0);
    T* get_device_pointer(const int stream_id = 0);
    T* get_host_pointer(const int stream_id = 0);
    void free_all();

  private:
    T* device_pointer;
    T* host_pointer;
    int one_stream_size;
    int one_stream_size_aligned;
    int stream_number;
    int total_size_aligned;
};

#endif // CUDA_TOOLS_CUH