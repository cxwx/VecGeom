/// \file cuda/Backend.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_BACKEND_CUDA_INTERFACE_H_
#define VECGEOM_BACKEND_CUDA_INTERFACE_H_

#include "VecGeom/base/Global.h"

#ifdef VECGEOM_ENABLE_CUDA

#include "driver_types.h" // Required for cudaError_t type
#include "cuda_runtime.h"

#include <type_traits>

namespace vecgeom {

#ifdef VECCORE_CUDA

inline namespace cuda {

template <typename DataClass, typename... ArgsTypes>
__global__ void ConstructOnGpu(DataClass *gpu_ptr, ArgsTypes... params)
{
  new (gpu_ptr) DataClass(params...);
}

template <typename DataClass, typename... ArgsTypes>
__global__ void ConstructArrayOnGpu(DataClass *gpu_ptr, size_t nElements, ArgsTypes... params)
{

  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  unsigned int idx = tid;
  while (idx < nElements) {
    new (gpu_ptr + idx) DataClass(params...);
    idx += blockDim.x * gridDim.x;
  }
}

template <typename DataClass, typename... ArgsTypes>
void Generic_CopyToGpu(DataClass *const gpu_ptr, ArgsTypes... params)
{
  ConstructOnGpu<<<1, 1>>>(gpu_ptr, params...);
}

} // End cuda namespace

#else

namespace cuda {

template <typename Type>
Type *AllocateOnDevice();
template <typename DataClass, typename... ArgsTypes>
void Generic_CopyToGpu(DataClass *const gpu_ptr, ArgsTypes... params);

} // End cuda namespace

#endif

#ifdef VECCORE_CUDA
namespace cxx {
#else
inline namespace cxx {
#endif

cudaError_t CudaCheckError(const cudaError_t err);

cudaError_t CudaCheckError();

void CudaAssertError(const cudaError_t err);

void CudaAssertError();

cudaError_t CudaMalloc(void **ptr, unsigned size);

cudaError_t CudaCopyToDevice(void *tgt, void const *src, unsigned size);

cudaError_t CudaCopyFromDevice(void *tgt, void const *src, unsigned size);

cudaError_t CudaFree(void *ptr);

template <typename Type>
Type *AllocateOnDevice()
{
  Type *ptr;
  vecgeom::cxx::CudaAssertError(vecgeom::cxx::CudaMalloc((void **)&ptr, sizeof(Type)));
  return ptr;
}

template <typename Type>
Type *AllocateOnGpu(const unsigned int size)
{
  Type *ptr;
  vecgeom::cxx::CudaAssertError(CudaMalloc((void **)&ptr, size));
  return ptr;
}

template <typename Type>
Type *AllocateOnGpu()
{
  return AllocateOnGpu<Type>(sizeof(Type));
}

template <typename Type>
void FreeFromGpu(Type *const ptr)
{
  vecgeom::cxx::CudaAssertError(CudaFree(ptr));
}

template <typename Type>
void CopyToGpu(Type const *const src, Type *const tgt, const unsigned size)
{
  vecgeom::cxx::CudaAssertError(CudaCopyToDevice(tgt, src, size));
}

template <typename Type>
void CopyToGpu(Type const *const src, Type *const tgt)
{
  CopyToGpu<Type>(src, tgt, sizeof(Type));
}

template <typename Type>
void CopyFromGpu(Type const *const src, Type *const tgt, const unsigned size)
{
  vecgeom::cxx::CudaAssertError(CudaCopyFromDevice(tgt, src, size));
}

class DevicePtrBase {
  void *fPtr;
#ifdef DEBUG_DEVICEPTR
  size_t fAllocatedSize;
  bool fIncremented;
#endif

protected:
  DevicePtrBase(const DevicePtrBase &orig)
      : fPtr(orig.fPtr)
#ifdef DEBUG_DEVICEPTR
        ,
        fAllocatedSize(0), fIncremented(false)
#endif
  {
  }

  DevicePtrBase &operator=(const DevicePtrBase &orig)
  {
    fPtr = orig.fPtr;
#ifdef DEBUG_DEVICEPTR
    fAllocatedSize = orig.fAllocatedSize;
    fIncremented   = orig.fIncremented;
#endif
    return *this;
  }

  void MemcpyToDevice(const void *what, unsigned long nbytes)
  {
    if (nbytes) vecgeom::cxx::CudaAssertError(cudaMemcpy(fPtr, what, nbytes, cudaMemcpyHostToDevice));
  }

  void MemcpyToHostAsync(void *where, unsigned long nbytes, cudaStream_t stream)
  {
    vecgeom::cxx::CudaAssertError(cudaMemcpyAsync(where, fPtr, nbytes, cudaMemcpyDeviceToHost, stream));
  }

  VECCORE_ATT_HOST_DEVICE
  void *GetPtr() const { return fPtr; }

  void Free()
  {
    vecgeom::cxx::CudaAssertError(vecgeom::cxx::CudaFree((void *)fPtr));
#ifdef DEBUG_DEVICEPTR
    fAllocatedSize = 0;
#endif
  }

  void Increment(long add)
  {
    fPtr = (char *)fPtr + add;
#ifdef DEBUG_DEVICEPTR
    if (add) fIncremented = true;
#endif
  }

public:
  DevicePtrBase()
      : fPtr(0)
#ifdef DEBUG_DEVICEPTR
        ,
        fAllocatedSize(0), fIncremented(0)
#endif
  {
  }

  explicit DevicePtrBase(void *input)
      : fPtr(input)
#ifdef DEBUG_DEVICEPTR
        ,
        fAllocatedSize(0), fIncremented(0)
#endif
  {
  }

  ~DevicePtrBase() { /* does not own content per se */}

  void Malloc(unsigned long size)
  {
    vecgeom::cxx::CudaAssertError(vecgeom::cxx::CudaMalloc((void **)&fPtr, size));
#ifdef DEBUG_DEVICEPTR
    fAllocatedSize = size;
#endif
  }
};

template <typename T>
class DevicePtr;

template <typename Type, typename Derived = DevicePtr<Type>>
class DevicePtrImpl : public DevicePtrBase {
protected:
  DevicePtrImpl(const DevicePtrImpl & /* orig */) = default;
  DevicePtrImpl &operator=(const DevicePtrImpl & /*orig*/) = default;
  DevicePtrImpl()                                          = default;
  explicit DevicePtrImpl(void *input) : DevicePtrBase(input) {}
  ~DevicePtrImpl() = default;

public:
  void Allocate(unsigned long nelems = 1) { Malloc(nelems * Derived::SizeOf()); }

  void Deallocate() { Free(); }

  void ToDevice(const Type *what, unsigned long nelems = 1) { MemcpyToDevice(what, nelems * Derived::SizeOf()); }
  void FromDevice(Type *where, cudaStream_t stream)
  {
    // Async since we pass a stream.
    MemcpyToHostAsync(where, Derived::SizeOf(), stream);
  }
  void FromDevice(Type *where, unsigned long nelems, cudaStream_t stream)
  {
    // Async since we pass a stream.
    MemcpyToHostAsync(where, nelems * Derived::SizeOf(), stream);
  }

  VECCORE_ATT_HOST_DEVICE
  Type *GetPtr() const { return reinterpret_cast<Type *>(DevicePtrBase::GetPtr()); }

  VECCORE_ATT_HOST_DEVICE
  operator Type *() const { return GetPtr(); }

  Derived &operator++() // prefix ++
  {
    Increment(Derived::SizeOf());
    return *(Derived *)this;
  }

  Derived operator++(int) // postfix ++
  {
    Derived tmp(*(Derived *)this);
    Increment(Derived::SizeOf());
    return tmp;
  }

  Derived &operator+=(long len) // prefix ++
  {
    Increment(len * Derived::SizeOf());
    return *(Derived *)this;
  }
};

template <typename Type>
class DevicePtr : public DevicePtrImpl<Type> {
public:
  DevicePtr()                  = default;
  DevicePtr(const DevicePtr &) = default;
  DevicePtr &operator=(const DevicePtr &orig) = default;

  // should be taking a DevicePtr<void*>
  explicit DevicePtr(void *input) : DevicePtrImpl<Type>(input) {}

  // Need to go via the explicit route accepting all conversion
  // because the regular c++ compilation
  // does not actually see the declaration for the cuda version
  // (and thus can not determine the inheritance).
  template <typename inputType>
  explicit DevicePtr(DevicePtr<inputType> const &input) : DevicePtrImpl<Type>((void *)input)
  {
  }

  // Disallow conversion from const to non-const.
  DevicePtr(DevicePtr<const Type> const &input,
            typename std::enable_if<!std::is_const<Type>::value, Type>::type * = nullptr) = delete;

#ifdef VECCORE_CUDA
  // Allows implicit conversion from DevicePtr<Derived> to DevicePtr<Base>
  template <typename inputType, typename std::enable_if<std::is_base_of<Type, inputType>::value>::type * = nullptr>
  DevicePtr(DevicePtr<inputType> const &input) : DevicePtrImpl<Type>(input.GetPtr())
  {
  }

  // Disallow conversion from const to non-const.
  template <typename inputType, typename std::enable_if<std::is_base_of<Type, inputType>::value>::type * = nullptr>
  DevicePtr(DevicePtr<const inputType> const &input) = delete;
#endif

#ifdef VECCORE_CUDA
  template <typename... ArgsTypes>
  void Construct(ArgsTypes... params) const
  {
    ConstructOnGpu<<<1, 1>>>(this->GetPtr(), params...);
  }

  template <typename... ArgsTypes>
  void ConstructArray(size_t nElements, ArgsTypes... params) const
  {
    ConstructArrayOnGpu<<<nElements, 1>>>(this->GetPtr(), nElements, params...);
  }

  static size_t SizeOf() { return sizeof(Type); }

#else
  template <typename... ArgsTypes>
  void Construct(ArgsTypes... params) const;
  template <typename... ArgsTypes>
  void ConstructArray(size_t nElements, ArgsTypes... params) const;

  static size_t SizeOf();
#endif
};

template <typename Type>
class DevicePtr<const Type> : private DevicePtrImpl<const Type> {
public:
  DevicePtr()                  = default;
  DevicePtr(const DevicePtr &) = default;
  DevicePtr &operator=(const DevicePtr &orig) = default;

  // should be taking a DevicePtr<void*>
  explicit DevicePtr(void *input) : DevicePtrBase(input) {}

  // Need to go via the explicit route accepting all conversion
  // because the regular c++ compilation
  // does not actually see the declaration for the cuda version
  // (and thus can not determine the inheritance).
  template <typename inputType>
  explicit DevicePtr(DevicePtr<inputType> const &input) : DevicePtrImpl<const Type>((void *)input)
  {
  }

  // Implicit conversion from non-const to const.
  DevicePtr(DevicePtr<typename std::remove_const<Type>::type> const &input) : DevicePtrImpl<const Type>((void *)input)
  {
  }

#ifdef VECCORE_CUDA
  // Allows implicit conversion from DevicePtr<Derived> to DevicePtr<Base>
  template <typename inputType, typename std::enable_if<std::is_base_of<Type, inputType>::value>::type * = nullptr>
  DevicePtr(DevicePtr<inputType> const &input) : DevicePtrImpl<const Type>(input.GetPtr())
  {
  }
#endif

  VECCORE_ATT_HOST_DEVICE
  const Type *GetPtr() const { return reinterpret_cast<const Type *>(DevicePtrBase::GetPtr()); }

  VECCORE_ATT_HOST_DEVICE
  operator const Type *() const { return GetPtr(); }

#ifdef VECCORE_CUDA
  template <typename DataClass, typename... ArgsTypes>
  void Construct(ArgsTypes... params) const
  {
    ConstructOnGpu<<<1, 1>>>(*(*this), params...);
  }

  template <typename... ArgsTypes>
  void ConstructArray(size_t nElements, ArgsTypes... params) const
  {
    ConstructArrayOnGpu<<<nElements, 1>>>(this->GetPtr(), nElements, params...);
  }

  static size_t SizeOf() { return sizeof(Type); }

#else
  template <typename... ArgsTypes>
  void Construct(ArgsTypes... params) const;
  template <typename... ArgsTypes>
  void ConstructArray(size_t nElements, ArgsTypes... params) const;

  static size_t SizeOf();
#endif
};

} // End cxx namespace

} // End vecgeom namespace

#endif // VECGEOM_ENABLE_CUDA

#endif // VECGEOM_BACKEND_CUDA_INTERFACE_H_
