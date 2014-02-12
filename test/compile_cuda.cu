#include "base/vector3d.h"
#include "base/soa3d.h"
#include "base/transformation_matrix.h"
#include "backend/cuda_backend.cuh"
#include "volumes/kernel/box_kernel.h"
#include "volumes/logical_volume.h"
#include "volumes/box.h"

using namespace vecgeom;

__device__
void compile_cuda() {
  CudaDouble scalar;
  Vector3D<double> scalar_v;
  Vector3D<CudaDouble > vector_v;
  SOA3D<CudaDouble> soa;
  TransformationMatrix<double> matrix;
  CudaBool output_inside;
  CudaDouble output_distance;
  BoxInside<translation::kOrigin, rotation::kIdentity, kCuda>(
    scalar_v, matrix, vector_v, &output_inside
  );
  BoxDistanceToIn<translation::kOrigin, rotation::kIdentity, kCuda>(
    scalar_v, matrix, vector_v, vector_v, scalar, &output_distance
  );
}