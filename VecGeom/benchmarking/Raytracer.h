/// \file Raytracer.h
/// \author Andrei Gheata (andrei.gheata@cern.ch)

#ifndef VECGEOM_RAYTRACER_H_
#define VECGEOM_RAYTRACER_H_

#include <list>
#include <VecGeom/base/Global.h>
#include <VecGeom/base/Vector3D.h>
#include <VecGeom/base/Transformation3D.h>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/navigation/NavigationState.h>
#include <VecGeom/navigation/NewSimpleNavigator.h>
#include <VecGeom/volumes/PlacedBox.h>
#include <VecGeom/volumes/PlacedVolume.h>
#include <VecGeom/volumes/utilities/VolumeUtilities.h>

namespace vecgeom {

enum ERTmodel { kRTxray, kRTspecular, kRTtransparent, kRTdiffuse };
enum ERTView { kRTVparallel, kRTVperspective };

union Color_t {
  static constexpr float kToFloat = 1. / 0xFF;
  unsigned int fColor; // color representation as unsigned integer
  struct {
    unsigned char red;
    unsigned char green;
    unsigned char blue;
    unsigned char alpha;
  } fComp;

  Color_t() : fColor(0) {}
  Color_t(unsigned int col) { fColor = col; }
  Color_t(unsigned char r, unsigned char g, unsigned char b, unsigned char a)
  {
    fComp.red   = r;
    fComp.green = g;
    fComp.blue  = b;
    fComp.alpha = a;
  }
  Color_t(float r, float g, float b, float a)
  {
    fComp.red   = 255 * r;
    fComp.green = 255 * g;
    fComp.blue  = b;
    fComp.alpha = 255 * a;
  }

  Color_t &operator+=(Color_t const &other)
  {
    if ((fComp.alpha == 0) && (other.fComp.alpha == 0)) return *this;
    float alpha = 1 - (1 - other.Alpha()) * (1 - Alpha()); // cannot be 0
    float red   = (other.Red() * other.Alpha() + Red() * Alpha() * (1 - other.Alpha())) / alpha;
    float green = (other.Green() * other.Alpha() + Green() * Alpha() * (1 - other.Alpha())) / alpha;
    float blue  = (other.Blue() * other.Alpha() + Blue() * Alpha() * (1 - other.Alpha())) / alpha;
    fComp.red   = red;
    fComp.green = green;
    fComp.blue  = blue;
    fComp.alpha = alpha;
    return *this;
  }

  Color_t &operator*=(float val)
  {
    using vecCore::math::Max;
    using vecCore::math::Min;
    float alpha = val * Alpha();
    alpha       = Max(alpha, 0.0f);
    alpha       = Min(alpha, 1.0f);
    fComp.alpha = 255 * alpha;
    return *this;
  }

  Color_t &operator/=(float val)
  {
    using vecCore::math::Max;
    using vecCore::math::Min;
    float alpha = Alpha() / val;
    alpha       = Max(alpha, 0.0f);
    alpha       = Min(alpha, 1.0f);
    fComp.alpha = 255 * alpha;
    return *this;
  }

  float Red() const { return kToFloat * fComp.red; }
  float Green() const { return kToFloat * fComp.green; }
  float Blue() const { return kToFloat * fComp.blue; }
  float Alpha() const { return kToFloat * fComp.alpha; }
};

VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE Color_t operator+(Color_t const &left, Color_t const &right)
{
  Color_t color(left);
  color += right;
  return color;
}

struct Ray_t {
  Vector3D<double> fPos;
  Vector3D<double> fDir;
  NavigationState *fCrtState   = nullptr;
  NavigationState *fNextState  = nullptr;
  VPlacedVolume const *fVolume = nullptr;
  int fNcrossed                = 0;
  Color_t fColor;
  bool fDone = false;
};

/// \brief Raytracing a logical volume content using a given model
///
/// In order to run a benchmark, a world volume must be provided to the
/// raytracer.
class Raytracer {

public:
  typedef typename std::vector<std::pair<Vector3D<Precision>, Vector3D<Precision>>> RayContainer;

private:
  using VPlacedVolume_t = cxx::VPlacedVolume const *;

  VPlacedVolume_t fWorld = nullptr;
  std::list<VPlacedVolume_t> fVolumes;
  int fVerbosity  = 0;         ///< Verbosity level
  int fNrays      = 0;         ///< Number of rays left to propagate
  int fSize_px    = 1024;      ///< Image pixel size in x
  int fSize_py    = 1024;      ///< Image pixel size in y
  double fScale   = 0;         ///< Scaling from pixels to world coordinates
  ERTmodel fModel = kRTxray;   ///< Selected RT model
  Vector3D<double> fSourceDir; ///< Light source direction
  double fShininess = 1.;      ///< Shininess exponent in the specular model
  Color_t fLightColor;         ///< Light color
  Color_t fObjColor;           ///< Object color
  Vector3D<double> fStart;     ///< Source ray rectangle center position
  Vector3D<double> fDir;       ///< Start direction of all rays in parallel view mode
  Vector3D<double> fUp;        ///< Up vector in the shooting rectangle plane
  Vector3D<double> fRight;     ///< Right vector in the shooting rectangle plane
  Vector3D<double> fLeftC;     ///< left-down corner of the ray shooting rectangle
  char *fNavStates = nullptr;  ///< Actual location in memory of the navigation states
  Ray_t *fRays     = nullptr;  ///< Data block for rays

public:
  /// \brief dummy constructor, requiring to set the camera parameters and world volume after
  Raytracer() {}

  /// \param world Mother volume containing daughters that will be benchmarked.
  /// \param source_position position of the source in world coordinates, rays are starting from a plane passing through
  /// it, normal vector pointing to the origin of the world reference frame \param up_vector the projection of this
  /// vector on the camera plane determines the 'up' direction of the image \param img_size_px image size on X in pixels
  /// \param img_size_py image size on Y in pixels
  Raytracer(VPlacedVolume const *world, Vector3D<double> const &source_position, Vector3D<double> const &up_vector,
            int img_size_px, int img_size_py);

  ~Raytracer();

  /// \brief set ray source parameters
  /// \param source_position position of the source in world coordinates, rays are starting from a plane passing through
  /// it, normal vector pointing to the origin of the world reference frame \param up_vector the projection of this
  /// vector on the camera plane determines the 'up' direction of the image
  void SetSource(Vector3D<double> const &source_position, Vector3D<double> const &up_vector);

  /// \brief set image size in pixels
  void SetImageSize(int img_size_px, int img_size_py);

  /// \brief set raytracing model
  void SetRTModel(ERTmodel model) { fModel = model; }

  /// \brief set raytracing model
  void SetLightSourceDir(Vector3D<double> const &pos) { fSourceDir = pos; }

  /// \brief set light color
  void SetLightColor(unsigned char r, unsigned char g, unsigned char b, unsigned char a)
  {
    fLightColor = Color_t(r, g, b, a);
  }

  /// \brief set object color
  void SetObjColor(unsigned char r, unsigned char g, unsigned char b, unsigned char a)
  {
    fObjColor = Color_t(r, g, b, a);
  }

  /// \brief Entry point to propagate all rays
  void PropagateRays();

  /// \brief Initialize a ray and propagate it inside the world volume if needed
  /// \param px x-axis pixel index
  /// \param py y-axis pixel index
  void StartRay(int iray);

  /// \brief Propagate a ray for a single step until next boundary
  /// \param px x-axis pixel index
  /// \param py y-axis pixel index
  /// \return 0 or 1, if ray has stopped or has to continue
  int PropagateOneStep(int iray);

  /// \brief Propagate all rays for a single step until reaching boundaries
  /// \return Number of surviving rays
  int PropagateAllOneStep();

  /// \return World whose daughters are benchmarked.
  VPlacedVolume_t GetWorld() const { return fWorld; }

  /// \param World volume containing daughters to be benchmarked.
  void SetWorld(VPlacedVolume_t world);

private:
  void ApplyRTmodel(Ray_t &ray, double step);
  void GenerateVolumePointers(VPlacedVolume_t vol);

#ifdef VECGEOM_ENABLE_CUDA
  /*
    void RunInsideCuda(Precision *posX, Precision *posY, Precision *posZ, bool *contains, Inside_t *inside);
    void RunToInCuda(Precision *posX, Precision *posY, Precision *posZ, Precision *dirX, Precision *dirY, Precision
    *dirZ, Precision *distances, Precision *safeties); void RunToOutCuda(Precision *posX, Precision *posY, Precision
    *posZ, Precision *dirX, Precision *dirY, Precision *dirZ, Precision *distances, Precision *safeties);
  */
  void GetVolumePointers(std::list<cxx::DevicePtr<cuda::VPlacedVolume>> &volumesGpu);
#endif
};

} // End namespace vecgeom

#endif // VECGEOM_RAYTRACER_H_
