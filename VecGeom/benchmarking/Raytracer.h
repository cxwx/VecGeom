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
#include <VecGeom/volumes/PlacedBox.h>
#include <VecGeom/volumes/PlacedVolume.h>
#include <VecGeom/volumes/utilities/VolumeUtilities.h>

namespace vecgeom {

enum ERTmodel { kRTxray = 0, kRTspecular, kRTtransparent, kRTdiffuse };
enum ERTView { kRTVparallel = 0, kRTVperspective };

union Color_t {
  static constexpr float kToFloat = 1. / 0xFF;
  unsigned int fColor; // color representation as unsigned integer
  struct {
    unsigned char alpha;
    unsigned char blue;
    unsigned char green;
    unsigned char red;
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
    fComp.red   = 255*red;
    fComp.green = 255*green;
    fComp.blue  = 255*blue;
    fComp.alpha = 255*alpha;
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

  void MultiplyLightChannel(float fact)
  {
    float hue, light, satur;
    GetHLS(hue, light, satur);
    SetHLS(hue, fact * light, satur);
  }

  void GetHLS(float &hue, float &light, float &satur) const
  {
    float rnorm, gnorm, bnorm, msum, mdiff;

    float minval = Min(Red(), Green(), Blue());
    float maxval = Max(Red(), Green(), Blue());

    rnorm = gnorm = bnorm = 0;
    mdiff = maxval - minval;
    msum  = maxval + minval;
    light = 0.5 * msum;
    if (maxval != minval) {
      rnorm = (maxval - Red())/mdiff;
      gnorm = (maxval - Green())/mdiff;
      bnorm = (maxval - Blue())/mdiff;
    } else {
      satur = hue = 0;
      return;
    }

    if (light < 0.5)
      satur = mdiff/msum;
    else
      satur = mdiff/(2.0 - msum);

    if (Red() == maxval)
      hue = 60.0 * (6.0 + bnorm - gnorm);
    else if (Green() == maxval)
      hue = 60.0 * (2.0 + rnorm - bnorm);
    else
      hue = 60.0 * (4.0 + gnorm - rnorm);

    if (hue > 360)
      hue = hue - 360;
  }

  void SetHLS(float hue, float light, float satur)
  {
    float rh, rl, rs, rm1, rm2;
    rh = rl = rs = 0;
    if (hue   > 0) { rh = hue;   if (rh > 360) rh = 360; }
    if (light > 0) { rl = light; if (rl > 1)   rl = 1; }
    if (satur > 0) { rs = satur; if (rs > 1)   rs = 1; }

    if (rl <= 0.5)
      rm2 = rl*(1.0 + rs);
    else
      rm2 = rl + rs - rl*rs;
    rm1 = 2.0*rl - rm2;

    if (!rs) { fComp.red = 255*rl; fComp.green = 255*rl; fComp.blue = 255*rl; return; }

    auto HLStoRGB1 = [](float rn1, float rn2, float huei) {
      float hue = huei;
      if (hue > 360) hue = hue - 360;
      if (hue < 0)   hue = hue + 360;
      if (hue < 60 ) return rn1 + (rn2-rn1)*hue/60;
      if (hue < 180) return rn2;
      if (hue < 240) return rn1 + (rn2-rn1)*(240-hue)/60;
      return rn1;
    };

    fComp.red   = 255 * HLStoRGB1(rm1, rm2, rh+120);
    fComp.green = 255 * HLStoRGB1(rm1, rm2, rh);
    fComp.blue  = 255 * HLStoRGB1(rm1, rm2, rh-120);
  }
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
  int fVisDepth   = 1;         ///< Visible geometry depth
  double fScale   = 0;         ///< Scaling from pixels to world coordinates
  ERTmodel fModel = kRTxray;   ///< Selected RT model
  Vector3D<double> fSourceDir; ///< Light source direction
  double fShininess = 1.;      ///< Shininess exponent in the specular model
  Color_t fLightColor;         ///< Light color
  Color_t fObjColor;           ///< Object color
  Vector3D<double> fStart;     ///< Screen position
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
  /// \param screen_position position of the screen in world coordinates, rays are starting from a plane passing through
  /// it, normal vector pointing to the origin of the world reference frame \param up_vector the projection of this
  /// vector on the camera plane determines the 'up' direction of the image \param img_size_px image size on X in pixels
  /// \param img_size_py image size on Y in pixels
  Raytracer(VPlacedVolume const *world, Vector3D<double> const &screen_position, Vector3D<double> const &up_vector,
            int img_size_px, int img_size_py, ERTmodel model);

  ~Raytracer();

  /// \brief set image size in pixels
  void SetImageSize(int img_size_px, int img_size_py);

  /// \brief set raytracing model
  void SetRTModel(ERTmodel model) { fModel = model; }

  /// \brief set raytracing model
  void SetLightSourceDir(Vector3D<double> const &dir) { if (fSourceDir.Mag2() > 0) {fSourceDir = dir; fSourceDir.Normalize();} }

  /// \brief set light color
  void SetLightColor(unsigned int col) { fLightColor = col; }

  /// \brief set object color
  void SetObjColor(unsigned int col) { fObjColor = col; }

  /// \brief set visible depth
  void SetVisDepth(int depth) { fVisDepth = depth; }

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
  void CreateNavigators();
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
