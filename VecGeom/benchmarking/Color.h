/// \file Color.h
/// \author Andrei Gheata

#ifndef VECGEOM_BENCHMARKING_COLOR_H_
#define VECGEOM_BENCHMARKING_COLOR_H_

#include <VecGeom/base/Global.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

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
    fComp.red   = 255 * red;
    fComp.green = 255 * green;
    fComp.blue  = 255 * blue;
    fComp.alpha = 255 * alpha;
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
  int GetColor() const { return fColor; }

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
    mdiff                 = maxval - minval;
    msum                  = maxval + minval;
    light                 = 0.5 * msum;
    if (maxval != minval) {
      rnorm = (maxval - Red()) / mdiff;
      gnorm = (maxval - Green()) / mdiff;
      bnorm = (maxval - Blue()) / mdiff;
    } else {
      satur = hue = 0;
      return;
    }

    if (light < 0.5)
      satur = mdiff / msum;
    else
      satur = mdiff / (2.0 - msum);

    if (Red() == maxval)
      hue = 60.0 * (6.0 + bnorm - gnorm);
    else if (Green() == maxval)
      hue = 60.0 * (2.0 + rnorm - bnorm);
    else
      hue = 60.0 * (4.0 + gnorm - rnorm);

    if (hue > 360) hue = hue - 360;
  }

  void SetHLS(float hue, float light, float satur)
  {
    float rh, rl, rs, rm1, rm2;
    rh = rl = rs = 0;
    if (hue > 0) {
      rh = hue;
      if (rh > 360) rh = 360;
    }
    if (light > 0) {
      rl = light;
      if (rl > 1) rl = 1;
    }
    if (satur > 0) {
      rs = satur;
      if (rs > 1) rs = 1;
    }

    if (rl <= 0.5)
      rm2 = rl * (1.0 + rs);
    else
      rm2 = rl + rs - rl * rs;
    rm1 = 2.0 * rl - rm2;

    if (!rs) {
      fComp.red   = 255 * rl;
      fComp.green = 255 * rl;
      fComp.blue  = 255 * rl;
      return;
    }

    auto HLStoRGB1 = [](float rn1, float rn2, float huei) {
      float hue = huei;
      if (hue > 360) hue = hue - 360;
      if (hue < 0) hue = hue + 360;
      if (hue < 60) return rn1 + (rn2 - rn1) * hue / 60;
      if (hue < 180) return rn2;
      if (hue < 240) return rn1 + (rn2 - rn1) * (240 - hue) / 60;
      return rn1;
    };

    fComp.red   = 255 * HLStoRGB1(rm1, rm2, rh + 120);
    fComp.green = 255 * HLStoRGB1(rm1, rm2, rh);
    fComp.blue  = 255 * HLStoRGB1(rm1, rm2, rh - 120);
  }
};

VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE Color_t operator+(Color_t const &left, Color_t const &right)
{
  Color_t color(left);
  color += right;
  return color;
}

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_BENCHMARKING_COLOR_H_
