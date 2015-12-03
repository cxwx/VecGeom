#include <cmath>
#include <limits>
#include <cstdio>
#include <cstdlib>

#include <malloc.h>

#include <VecCore>
#include "Timer.h"

using namespace vecCore;

static const int N = (8 * 1024 * 1024);

// solve ax2 + bx + c = 0, return number of roots found

template <typename T> int quadsolve(T a, T b, T c, T &x1, T &x2) {
  T delta = b * b - 4.0 * a * c;

  if (delta < 0.0)
    return 0;

  if (delta < std::numeric_limits<T>::epsilon()) {
    x1 = x2 = -0.5 * b / a;
    return 1;
  }

  if (b >= 0.0) {
    x1 = -0.5 * (b + std::sqrt(delta)) / a;
    x2 = c / (a * x1);
  } else {
    x2 = -0.5 * (b - std::sqrt(delta)) / a;
    x1 = c / (a * x2);
  }

  return 2;
}

// SIMD version

// solve ax2 + bx + c = 0, return number of roots found
// solves many equation simultaneously, depending on value of Float_v::Size

template <class Backend>
void quadsolve_vc(typename Backend::Float_v const &a,
                  typename Backend::Float_v const &b,
                  typename Backend::Float_v const &c,
				  typename Backend::Float_v &x1, typename Backend::Float_v &x2,
                  typename Backend::Int_v &roots) {
  using Int_v = typename Backend::Int_v;
  using Float_t = typename Backend::Float_t;
  using Float_v = typename Backend::Float_v;

  roots = 0;
  Float_v epsilon = Float_v(std::numeric_limits<Float_t>::epsilon());
  Float_v delta = b * b - Float_v(4.0) * a * c;

  typename Float_v::Mask no_roots(delta < Float_v(0.0));
  typename Float_v::Mask two_roots(delta >= epsilon);

  roots((typename Int_v::Mask)two_roots) = 2;

  typename Float_v::Mask mask = (b >= Float_v(0.0));

  Float_v sign = Blend(mask, Float_v(-1.0), Float_v(1.0));
  Float_v root = (-b + sign * math::Sqrt(delta)) / (Float_v(2.0) * a);

  x1(two_roots && mask) = root;
  x2(two_roots && !mask) = root;

  root = c / (a * Blend(mask, x1, x2));

  x2(two_roots && mask) = root;
  x1(two_roots && !mask) = root;

  mask = !(no_roots || two_roots);

  if (IsEmpty(mask))
    return;

  roots((typename Int_v::Mask)mask) = 1;

  root = -Float_v(0.5) * b / a;
  x1(mask) = root;
  x2(mask) = root;
}

int main(int argc, char *argv[]) {
  float *a = (float *)memalign(64, N * sizeof(float));
  float *b = (float *)memalign(64, N * sizeof(float));
  float *c = (float *)memalign(64, N * sizeof(float));

  int *roots = (int *)memalign(64, N * sizeof(int));
  float *x1 = (float *)memalign(64, N * sizeof(float));
  float *x2 = (float *)memalign(64, N * sizeof(float));

  srand48(time(NULL));
  int index = (int)((N - 100) * drand48());

  for (int i = 0; i < N; i++) {
    a[i] = 10.0 * (drand48() - 0.5);
    b[i] = 10.0 * (drand48() - 0.5);
    c[i] = 50.0 * (drand48() - 0.5);
    x1[i] = 0.0;
    x2[i] = 0.0;
    roots[i] = 0;
  }

  Timer<> timer;

  for (int i = 0; i < N; i++) {
    roots[i] = quadsolve(a[i], b[i], c[i], x1[i], x2[i]);
  }

  double t = timer.Elapsed();

  // print random result to ensure scalar and vector backends give same result

  for (int i = index; i < index + 10; i++) {
    printf("%d: a = % 8.3f, b = % 8.3f, c = % 8.3f, roots = %d, x1 = % 8.3f, x2 = % 8.3f\n", i, a[i], b[i], c[i],
           roots[i], roots[i] > 0 ? x1[i] : 0, roots[i] > 1 ? x2[i] : 0);
  }

  printf("\nelapsed time = %.3lfms (scalar code)\n", t / 1.0e6);

  timer.Start();

  for (int i = 0; i < N; i += backend::Vector<float>::Float_v::Size) {
    quadsolve_vc<backend::Vector<float>>(
        (backend::Vector<float>::Float_v &)(a[i]), (backend::Vector<float>::Float_v &)(b[i]),
        (backend::Vector<float>::Float_v &)(c[i]), (backend::Vector<float>::Float_v &)(x1[i]),
        (backend::Vector<float>::Float_v &)(x2[i]), (backend::Vector<float>::Int_v &)(roots[i]));
  }

  t = timer.Elapsed();

  for (int i = index; i < index + 10; i++) {
    printf("%d: a = % 8.3f, b = % 8.3f, c = % 8.3f, roots = %d, x1 = % 8.3f, x2 = % 8.3f\n", i, a[i], b[i], c[i],
           roots[i], roots[i] > 0 ? x1[i] : 0, roots[i] > 1 ? x2[i] : 0);
  }

  printf("\nelapsed time = %.3lfms (vector backend)\n", t / 1.0e6);

  return 0;
}
