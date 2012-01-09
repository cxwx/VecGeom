#ifndef PROPAGATOR_GLOBALS
#define PROPAGATOR_GLOBALS

#include "TStopwatch.h"
#include "GeantPropagator.h"

R__EXTERN GeantPropagator   *gPropagator;     // propagator class

// types
enum Species_t {kHadron, kLepton};

#endif
