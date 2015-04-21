/*
 *  ABBoxBenchmark.cpp
 *
 *  simple benchmark to play with aligned bounding boxes
 *  which could be a simple but efficient form of "voxelization"
 *
 *  bounding boxes are setup in some form of regular 3D grid of size N x N x N
 *
 *  benchmark:
 *  a) intersecting bounding boxes without caching of inverseray + no vectorization
 *  b) intersecting bounding boxes with caching of inverseray + no vectorization
 *  c) intersecting bounding boxes with caching + vectorization
 *
 *  Created on: 20.04.2015
 *      Author: swenzel
 */

#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"
#include "benchmarking/Benchmarker.h"
#include "volumes/kernel/BoxImplementation.h"
#include "volumes/utilities/VolumeUtilities.h"
#include "management/GeoManager.h"
#include "ArgParser.h"
#include "base/SOA3D.h"
#include "base/Stopwatch.h"
#include <iostream>
#include <vector>
#include <backend/vc/Backend.h>
#include <backend/Backend.h>

using namespace vecgeom;

int benchNoCachingNoVector( Vector3D<Precision> const & point, Vector3D<Precision> const & dir, std::vector<Vector3D<Precision> > const & );
int benchCachingNoVector(   Vector3D<Precision> const & point, Vector3D<Precision> const & dir, std::vector<Vector3D<Precision> > const & );
int benchCachingAndVector(  Vector3D<Precision> const & point, Vector3D<Precision> const & dir, Vector3D<kVc::precision_v> const *, int number );

#define N 20        // boxes per dimension
#define SZ 1000      // samplesize
double delta = 0.5; // if delta > 1 the boxes will overlap

// #define INNERTIMER

int main(){
    // number of boxes
    int numberofboxes = N*N*N;

    int code = (2 << 1) + (2 << 0) + (2 << 1);
    std::cerr << code << "\n";

    // setup AOS form of boxes
    std::vector<Vector3D<Precision> > uppercorners(numberofboxes);
    std::vector<Vector3D<Precision> > lowercorners(numberofboxes);

    // setup same in mixed array of corners ... upper-lower-upper-lower ...
    std::vector<Vector3D<Precision> > corners(2*numberofboxes);

    // setup SOA form of boxes -- the memory layout should probably rather be SOA6D
    Vector3D<kVc::precision_v> * VcCorners =  new Vector3D<kVc::precision_v>[ 2*numberofboxes/kVc::precision_v::Size ];

    int counter1=0;
    int counter2=0;
    for(int i=0; i<N; ++i){
        for(int j=0; j<N; ++j){
            for(int k=0; k<N; ++k){
                lowercorners[counter1] = Vector3D<Precision>( i, j, k );
                uppercorners[counter1] = Vector3D<Precision>( i + delta, j + delta, k + delta );

                corners[counter2] = lowercorners[counter1];
                counter2++;
                corners[counter2] = uppercorners[counter1];
                counter2++;

                counter1++;
            }
        }
    }

    // print boxes
    for(int i=0; i<numberofboxes; ++i)
    {
        // std::cerr << "# " << i << " lower " << lowercorners[i] << " " << uppercorners[i] << "\n";
    }

    // set up VcCorners
    counter2=0;
    for(int i=0; i < numberofboxes; i+=kVc::precision_v::Size)
    {
        Vector3D<kVc::precision_v> lower;
        Vector3D<kVc::precision_v> upper;
        // assign by components
        for( int k=0;k<kVc::precision_v::Size;++k ){
            lower.x()[k] = lowercorners[i+k].x();
            lower.y()[k] = lowercorners[i+k].y();
            lower.z()[k] = lowercorners[i+k].z();
            upper.x()[k] = uppercorners[i+k].x();
            upper.y()[k] = uppercorners[i+k].y();
            upper.z()[k] = uppercorners[i+k].z();
        }
        //std::cerr << lower << "\n";
        //std::cerr << upper << "\n";
        VcCorners[counter2++] = lower;
        VcCorners[counter2++] = upper;
    }
    std::cerr << "assigned " << counter2 << "Vc vectors\n";

    // constructing samples
    std::vector<Vector3D<Precision> > points(SZ);
    std::vector<Vector3D<Precision> > directions(SZ);
    for(int i=0;i<SZ;++i)
    {
        points[i] = Vector3D<Precision>( N*delta + 0.1, N*delta + 0.1, N*delta + 0.1 );
        directions[i] = volumeUtilities::SampleDirection();
    }

    Stopwatch timer;
    int hits=0;

    timer.Start();
    for(int i=0;i<SZ;++i){
        hits+=benchCachingNoVector( points[i], directions[i], corners );
    }
    timer.Stop();
    std::cerr << "Cached times and hit " << timer.Elapsed() << " " << hits << "\n";

    hits=0;
    timer.Start();
    for(int i=0;i<SZ;++i){
        hits+=benchNoCachingNoVector( points[i], directions[i], corners );
    }
    timer.Stop();
    std::cerr << "Ordinary times and hit " << timer.Elapsed() << " " << hits << "\n";

    hits=0;
    timer.Start();
    for(int i=0;i<SZ;++i){
        hits+=benchCachingAndVector( points[i], directions[i], VcCorners, numberofboxes/kVc::precision_v::Size );
    }
    timer.Stop();
    std::cerr << "Vector times and hit " << timer.Elapsed() << " " << hits << "\n";

    return 0;
}


int benchNoCachingNoVector( Vector3D<Precision> const & point, Vector3D<Precision> const & dir, std::vector<Vector3D<Precision> > const & corners ){
#ifdef INNERTIMER
    Stopwatch timer;
    timer.Start();
#endif
    int vecsize = corners.size() / 2;
    int hitcount = 0;
    for( auto box = 0; box < vecsize; ++box ){
         double distance = BoxImplementation<translation::kIdentity, rotation::kIdentity>::Intersect(
                &corners[2*box],
                point,
                dir,
                0, vecgeom::kInfinity );
         if( distance < vecgeom::kInfinity ) hitcount++;
    }
#ifdef INNERTIMER
    timer.Stop();
    std::cerr << "# ORDINARY hitting " << hitcount << "\n";
    std::cerr << "# ORDINARY timer " << timer.Elapsed() << "\n";
#endif
    return hitcount;
}

int benchCachingNoVector( Vector3D<Precision> const & point, Vector3D<Precision> const & dir, std::vector<Vector3D<Precision> > const & corners ){
#ifdef INNERTIMER
    Stopwatch timer;
    timer.Start();
#endif

    Vector3D<Precision> invdir(1./dir.x(), 1./dir.y(), 1./dir.z());
    int vecsize = corners.size() / 2;
    int hitcount = 0;
    int sign[3]; sign[0] = invdir.x() < 0; sign[1] = invdir.y() < 0; sign[2] = invdir.z() < 0;
    // interpret as binary number and do a switch statement
    // do a big switch statement here
   // int code = 2 << size[0] + 2 << size[1] + 2 << size[2];
    for( auto box = 0; box < vecsize; ++box ){
         double distance = BoxImplementation<translation::kIdentity, rotation::kIdentity>::IntersectCachedKernel2<kScalar>(
            &corners[2*box],
            point,
           invdir,
           sign[0],sign[1],sign[2],
            0, vecgeom::kInfinity );
            if( distance < vecgeom::kInfinity ) hitcount++;
        }

    //    switch( size[0] + size[1] + size[2] ){
//    case 0: {
//        for( auto box = 0; box < vecsize; ++box ){
//        double distance = BoxImplementation<translation::kIdentity, rotation::kIdentity>::IntersectCachedKernel<kScalar,0,0,0>(
//           &corners[2*box],
//           point,
//           invdir,
//           0, vecgeom::kInfinity );
//           if( distance < vecgeom::kInfinity ) hitcount++;
//         }       break; }
//    case 3: {
//        for( auto box = 0; box < vecsize; ++box ){
//                double distance = BoxImplementation<translation::kIdentity, rotation::kIdentity>::IntersectCachedKernel<kScalar,1,1,1>(
//                   &corners[2*box],
//                   point,
//                   invdir,
//                   0, vecgeom::kInfinity );
//                   if( distance < vecgeom::kInfinity ) hitcount++;
//                 }       break; }
//    default : std::cerr << "DEFAULT CALLED\n";
//    }
#ifdef INNERTIMER
    timer.Stop();
    std::cerr << "# CACHED hitting " << hitcount << "\n";
    std::cerr << "# CACHED timer " << timer.Elapsed() << "\n";
#endif
    return hitcount;
}


int benchCachingAndVector( Vector3D<Precision> const & point, Vector3D<Precision> const & dir, Vector3D<kVc::precision_v> const * corners, int vecsize ){
#ifdef INNERTIMER
    Stopwatch timer;
    timer.Start();
#endif

    Vector3D<Precision> invdir(1./dir.x(), 1./dir.y(), 1./dir.z());
    int hitcount = 0;


    int sign[3]; sign[0] = invdir.x() < 0; sign[1] = invdir.y() < 0; sign[2] = invdir.z() < 0;
    for( auto box = 0; box < vecsize; ++box ){
            kVc::precision_v distance = BoxImplementation<translation::kIdentity, rotation::kIdentity>::IntersectCachedKernel2<kVc>(
               &corners[2*box],
               point,
               invdir,sign[0],sign[1],sign[2],
               0, vecgeom::kInfinity );
               kVc::bool_v hit = distance < vecgeom::kInfinity;
               //std::cerr << hit << "\n";
               hitcount += hit.count();
    }
    // interpret as binary number and do a switch statement
    // do a big switch statement here
//    switch( size[0] + size[1] + size[2] ){
//    case 0: {
//        for( auto box = 0; box < vecsize; ++box ){
//        kVc::precision_v distance = BoxImplementation<translation::kIdentity, rotation::kIdentity>::IntersectCachedKernel<kVc,0,0,0>(
//           &corners[2*box],
//           point,
//           invdir,
//           0, vecgeom::kInfinity );
//           kVc::bool_v hit = distance < vecgeom::kInfinity;
//           //std::cerr << hit << "\n";
//           hitcount += hit.count();
//        }       break; }
//    case 3: {
//        for( auto box = 0; box < vecsize; ++box ){
//          kVc::precision_v distance = BoxImplementation<translation::kIdentity, rotation::kIdentity>::IntersectCachedKernel<kVc,1,1,1>(
//          &corners[2*box],
//          point,
//          invdir,
//          0, vecgeom::kInfinity );
//          kVc::bool_v hit = distance < vecgeom::kInfinity;
//          //std::cerr << hit << "\n";
//          hitcount += hit.count();
//    }       break; }
//    default : std::cerr << "DEFAULT CALLED\n";
//    }
#ifdef INNERTIMER
    timer.Stop();
    std::cerr << "# VECTOR hitting " << hitcount << "\n";
    std::cerr << "# VECTOR timer " << timer.Elapsed() << "\n";
#endif
    return hitcount;
}

