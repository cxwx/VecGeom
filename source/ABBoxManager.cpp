/*
 * ABBoxManager.cpp
 *
 *  Created on: 24.04.2015
 *      Author: swenzel
 */

#include "navigation/ABBoxNavigator.h"
#include "volumes/UnplacedBox.h"
#include "utilities/Visualizer.h"
#undef NDEBUG
#include <cassert>

namespace vecgeom {
inline namespace cxx {

void ABBoxManager::ComputeABBox( VPlacedVolume const * pvol, ABBox_t * lowerc, ABBox_t * upperc ) {
        // idea: take the 8 corners of the bounding box in the reference frame of pvol
        // transform those corners and keep track of minimum and maximum extent
        // TODO: could make this code shorter with a more complex Vector3D class
        Vector3D<Precision> lower, upper;
        pvol->Extent( lower, upper );
        Vector3D<Precision> delta = upper-lower;
        Precision minx,miny,minz,maxx,maxy,maxz;
        minx = kInfinity;
        miny = kInfinity;
        minz = kInfinity;
        maxx = -kInfinity;
        maxy = -kInfinity;
        maxz = -kInfinity;
        Transformation3D const * transf = pvol->GetTransformation();
        for(int x=0;x<=1;++x)
            for(int y=0;y<=1;++y)
                for(int z=0;z<=1;++z){
                      Vector3D<Precision> corner;
                      corner.x() = lower.x() + x*delta.x();
                      corner.y() = lower.y() + y*delta.y();
                      corner.z() = lower.z() + z*delta.z();
                      Vector3D<Precision> transformedcorner =
                        transf->InverseTransform( corner );
                      minx = std::min(minx, transformedcorner.x());
                      miny = std::min(miny, transformedcorner.y());
                      minz = std::min(minz, transformedcorner.z());
                      maxx = std::max(maxx, transformedcorner.x());
                      maxy = std::max(maxy, transformedcorner.y());
                      maxz = std::max(maxz, transformedcorner.z());
                }
        *lowerc = Vector3D<Precision>(minx - 1E-3 ,miny - 1E-3, minz - 1E-3);
        *upperc = Vector3D<Precision>(maxx + 1E-3 ,maxy +  1E-3,maxz + 1E-3 );

#ifdef CHECK
        // do some tests on this stuff
        delta = (*upperc - *lowerc)/2.;
        Vector3D<Precision> boxtranslation = (*lowerc + *upperc)/2.;
        UnplacedBox box(delta);
        Transformation3D tr( boxtranslation.x(), boxtranslation.y(), boxtranslation.z() );
        VPlacedVolume const * boxplaced = LogicalVolume("",&box).Place(&tr);
        // no point on the surface of the aligned box should be inside the volume
        std::cerr << "lower " << *lowerc;
        std::cerr << "upper " << *upperc;
        int contains = 0;
        for(int i=0;i<10000;++i)
        {
            Vector3D<Precision> p =  box.GetPointOnSurface() + boxtranslation;
            std::cerr << *lowerc << " " << * upperc << " " << p << "\n";
            if( pvol->Contains( p ) ) contains++;
        }
        if( contains > 10){
            Visualizer visualizer;
            visualizer.AddVolume(*pvol, *pvol->GetTransformation());
            visualizer.AddVolume(*boxplaced, tr );
            visualizer.Show();
        }
        std::cerr << "## wrong points " << contains << "\n";
#endif
}

void ABBoxManager::InitABBoxes( LogicalVolume const * lvol ){
        if( fVolToABBoxesMap.find(lvol) != fVolToABBoxesMap.end() )
        {
            // remove old boxes first
            RemoveABBoxes(lvol);
        }
        int ndaughters = lvol->daughtersp()->size();
        // is this insertion correct ?
        ABBox_t * boxes = new ABBox_t[ 2*ndaughters ];
        fVolToABBoxesMap[lvol] = boxes;

        Visualizer visualizer;
        // calculate boxes by iterating over daughters
        for(int d=0;d<ndaughters;++d){
            auto pvol = lvol->daughtersp()->operator [](d);
            ComputeABBox( pvol, &boxes[2*d], &boxes[2*d+1] );
#ifdef CHECK
            // do some tests on this stuff
            Vector3D<Precision> lower = boxes[2*d];
            Vector3D<Precision> upper = boxes[2*d+1];

            Vector3D<Precision> delta = (upper - lower)/2.;
           Vector3D<Precision> boxtranslation = (lower + upper)/2.;
           UnplacedBox box(delta);
           Transformation3D tr( boxtranslation.x(), boxtranslation.y(), boxtranslation.z() );
           VPlacedVolume const * boxplaced = LogicalVolume("",&box).Place(&tr);
//                   int contains = 0;
//                   for(int i=0;i<10000;++i)
//                   {
//                       Vector3D<Precision> p =  box.GetPointOnSurface() + boxtranslation;
//                       std::cerr << *lowerc << " " << * upperc << " " << p << "\n";
//                       if( pvol->Contains( p ) ) contains++;
//                   }
//                   if( contains > 10){
            visualizer.AddVolume(*pvol, *pvol->GetTransformation());
            visualizer.AddVolume(*boxplaced, tr );
#endif
        }
     //    visualizer.Show();
}



void ABBoxManager::RemoveABBoxes( LogicalVolume const * lvol){
        if( fVolToABBoxesMap.find(lvol) != fVolToABBoxesMap.end() ) {
            delete[] fVolToABBoxesMap[lvol];
            fVolToABBoxesMap.erase(lvol);
        }
    }
}




}
