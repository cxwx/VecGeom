/*
 * volume_utilities.h
 *
 *  Created on: Mar 24, 2014
 *      Author: swenzel
 */

#ifndef VOLUME_UTILITIES_H_
#define VOLUME_UTILITIES_H_

#include "base/vector3d.h"
#include "base/global.h"
#include "base/rng.h"

namespace VECGEOM_NAMESPACE
{
	namespace volumeutilities
	{

		using namespace vecgeom;
		bool IsHittingVolume(Vector3D<Precision> const &point,
							 Vector3D<Precision> const &dir,
                             VPlacedVolume const &volume)
		{
			return volume.DistanceToIn(point, dir, kInfinity) < kInfinity;
		}

		static Vector3D<Precision> SamplePoint(Vector3D<Precision> const &size,
                                         const Precision scale = 1) {
			const Vector3D<Precision> ret(
					scale * (1. - 2. * RNG::Instance().uniform()) * size[0],
					scale * (1. - 2. * RNG::Instance().uniform()) * size[1],
					scale * (1. - 2. * RNG::Instance().uniform()) * size[2]
			);
			return ret;
		}


		static Vector3D<Precision> SampleDirection(){
			Vector3D<Precision> dir(
					(1. - 2. * RNG::Instance().uniform()),
					(1. - 2. * RNG::Instance().uniform()),
					(1. - 2. * RNG::Instance().uniform())
			);

			const Precision inverse_norm =
					1. / std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
			dir[0] *= inverse_norm;
			dir[1] *= inverse_norm;
			dir[2] *= inverse_norm;

			return dir;
		}


		template<typename TrackContainer>
		static void FillRandomDirections(TrackContainer & dirs){
			const int size = dirs.size();
			for (int i = 0; i < size; ++i) {
				const Vector3D<Precision> temp = SampleDirection();
				// should use push_back here; or dirs[i]=temps

				dirs.x(i) = temp[0];
				dirs.y(i) = temp[1];
				dirs.z(i) = temp[2];
			}
		}

		template<typename TrackContainer>
		static void FillBiasedDirections(VPlacedVolume const &volume,
                                   	   TrackContainer const &points,
                                   	   const Precision bias,
                                   	   TrackContainer & dirs) {
			assert(bias >= 0. && bias <= 1.);

			const int size = dirs.size();
			int n_hits = 0;
			std::vector<bool> hit(size, false);
			int h;

			// Randomize points
			FillRandomDirections(dirs);

			// Check hits
			for (int i = 0; i < size; ++i) {
				for (Iterator<Daughter> j = volume.daughters().begin();
						j != volume.daughters().end(); ++j) {
					if (IsHittingVolume(points[i], dirs[i], **j)) {
						n_hits++;
						hit[i] = true;
					}
				}
			}

			// Remove hits until threshold
			while (static_cast<Precision>(n_hits)/static_cast<Precision>(size) >= bias) {
				h = static_cast<int>(
						static_cast<Precision>(size) * RNG::Instance().uniform()
				);
				while (hit[h]) {
					dirs.Set(h, SampleDirection());
					for (Iterator<Daughter> i = volume.daughters().begin();
							i != volume.daughters().end(); ++i) {
						if (!IsHittingVolume(points[h], dirs[h], **i)) {
							n_hits--;
							hit[h] = false;
							break;
						}
					}
				}
			}


			// Add hits until threshold
			while (static_cast<Precision>(n_hits)/static_cast<Precision>(size) < bias) {
				h = static_cast<int>(
						static_cast<Precision>(size) * RNG::Instance().uniform()
				);
				while (!hit[h]) {
					dirs.Set(h, SampleDirection());
					for (Iterator<Daughter> i = volume.daughters().begin();
								i != volume.daughters().end(); ++i) {
						if (IsHittingVolume(points[h], dirs[h], **i)) {
							n_hits++;
							hit[h] = true;
							break;
						}
					}
				}
			}
		}


		template<typename TrackContainer>
		static void FillBiasedDirections(LogicalVolume const &volume,
                                   	   TrackContainer const &points,
                                   	   const Precision bias,
                                   	   TrackContainer & dirs)
		{
			VPlacedVolume const *const placed = volume.Place();
			FillBiasedDirections(*placed, points, bias, dirs);
			delete placed;
		}

		template<typename TrackContainer>
		static void FillUncontainedPoints(VPlacedVolume const &volume,
                                    	TrackContainer & points) {
			const int size = points.size();
			const Vector3D<Precision> dim = volume.bounding_box()->dimensions();
			for (int i = 0; i < size; ++i) {
				bool contained;
				do {
					points.Set(i, SamplePoint(dim));
					contained = false;
					for (Iterator<Daughter> j = volume.daughters().begin();
							j != volume.daughters().end(); ++j) {
						if ((*j)->Inside( points[i] )) {
							contained = true;
							break;
						}
					}
				} while (contained);
			}
		}


		template<typename TrackContainer>
		static void FillUncontainedPoints(LogicalVolume const &volume,
                                    	TrackContainer & points)
		{
			VPlacedVolume const *const placed = volume.Place();
			FillUncontainedPoints(*placed, points);
			delete placed;
		}

		template <typename TrackContainer>
		static void FillRandomPoints(VPlacedVolume const &volume,
															   TrackContainer &points) {
			const int size = points.size();
			const Vector3D<Precision> dim = volume.bounding_box()->dimensions();
			for (int i = 0; i < size; ++i) {
				Vector3D<Precision> point;
				do {
					point = SamplePoint(dim);
				} while (!volume.Inside(point));
				points.Set(i, point);
			}
		}

	} // end namespace volumeutilities
} // end namespace VEGGEOM

#endif /* VOLUME_UTILITIES_H_ */