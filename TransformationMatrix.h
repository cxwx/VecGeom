/*
 * TransformationMatrix.h
 *
 *  Created on: Nov 8, 2013
 *      Author: swenzel
 */

#ifndef TRANSFORMATIONMATRIX_H_
#define TRANSFORMATIONMATRIX_H_

#include "Vector3D.h"
#include <cmath>

typedef int TranslationIdType;
typedef int RotationIdType;

class TransformationMatrix
{
private:
    // this can be varied depending on template specialization
	double trans[3];
	double rot[9];

	template<RotationIdType rid=-1>
	inline
	void
	emitrotationcode(Vector3D const &, Vector3D &) const;


public:
	// constructor
	TransformationMatrix(double const *t, double const *r)
	{
	    trans[0]=t[0];
	    trans[1]=t[1];
		trans[2]=t[2];
		for(int i=0;i<9;i++)
			rot[i]=r[i];
		// we need to check more stuff ( for instance that product along diagonal is +1)
	}

	inline
	static TranslationIdType GetTranslationIdType(double const *t)
	{
		if( t[0]==0 && t[1]==0 && t[2]==0 )
			return 0;
		else
			return 1;
	}

	inline
	static RotationIdType GetRotationIdType(double const *r)
	{
		// this function has to be generated by some script
		if( r[0]==1 && r[1]==0 && r[2]==0 && r[3]==0 && r[4]==1 && r[5]==0 && r[6]==0 && r[7]==0 && r[8]==1 )
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}


	template <TranslationIdType tid=-1, RotationIdType rid=-1>
	inline
	void
	MasterToLocal(Vector3D const &, Vector3D &) const;

	// to transform real vectors, we don't care about translation
	template <RotationIdType rid=-1>
	inline
	void
	MasterToLocalVec(Vector3D const &, Vector3D &) const;

	// why not ??
	// inline Vector3D MasterToLocal(Vector3D const &) const;

	// for a Vc vector: not so easy
	// inline MasterToLocal(Vc::double) const;


	friend class PhysicalVolume;
	friend class PhysicalBox;
};


template <TranslationIdType tid=-1, RotationIdType rid=-1>
void
TransformationMatrix::MasterToLocal(Vector3D const & master, Vector3D & local) const
{
	if( tid==0 && rid ==0 ) // this means identity
	{
		local.x = master.x;
		local.y = master.y;
		local.z = master.z;
	}
	else if( tid != 0 && rid == 0 ) // tid == 1 means we have
	{
		local.x = master.x + trans[0];
		local.y = master.y + trans[1];
		local.z = master.z + trans[2];
	}
	else if( tid == 0 && rid!=0 ) // pure rotation
	{
		double mt0, mt1, mt2;
		mt0 = master.x;
		mt1 = master.y;
		mt2 = master.z;
		// try to inline code for rotations ( don't repeat yourself )
		local.x = mt0*rot[0] + mt1*rot[3] + mt2*rot[6];
		local.y = mt0*rot[1] + mt1*rot[4] + mt2*rot[7];
		local.z = mt0*rot[2] + mt1*rot[5] + mt2*rot[8];
	}
	else( tid != 0 && rid!=0 ) // both rotation and translation
	{
		double mt0, mt1, mt2;
		mt0 = master.x + trans[0];
		mt1 = master.y + trans[1];
		mt2 = master.z + trans[2];
		local.x = mt0*rot[0] + mt1*rot[3] + mt2*rot[6];
		local.y = mt0*rot[1] + mt1*rot[4] + mt2*rot[7];
		local.z = mt0*rot[2] + mt1*rot[5] + mt2*rot[8];
	}
}

template <RotationIdType rid=-1>
void
TransformationMatrix::MasterToLocalVec(Vector3D const & master, Vector3D & local ) const
{
	MasterToLocal<0, rid>(master, local);
}



#endif /* TRANSFORMATIONMATRIX_H_ */
