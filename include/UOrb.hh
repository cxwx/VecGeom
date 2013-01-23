/**
* @file
* @author  John Doe <jdoe@example.com>
* @version 0.1
*
* @section LICENSE
*
* @section DESCRIPTION
*
* A simple Orb defined by half-lengths on the three axis. The center of the Orb matches the origin of the local reference frame.
*/

#ifndef USOLIDS_UOrb
#define USOLIDS_UOrb

#include "VUSolid.hh"
#include "UUtils.hh"

class UOrb : public VUSolid
{

public:
	UOrb() : VUSolid(), fR(0), fRTolerance(0) {}
	UOrb(const std::string &name, double pRmax);
	~UOrb() {}

	// Navigation methods
	EnumInside     Inside (const UVector3 &aPo6int) const;   

	double  SafetyFromInside ( const UVector3 &aPoint, 
		bool aAccurate=false) const;
	double  SafetyFromOutside( const UVector3 &aPoint, 
		bool aAccurate=false) const;
	double  DistanceToIn     ( const UVector3 &aPoint, 
		const UVector3 &aDirection,
		// UVector3       &aNormalVector,
		double aPstep = UUtils::kInfinity) const;                               

	double DistanceToOut     ( const UVector3 &aPoint,
		const UVector3 &aDirection,
		UVector3       &aNormalVector, 
		bool           &aConvex,
		double aPstep = UUtils::kInfinity) const;

	bool Normal ( const UVector3& aPoint, UVector3 &aNormal ) const; 
//	virtual void Extent ( EAxisType aAxis, double &aMin, double &aMax ) const;
	void Extent (UVector3 &aMin, UVector3 &aMax) const; 
	double Capacity();
	double SurfaceArea();
	UGeometryType GetEntityType() const { return "Orb";}
	void ComputeBBox(UBBox * /*aBox*/, bool /*aStore = false*/) {}

	//G4Visualisation
	void GetParametersList(int /*aNumber*/, double * /*aArray*/) const{} 
	UPolyhedron* GetPolyhedron() const{return CreatePolyhedron(); }

	inline VUSolid* Clone() const
	{
		return new UOrb(GetName(), fR);
	}

	double GetRadialTolerance()
	{
		return fRTolerance;
	}

	UVector3 GetPointOnSurface() const;

	std::ostream& StreamInfo(std::ostream& os) const;

	UPolyhedron* CreatePolyhedron () const;

private:  
	double fR;
	double fRTolerance;

	double DistanceToOutForOutsidePoints(const UVector3 &p, const UVector3 &v, UVector3 &n) const;

};
#endif
