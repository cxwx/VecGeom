#include "UMultiUnion.hh"

#include <iostream>
#include "UUtils.hh"
#include "UVoxelFinder.hh"
#include <sstream>

using namespace std;

//______________________________________________________________________________       
UMultiUnion::UMultiUnion(const char *name)
{
   SetName(name);
   fNodes  = 0; 
   fVoxels = 0;  
}

//______________________________________________________________________________       
UMultiUnion::~UMultiUnion()
{
   if(fNodes) delete fNodes;
   if(fVoxels) delete fVoxels;
}

//______________________________________________________________________________       
void UMultiUnion::AddNode(VUSolid *solid, UTransform3D *trans)
{
   UNode* node = new UNode(solid,trans);   

   if (!fNodes)
   {
      fNodes = new vector<UNode*>;
   }
   
   fNodes->push_back(node);
}

//______________________________________________________________________________       
double UMultiUnion::Capacity()
{
   // Capacity computes the cubic volume of the "UMultiUnion" structure using random points

   // Random initialization:
   srand(time(NULL));

   double* extentMin = new double[3];
   double* extentMax = new double[3];
   UVector3 tempPoint;
   double dX, dY, dZ, oX, oY, oZ;
   int iGenerated = 0, iInside = 0;
   
   this->Extent(extentMin,extentMax);

   dX = (extentMax[0] - extentMin[0])/2;
   dY = (extentMax[1] - extentMin[1])/2;
   dZ = (extentMax[2] - extentMin[2])/2;
   
   oX = (extentMax[0] + extentMin[0])/2;
   oY = (extentMax[1] + extentMin[1])/2;
   oZ = (extentMax[2] + extentMin[2])/2;     
   
   double vbox = (2*dX)*(2*dY)*(2*dZ);
   
   while(iInside < 100000)
   {
      tempPoint.Set(oX - dX + 2*dX*(rand()/(double)RAND_MAX),
                    oY - dY + 2*dY*(rand()/(double)RAND_MAX),
                    oZ - dZ + 2*dZ*(rand()/(double)RAND_MAX));
      iGenerated++;
      if(this->Inside(tempPoint) != eOutside) iInside++;
   }
   double capacity = iInside*vbox/iGenerated;
   return capacity;      
}

//______________________________________________________________________________
void UMultiUnion::ComputeBBox (UBBox */*aBox*/, bool /*aStore*/)
{
// Computes bounding box.
   cout << "ComputeBBox - Not implemented" << endl;
}   

//______________________________________________________________________________
double UMultiUnion::DistanceToIn(const UVector3 &aPoint, 
                          const UVector3 &aDirection, 
                       // UVector3 &aNormal, 
                          double aPstep) const
{
// Computes distance from a point presumably outside the solid to the solid 
// surface. Ignores first surface if the point is actually inside. Early return
// infinity in case the safety to any surface is found greater than the proposed
// step aPstep.
// The normal vector to the crossed surface is filled only in case the box is 
// crossed, otherwise aNormal.IsNull() is true.
   cout << "DistanceToIn - Not implemented" << endl;
   return 0.;
}     

//______________________________________________________________________________
double UMultiUnion::DistanceToOut(const UVector3 &aPoint, const UVector3 &aDirection,
			       UVector3 &aNormal,
			       bool     &convex,
                double   aPstep) const
{
// Computes distance from a point presumably intside the solid to the solid 
// surface. Ignores first surface along each axis systematically (for points
// inside or outside. Early returns zero in case the second surface is behind
// the starting point.
// o The proposed step is ignored.
// o The normal vector to the crossed surface is always filled.

   cout << "DistanceToOut - Not implemented" << endl;
   return 0.;
}  

//______________________________________________________________________________
VUSolid::EnumInside UMultiUnion::Inside(const UVector3 &aPoint) const
{
// Classify point location with respect to solid:
//  o eInside       - inside the solid
//  o eSurface      - close to surface within tolerance
//  o eOutside      - outside the solid

// Hitherto, it is considered that:
//        - only parallelepipedic nodes can be added to the container

   // Implementation using voxelisation techniques:
   // ---------------------------------------------
   int iIndex;
   vector<int> vectorOutcome;   
   VUSolid *tempSolid = 0;
   UTransform3D *tempTransform = 0;
   UVector3 tempPoint, tempPointConv;
   VUSolid::EnumInside tempInside = eOutside;
   int countSurface = 0;//, countInside = 0, countOutside = 0;
         
   vectorOutcome = fVoxels -> GetCandidatesVoxelArray(aPoint); 

   for(iIndex = 0 ; iIndex < (int)vectorOutcome.size() ; iIndex++)
   {
      tempSolid = ((*fNodes)[vectorOutcome[iIndex]])->fSolid;
      tempTransform = ((*fNodes)[vectorOutcome[iIndex]])->fTransform;  
            
      // The coordinates of the point are modified so as to fit the intrinsic solid local frame:
      tempPoint.Set(aPoint.x,aPoint.y,aPoint.z);   
      tempPointConv = tempTransform->LocalPoint(tempPoint);
     
      tempInside = tempSolid->Inside(tempPointConv);        
      
      if(tempInside == eSurface) countSurface++; 
      
      if(tempInside == eInside) return eInside;      
   }       
   if(countSurface != 0) return eSurface;
   return eOutside;   
}

//______________________________________________________________________________ 
void UMultiUnion::Extent(EAxisType aAxis,double &aMin,double &aMax)
{
// Determines the bounding box for the considered instance of "UMultipleUnion"
   int carNodes = fNodes->size();  
   int iIndex = 0;
   double mini = 0, maxi = 0;  

   VUSolid *tempSolid = 0;
   UTransform3D *tempTransform = 0;
  
   double* arrMin = new double[3];
   double* arrMax = new double[3];
   
   UVector3 min, max;
    
   // Loop on the nodes:
   for(iIndex = 0 ; iIndex < carNodes ; iIndex++)
   {
      tempSolid = ((*fNodes)[iIndex])->fSolid;
      tempTransform = ((*fNodes)[iIndex])->fTransform;
      
      tempSolid->Extent(arrMin, arrMax);
      min.Set(arrMin[0],arrMin[1],arrMin[2]);      
      max.Set(arrMax[0],arrMax[1],arrMax[2]);           
      UUtils::TransformLimits(min, max, tempTransform);
     
      if(iIndex == 0)
      {
         switch(aAxis)
         {
            case eXaxis:
            {
               mini = min.x;
               maxi = max.x;            
               break;
            }
            case eYaxis:
            {
               mini = min.y;
               maxi = max.y;            
               break;
            }
            case eZaxis:
            {
               mini = min.z;
               maxi = max.z;            
               break;
            }                       
         }
         continue;
      }            
     
      // Deternine the min/max on the considered axis:
      switch(aAxis)
      {
         case eXaxis:
         {
            if(min.x < mini)
               mini = min.x;
            if(max.x > maxi)
               maxi = max.x;
            break;
         }
         case eYaxis:
         {
            if(min.y < mini)
               mini = min.y;
            if(max.y > maxi)
               maxi = max.y;
            break;
         }
         case eZaxis:
         {
            if(min.z < mini)
               mini = min.z;
            if(max.z > maxi)
               maxi = max.z;
            break;
         }                        
      }                 
   }
   aMin = mini;
   aMax = maxi;
}

//______________________________________________________________________________ 
void UMultiUnion::Extent(double aMin[3],double aMax[3])
{
   double min = 0,max = 0;
   this->Extent(eXaxis,min,max);
   aMin[0] = min; aMax[0] = max;
   this->Extent(eYaxis,min,max);
   aMin[1] = min; aMax[1] = max;
   this->Extent(eZaxis,min,max);
   aMin[2] = min; aMax[2] = max;      
}

//______________________________________________________________________________
bool UMultiUnion::Normal(const UVector3& aPoint, UVector3 &aNormal)
{
// Computes the normal on a surface and returns it as a unit vector
//   In case a point is further than tolerance_normal from a surface, set validNormal=false
//   Must return a valid vector. (even if the point is not on the surface.)
//
//   On an edge or corner, provide an average normal of all facets within tolerance
// NOTE: the tolerance value used in here is not yet the global surface
//     tolerance - we will have to revise this value - TODO

   int iIndex;
   vector<int> vectorOutcome;   
   VUSolid *tempSolid = 0;
   UTransform3D *tempTransform = 0;
   UVector3 tempPoint, tempPointConv, outcomeNormal;
   UVector3 temp1, temp2, temp3;     
       
   vectorOutcome = fVoxels -> GetCandidatesVoxelArray(aPoint); 

   if(this->Inside(aPoint) != eSurface)
   {
      return false;
   }
   else
   {
      for(iIndex = 0 ; iIndex < (int)vectorOutcome.size() ; iIndex++)
      {
         tempSolid = ((*fNodes)[vectorOutcome[iIndex]])->fSolid;
         tempTransform = ((*fNodes)[vectorOutcome[iIndex]])->fTransform;  
               
         // The coordinates of the point are modified so as to fit the intrinsic solid local frame:
         tempPoint.Set(aPoint.x,aPoint.y,aPoint.z);
         tempPointConv = tempTransform->LocalPoint(tempPoint);
         
         if(tempSolid->Inside(tempPointConv) == eSurface)
         {
            if(tempSolid->Normal(tempPointConv,outcomeNormal) != true) continue;
            
            temp1 = (tempTransform->GlobalPoint(outcomeNormal));
            temp2 = (tempTransform->GlobalPoint(tempPointConv));    
            temp3.Set(-temp1.x+temp2.x,-temp1.y+temp2.y,-temp1.z+temp2.z)  ;   
              
            aNormal = temp3.Unit();          
            return true;  
         }
      } 
      return false;
   }
}

//______________________________________________________________________________ 
double UMultiUnion::SafetyFromInside(const UVector3 aPoint, bool aAccurate) const
{
   // Estimates isotropic distance to the surface of the solid. This must
   // be either accurate or an underestimate. 
   //  Two modes: - default/fast mode, sacrificing accuracy for speed
   //             - "precise" mode,  requests accurate value if available.   

//   int iIndex;
   vector<int> vectorOutcome;
   int numberCandidates;
     
   VUSolid *tempSolid = 0;
   UTransform3D *tempTransform = 0;
   
   UVector3 tempPointConv;      
         
   vectorOutcome = fVoxels -> GetCandidatesVoxelArray(aPoint); 
   
   // Determination of the number of candidates detected for the considered point:
   numberCandidates = (int)vectorOutcome.size();

   // If only one candidate is detected, a simple algorithm of "SafetyFromInside" is applied:
   if(numberCandidates == 1)
   {
      tempSolid = ((*fNodes)[vectorOutcome[0]])->fSolid;
      tempTransform = ((*fNodes)[vectorOutcome[0]])->fTransform;
           
      tempPointConv = tempTransform->LocalPoint(aPoint);
      
      return tempSolid->SafetyFromInside(tempPointConv,aAccurate);
   }
   // If the number of candidates is greater than one, there is overlapping (or at least contact):
   else if(numberCandidates > 1)
   {
      // Boolean union - TODO
      return 0.;
   }
   return 0.;   
}

//______________________________________________________________________________
double UMultiUnion::SafetyFromOutside(const UVector3 aPoint, bool aAccurate) const
{
   // Estimates the isotropic safety from a point outside the current solid to any 
   // of its surfaces. The algorithm may be accurate or should provide a fast 
   // underestimate.
   
   cout << "SafetyFromOutside - Not implemented" << endl;
   return 0.;
}

//______________________________________________________________________________       
double UMultiUnion::SurfaceArea()
{
// Computes analytically the surface area.
   cout << "SurfaceArea - Not implemented" << endl;
   return 0.;
}     

//______________________________________________________________________________       
int UMultiUnion::GetNumNodes() const
{
   return(fNodes->size());
}

//______________________________________________________________________________       
/*const*/ VUSolid* UMultiUnion::GetSolid(int index) /*const*/
{
   return(((*fNodes)[index])->fSolid);
}

//______________________________________________________________________________       
const UTransform3D* UMultiUnion::GetTransform(int index) const
{
   return(((*fNodes)[index])->fTransform);
}

//______________________________________________________________________________       
void UMultiUnion::SetVoxelFinder(UVoxelFinder* finder)
{
   fVoxels = finder;
}

//______________________________________________________________________________       
void UMultiUnion::Voxelize()
{
   fVoxels = new UVoxelFinder(this);
   fVoxels -> Voxelize();
}
