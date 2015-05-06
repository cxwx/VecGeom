
#include "LogicalVolumeHandler.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4TransportationManager.hh"

//#include <cstdio>

LogicalVolumeHandler* LogicalVolumeHandler::fgInstance = 0;

LogicalVolumeHandler* LogicalVolumeHandler::Instance() {
  if(!fgInstance)
    fgInstance = new LogicalVolumeHandler();

 return fgInstance;
}


LogicalVolumeHandler::LogicalVolumeHandler() {   
   Initialize();
}

void  LogicalVolumeHandler::Initialize(){
 std::cout<< "####  Checking logical volumes that are used in the detector .... "<<std::endl;
 G4LogicalVolume *lworld = G4TransportationManager::GetTransportationManager()
                           -> GetNavigatorForTracking() 
                           -> GetWorldVolume()
                           -> GetLogicalVolume();  // world logical volume pointer
 G4PhysicalVolumeStore *pvolstore = G4PhysicalVolumeStore::GetInstance();
 G4LogicalVolumeStore  *lvolstore = G4LogicalVolumeStore::GetInstance();
 G4int thesize = lvolstore->size();  // number of logical volumes in the store 

 // create a mask vector with size equal to the number of G4LogicalVolume-s in the store 
 // and initialise each elements to FALSE
 fIsLogicalVolumeUsedMask.resize(0);
 fIsLogicalVolumeUsedMask.resize(thesize,FALSE);

 // set a map from G4LogicalVolume* -> its index in the volume store
 fLogicVolToIndex.clear();
 for(G4int i = 0; i< thesize; ++i)
  fLogicVolToIndex[(*lvolstore)[i]] = i;


 // write this into file 
 FILE* f;
 f=fopen("/home/agheata/geant/bin/usedlogvolmask","r");
 for(unsigned int i = 0; i < fIsLogicalVolumeUsedMask.size(); ++i ){
    G4int the = 0;
    fscanf(f,"%d",&the);
    if(the) fIsLogicalVolumeUsedMask[i]=TRUE;
 }
 fclose(f);

/*
 // loop over the G4PhysicalVolume-s, get their G4LogicalVolume* and see if they are used in 
 // the geometry or not; set a mask value to TRUE if they are used 
 for(unsigned int ivol = 0; ivol < pvolstore->size(); ++ivol)
    if(lworld->IsAncestor((*pvolstore)[ivol]))
            fIsLogicalVolumeUsedMask[fLogicVolToIndex.find((*pvolstore)[ivol]->GetLogicalVolume())->second] = TRUE;

 // the world logical volume is always used 
 fIsLogicalVolumeUsedMask[fLogicVolToIndex.find(lworld)->second] = TRUE; // world is used 
*/

 // count number of used logical volumes 
 G4int numUsedLogicalVols = 0;
 for(unsigned int i = 0; i < fIsLogicalVolumeUsedMask.size(); ++i )
    if(fIsLogicalVolumeUsedMask[i]) ++numUsedLogicalVols;

/*
 // write this into file 
 FILE* f;
 f=fopen("usedlogvolmask","w");
 for(unsigned int i = 0; i < fIsLogicalVolumeUsedMask.size(); ++i ){
    if(i!=0 && i%100 == 0) fprintf(f,"\n");
    if(fIsLogicalVolumeUsedMask[i]) fprintf(f,"1 ");
    else fprintf(f,"0 ");
 }
 fclose(f);
*/

 std::cout<< " ####   Number of logical volumes used = "
          << numUsedLogicalVols << "  out of the " 
          << thesize << std::endl; 
}


G4int  LogicalVolumeHandler::GetIndex(G4LogicalVolume *lvolp){return fLogicVolToIndex.find(lvolp)->second;}

G4bool LogicalVolumeHandler::IsUsed(G4LogicalVolume *lvolp){ 
  return fIsLogicalVolumeUsedMask[fLogicVolToIndex.find(lvolp)->second];
}

G4bool LogicalVolumeHandler::IsUsed(G4int indx){ 
  return fIsLogicalVolumeUsedMask[indx];
}

