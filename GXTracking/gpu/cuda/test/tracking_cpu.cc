#include <iostream>
#include <fstream>

#include <cuda.h>
#include <curand.h>
#include "stopwatch.h"

#include <stdio.h>
#include <math.h>

#include "tracking_kernel.h"
#include "random_kernel.h"

#include "GPConstants.h"
#include "GXFieldMap.h"
#include "GXFieldMapData.h"

#include "GXTrack.h"
#include "GPThreeVector.h"

//Geometry
#include "GPVGeometry.h"
#include "GPUserGeometry.h"
#include "GPSimpleEcal.h"
#include "GPSimpleCMS.h"
#include "GPVPhysicalVolume.h"

//EMPhysics
#include "GPConstants.h"
#include "GXEMPhysicsUtils.h"

#include "tracking_kernel.cu"

using namespace std;

int main (int argc, char* argv[]) 
{
  //argument
  int kernelType = 1; //1. general traportation kernel
  //2. electron/photon kernels
  int flagTest = 0; //1 test mode
  int ndivarg = -1;

  if(argc >= 2) kernelType = atoi(argv[1]);
  if(argc >= 3) flagTest = atoi(argv[2]);
  if(argc >= 4) ndivarg = atoi(argv[3]);

  if(kernelType >2 ) {
    std::cout << "Wrong kernelType: Bailing out ... "<< std::endl;
    return 0;
  }
  
  //check whether device overlap is available
  
  //default number of blocks and threads
  
  int theNBlocks  =  32;
  int theNThreads = 128;
  //int theNBlocks  = 128;
  //int theNThreads = 448;
  
  //reconfigure 
  if(flagTest>0) {
    theNBlocks = 2;
    theNThreads = 10;
  }
  
  char* cudaNBlocks = getenv("GP_CUDA_NBLOCKS");
  char* cudaNThreads = getenv("GP_CUDA_NTHREADS");
  
  if(cudaNBlocks) theNBlocks = atoi(cudaNBlocks);
  if(cudaNThreads) theNThreads = atoi(cudaNThreads);
  
  std::cout << "... tracking_kernel<<<" << theNBlocks << "," 
  << theNThreads <<">>> (...) ..." << std::endl;
  
  //1. Construct geometry
  int nphi = 4;
  int nz   = 3;
  double density = 8.28;
  
  GPVGeometry *geom = new GPSimpleEcal(nphi,nz,density);
  //  GPVGeometry *geom = new GPSimpleCMS();
  
  geom->create();
  
  GPVGeometry::byte *geom_h = (GPVGeometry::byte *) malloc (geom->size()) ;
  geom->relocate( geom_h );
  memcpy(geom_h,geom->getBuffer(),geom->size());  
  
  
  //2. Read magnetic field map
  
  GXFieldMap *bmap_d;
  GXFieldMap *bmap_h;
  GXFieldMap** fieldMap;
  
  fieldMap = (GXFieldMap **) malloc (nbinZ*sizeof(GXFieldMap *));
  for (int j = 0 ; j < nbinZ ; j++) {
    fieldMap[j] = (GXFieldMap *) malloc (nbinR*sizeof(GXFieldMap));
  } 
  
  const char* fieldMapFile = getenv("GP_BFIELD_MAP");
  fieldMapFile = (fieldMapFile) ? fieldMapFile : "data/cmsExp.mag.3_8T";
  
  std::ifstream ifile(fieldMapFile, ios::in | ios::binary | ios::ate);
  
  if (ifile.is_open()) {
    
    //field map structure
    GXFieldMapData fd;
    
    ifstream::pos_type fsize = ifile.tellg();
    size_t dsize = sizeof(GXFieldMapData);    
    
    long int ngrid = fsize/dsize;
    ifile.seekg (0, ios::beg);
    
    std::cout << "... transportation ... Loading magnetic field map: " 
    << fieldMapFile << std::endl;
    
    for(int i = 0 ; i < ngrid ; i++) {
      ifile.read((char *)&fd, sizeof(GXFieldMapData));
      
      //check validity of input data
      if(abs(fd.iz) > noffZ || fd.ir > nbinR) {
        std::cout << " Field Map Array Out of Range" << std::endl;
      }
      else {
        fieldMap[noffZ+fd.iz][fd.ir].Bz = fd.Bz; 
        fieldMap[noffZ+fd.iz][fd.ir].Br = fd.Br;
      }
    }
    ifile.close();
  }
  
  //3. Create magnetic field on the device
  
  printf("Creating magnetic field map on the GPU device\n");
  //prepare fieldMap array: fieldMap[nbinZ][nbinR];
  
  bmap_h = (GXFieldMap *) malloc (nbinZ*nbinR*sizeof(GXFieldMap));
  
  for (int i = 0 ; i < nbinZ ; i++) {
    for (int j = 0 ; j < nbinR ; j++) {
      bmap_h[i+j*nbinZ].Bz = fieldMap[i][j].Bz;
      bmap_h[i+j*nbinZ].Br = fieldMap[i][j].Br;
    }
  }
  
  printf("Copying data from host to device\n");
  
  
  // 4. Prepare EM physics tables
  bool useSpline = true;
  
  GPPhysicsTable eBrem_table;
  readTable(&eBrem_table,"data/Lambda.eBrem.e-.asc");
  int nv = eBrem_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eBrem_table.physicsVectors[j].SetSpline(useSpline);
  }
  //  GPPhysicsTable_Print(&eBrem_table);
  
  GPPhysicsTable eIoni_table;
  readTable(&eIoni_table,"data/Lambda.eIoni.e-.asc");
  nv = eIoni_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eIoni_table.physicsVectors[j].SetSpline(useSpline);
  }
  //  GPPhysicsTable_Print(&eIoni_table);
  
  GPPhysicsTable msc_table;
  readTable(&msc_table,"data/Lambda.msc.e-.asc");
  nv = msc_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    msc_table.physicsVectors[j].SetSpline(useSpline);
  }
  //  GPPhysicsTable_Print(&msc_table);
  
  printf("Copying GPPhysicsTable data from host to device\n");
  
  GPPhysicsTable* eBrem_table_d;
  GPPhysicsTable* eIoni_table_d;
  GPPhysicsTable* msc_table_d;
  
  
  //5. Get the list of tracks from the stack
  printf("Preparing a bundle tracks at host\n");
  
  //event-loop
  int nevent = 10;
  
  //initialize the stream
  
  
  int nTracks = 4096*24;
  //int nTracks = 100000;
  GXTrack *track_h = 0;
  GXTrack *track_c = 0;
  
  track_h = new GXTrack[nTracks];
  track_c = new GXTrack[nTracks];
  
  GXTrack *electron_h = 0;
  GXTrack *electron_c = 0;


  GXTrack *photon_h = 0;
  GXTrack *photon_c = 0;

  
  if (kernelType == 2) {
    int nPhotons = nTracks;
    int nElectrons = nTracks;

    electron_h = new GXTrack[nElectrons];
    photon_h = new GXTrack[nPhotons];
    electron_c = new GXTrack[nElectrons];
    photon_c = new GXTrack[nPhotons]
;
  }
  double *results = 0;
  
  for (int i = 0 ; i < nevent ; i++) {
    
    int nDiv = ndivarg > 0 ? ndivarg : 24;
    
    int Nchunk  =       0;
    int NchunkG  =      0;
    int NchunkE  =      0;
    
    if(flagTest>0) {
      nTracks = 20;
    }
    
    printf("Event %d> Transportation in kernel for %d tracks\n",i,nTracks);
    
    //populate tracks with track position, momentum, energy, steplength
    
    int nPhotons = 0;    
    int nElectrons = 0;    
    
    //photon probability: default 20%
    G4double photonFraction = 0.2;
    
    //randomize position and momentum    
    for(size_t i = 0 ; i < nTracks ; i++){
      //barrel only
      track_h[i].x      = track_c[i].x     = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].y      = track_c[i].y     = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].z      = track_c[i].z     = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].s      = track_c[i].s     = 1.0+1*(2.0*rand()/RAND_MAX-1.0);
      
      track_h[i].px     = track_c[i].px    = 1.+1000*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].py     = track_c[i].py    = 1.+1000*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].pz     = track_c[i].pz    = 1.+1000*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].q      = track_c[i].q     = -1.0;
      
      if( 1.0*rand()/(float)RAND_MAX < photonFraction) {
        track_h[i].q      = track_c[i].q     = 0.0;
        nPhotons++;
      }
      else {
        track_h[i].q      = track_c[i].q     = -1.0;
        nElectrons++;
      }
    }
    
    if(kernelType ==2) {
      int j = 0;
      int k = 0;
      
      for(size_t i = 0 ; i < nTracks ; i++){
        if(track_h[i].q == 0.0) {
          photon_h[j].x      = photon_c[j].x     = track_h[i].x  ;
          photon_h[j].y      = photon_c[j].y     = track_h[i].y  ;
          photon_h[j].z      = photon_c[j].z     = track_h[i].z  ;
          photon_h[j].s      = photon_c[j].s     = track_h[i].s  ;
          
          photon_h[j].px     = photon_c[j].px    = track_h[i].px ;
          photon_h[j].py     = photon_c[j].py    = track_h[i].py ;
          photon_h[j].pz     = photon_c[j].pz    = track_h[i].pz ;
          photon_h[j].q      = photon_c[j].q     = track_h[i].q ;
          j++;  
        }
        else {
          electron_h[k].x      = electron_c[k].x     = track_h[i].x  ;
          electron_h[k].y      = electron_c[k].y     = track_h[i].y  ;
          electron_h[k].z      = electron_c[k].z     = track_h[i].z  ;
          electron_h[k].s      = electron_c[k].s     = track_h[i].s  ;
          
          electron_h[k].px     = electron_c[k].px    = track_h[i].px ;
          electron_h[k].py     = electron_c[k].py    = track_h[i].py ;
          electron_h[k].pz     = electron_c[k].pz    = track_h[i].pz ;
          electron_h[k].q      = electron_c[k].q     = track_h[i].q ;
          k++;
        }
      }
    }
    
    printf(" (nElectron,nPhoton) = (%d,%d)\n",nElectrons,nPhotons);
    
    float elapsedTimeAlloc = 0.0;
    float elapsedTimeGPU = 0.0;
    
    printf("Copying track data from host to device\n");
    
    //start time memory allocation on GPU
    
    GXTrack *track_d0;
    GXTrack *track_d1;
    
    GXTrack *photon_d0;
    GXTrack *photon_d1;
    
    GXTrack *electron_d0;
    GXTrack *electron_d1;
    
    if(kernelType ==1) {
       // Nchunk = ( nTracks%nDiv != 0 ) ? (nTracks/nDiv + 1 ) : (nTracks/nDiv);
       int nthreads_total = theNBlocks * theNThreads;
       Nchunk = (nTracks/nDiv);
       if (Nchunk >= nthreads_total) {
          Nchunk -= ( (nTracks/nDiv) % nthreads_total); // align on multiple of the grid size
       } else {
          fprintf(stderr,"Warning: too many division (%d) resulting in low chunk size %d (< %d)\n",nDiv,Nchunk,nthreads_total);
       }
       
    } 
    else if(kernelType ==2) {
       
       int nthreads_total = theNBlocks * theNThreads;
       
       NchunkG = (nPhotons/nDiv);
       if (NchunkG >= nthreads_total) {
          NchunkG -= ( (nPhotons/nDiv) % nthreads_total); // align on multiple of the grid size
       } else {
          fprintf(stderr,"Warning: too many division for photons (%d) resulting in low chunk size %d (< %d)\n",nDiv,NchunkG,nthreads_total);
       }

       NchunkE = (nElectrons/nDiv);
       if (NchunkE >= nthreads_total) {
          NchunkE -= ( (nElectrons/nDiv) % nthreads_total); // align on multiple of the grid size
       } else {
          fprintf(stderr,"Warning: too many division for electrons (%d) resulting in low chunk size %d (< %d)\n",nDiv,NchunkE,nthreads_total);
       }
       
    }
    
    //start steam loop
    
    //prepare random engines on the device
        
    StopWatch timer;
    timer.start();
    
    //<<<---------------------------------------------------------->>>
    if(kernelType ==1) {
      tracking_cpu((GPGeomManager*)geom_h,bmap_h,track_c,
                   &eBrem_table,&eIoni_table,&msc_table,
                   nTracks);
    }
    else if(kernelType ==2 ) {
      tracking_photon_cpu((GPGeomManager*)geom_h,bmap_h,
                          photon_c,
                          &eBrem_table,&eIoni_table,&msc_table,
                          nPhotons);
      tracking_electron_cpu((GPGeomManager*)geom_h,bmap_h,
                            electron_c,
                            &eBrem_table,&eIoni_table,&msc_table,
                            nElectrons);
    }
    //<<<---------------------------------------------------------->>>
    
    timer.stop();
    // cudaEventRecord (stop,0);
    // cudaEventSynchronize (stop);
    
    float elapsedTime2 = timer.getTime();;
    // cudaEventElapsedTime (&elapsedTime2,start,stop);
    printf("Time Elapsed on CPU : %6.3f ms\n",elapsedTime2);
    
    //Compare results from GPU and CPU
    /*
     long int miscount = 0;
     for (size_t i = 0 ; i < nTracks ; i++) {
     if( track_h[i].x != track_c[i].x || track_h[i].y != track_c[i].y || 
     track_h[i].z != track_c[i].z ) {
     std::cout << "Compare results: check track position from gpu cpu " 
     << "(" << track_h[i].x << "," << track_h[i].y << "," 
     <<track_h[i].z << ")" << " (" << track_c[i].x << "," 
     << track_c[i].y <<","<<track_c[i].z << ")" << std::endl;
     miscount++;
     }
     }
     
     printf("Number of Cuda MisMatched Tracks = %ld\n",miscount);
     */
    //clean up: destory cuda event and free memory on device and host
  }
  delete track_h;
  delete electron_h;
  delete photon_h;

  delete track_c;
  delete electron_c;
  delete photon_c;
}