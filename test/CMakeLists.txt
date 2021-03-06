################################################################################
## TESTING
################################################################################
# Helper Macros/Functions
macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE})
    target_link_libraries(${TARGET_NAME} ${VECGEOM_LIBRARIES} ${USERKERNELLIB})
    if(MIC)
      get_target_property(current_property ${TARGET_NAME} LINK_FLAGS)
      if(NOT current_property)
        set_target_properties(${TARGET_NAME} PROPERTIES LINK_FLAGS "-B/usr/linux-k1om-4.7/linux-k1om/usr/lib64/k1om-mpss-linux/4.7.4 -B/usr/linux-k1om-4.7/linux-k1om/lib64/4.7.4")
      else()
        set_target_properties(${TARGET_NAME} PROPERTIES LINK_FLAGS "${current_property} -B/usr/linux-k1om-4.7/linux-k1om/usr/lib64/k1om-mpss-linux/4.7.4 -B/usr/linux-k1om-4.7/linux-k1om/lib64/4.7.4")
      endif()
    endif()
  endforeach()
endmacro()

function(add_to_ctest EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    if(CTEST)
      get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
      add_test(${TARGET_NAME} ${TARGET_NAME})
    endif()
  endforeach()
endfunction()

# define a function downloading from a URL into a local file LOCALFILE
function(FILE_DOWNLOAD FILE_URL LOCALFILE )
  if(APPLE)
    execute_process(COMMAND curl -o ${LOCALFILE} ${FILE_URL})
  else()
    execute_process(COMMAND wget -q ${FILE_URL} -O ${LOCALFILE})
  endif()
endfunction()

# define a function checking md5 hashes
# result is stored in MD5MATCHES ( 1 == true, 0 == false )
function(CHECKMD5 FILETOCHECK EXPECTEDMD5HASH MD5MATCHES)
  if(APPLE)
    execute_process(COMMAND md5 ${FILETOCHECK} OUTPUT_VARIABLE MD5SUM)
    string(LENGTH ${MD5SUM} MD5LENGTH)
    MATH(EXPR START "${MD5LENGTH} - 33")
    string(SUBSTRING ${MD5SUM} ${START} 32 MD5SUM)
  else()
    execute_process(COMMAND md5sum ${FILETOCHECK} OUTPUT_VARIABLE MD5SUM)
    string(SUBSTRING ${MD5SUM} 0 32 MD5SUM)
  endif()
  if(MD5SUM STREQUAL EXPECTEDMD5HASH)
    set(${MD5MATCHES} 1 PARENT_SCOPE)
  else()
    set(${MD5MATCHES} 0 PARENT_SCOPE)
  endif()
endfunction()

# actual function for managing the download
function(DOWNLOAD_IF_NOT_INSTALLED FILE_URL LOCALFILE TARGETPATH MD5HASH )
  find_file(FOUNDFILE ${LOCALFILE} ${TARGETPATH} )
  if(FOUNDFILE STREQUAL "FOUNDFILE-NOTFOUND")
    message(STATUS "need download of ${LOCALFILE} since not found")
    set(NEEDTODOWNLOAD 1)
  else()
    message(STATUS "found existing file ${LOCALFILE}")
    CHECKMD5( ${FOUNDFILE} ${MD5HASH} MD5CORRECT )
    if(${MD5CORRECT} STREQUAL "1")
      set(NEEDTODOWNLOAD 0)
    else( )
      message(STATUS "hash ${MD5HASH} not correct for file ${FOUNDFILE} ${MD5CORRECT}" )
      set(NEEDTODOWNLOAD 1)
    endif()
  endif()

  if(${NEEDTODOWNLOAD} STREQUAL 1)
    message(STATUS " downloading ... ")
    set(DOWNLOADLOCATION "${TARGETPATH}/${LOCALFILE}")
    FILE_DOWNLOAD( ${FILE_URL} ${DOWNLOADLOCATION} )
  endif()
  # in principle have to check now if download succeeded and has right MD5
  # TOBEDONE

  # this is annoying but we have to clear FOUNDFILE SINCE THIS IS TREATED LIKE A STATIC VARIABLE
  unset(FOUNDFILE CACHE)
endfunction()
################################################################################

################################################################################
# Download data if required
if(DATA_DOWNLOAD)
  message(STATUS "Downloading data files")

  download_if_not_installed("http://mnovak.web.cern.ch/mnovak/data/cms2015.root"
    "cms2015.root"
    "${PROJECT_SOURCE_DIR}/test/cmstestdata/"
    "fb3ede867532b372c5e6f7138d00c07e")
  download_if_not_installed("http://mnovak.web.cern.ch/mnovak/data/ExN03.root"
    "ExN03.root"
    "${CMAKE_CURRENT_BINARY_DIR}"
    "b6b0cfdd5035117171bfe1b7f1f40c3f")
endif()

################################################################################
# Define executables
################################################################################
# - Dedicated Benchmark tests
if(VECGEOM_BENCHMARK)
  add_subdirectory(benchmark)
  build_executables("${TEST_EXECUTABLES_BENCHMARK}")
  add_to_ctest("${CTESTS_BENCHMARK}")
endif()

################################################################################
# - Static Analysis tests
if(STATIC_ANALYSIS)
  add_subdirectory(static_analysis)
endif()

################################################################################
# - Core Tests
if(CTEST)
  set(TEST_EXECUTABLES_CORE
    core/TestVoxelHashMap.cpp
    core/SafetyEstimatorTest.cpp
    core/ContainerTest.cpp
    core/create_geometry.cpp
    core/BitSetTest.cpp
    core/PlanesTest.cpp
    core/QuadrilateralTest.cpp
    # core/SOATest.cpp
    core/Transformation3DTest.cpp
    # core/boolminustest.cpp
    # core/boolminustest2.cpp
    core/PhiWedgeTest.cpp
    core/ThetaConeTest.cpp
    core/TestConvexity.cpp
    core/BooleanConvexityTest.cpp
    unit_tests/TestEstimateSurfaceArea.cpp
    unit_tests/TestVecGeomPolycone.cpp
    core/TestSExtru.cpp
    # to be enabled when running cleanly
    unit_tests/TestBooleans.cpp
    core/AssemblyTest.cpp
    core/TestMakeInstance.cpp
    #  services/CompNavStatePools.cpp
    core/TestMaskedAssign.cpp
    core/TestVector.cpp
    core/TestMap.cpp
    core/TestHybridBVH.cpp
    core/TestEarlyReturns.cpp
    core/CreateTessels.cpp
    core/CreateExtruded.cpp)

  if(USE_INDEXEDNAVSTATES)
    list(APPEND TEST_EXECUTABLES_CORE core/TestNavigationState.cpp)
  endif()

  if(NOT VECGEOM_ENABLE_CUDA)
    list(APPEND TEST_EXECUTABLES_CORE core/testVectorSafety.cpp)
  endif()

  build_executables("${TEST_EXECUTABLES_CORE}")
  add_to_ctest("${TEST_EXECUTABLES_CORE}")
endif()

################################################################################
# - ROOT Tests
if(VECGEOM_ROOT)
  set(TEST_EXECUTABLES_ROOT
    core/TestNavigationStatePool.cpp
    core/BoxBoxIntersectionTest.cpp
    core/SplittedABBox.cpp
    core/TestRegions.cpp
    root/root_geometry.cpp
    root/complex_test1.cpp
    root/E03Test.cpp
    root/ImportFromRootFileTest.cpp
    root/ImportTGeoPgon.cpp
    root/TestExportToROOT.cpp
    root/ExitingOrEntering.cpp
    root/AssemblyExample.cpp
    # higher level benchmarks or executables
    globalbenchmarks/LocatePointsBenchmark.cpp
    globalbenchmarks/SafetyKernelBenchmarker.cpp
    globalbenchmarks/NavigationKernelBenchmarker.cpp
    globalbenchmarks/TraceTrack.cpp
    globalbenchmarks/XRayBenchmarkFromROOTFile.cpp
    # ???
    ${PROJECT_SOURCE_DIR}/services/NavigationSpecializerTest.cpp
    ${PROJECT_SOURCE_DIR}/services/LibraryGenerator.cpp)

  if(VECGEOM_BENCHMARK)
    list(APPEND TEST_EXECUTABLES_ROOT
      root/BenchmarkShapeFromROOTFile.cpp
      root/BenchmarkShapeFromROOTFile_WithVisualization.cpp
      root/CompareDistances.cpp
      root/GenerateSurfacePoints.cpp
      shape_tester/shapeDebug.cpp)
  endif()

  if(VECGEOM_EMBREE)
    list(APPEND TEST_EXECUTABLES_ROOT core/EmbreeManagerTest.cpp)
  endif()

  # - separate list for ROOT UNIT tests
  set(TEST_UNITTESTEXECUTABLES_ROOT
    root/complex_test1.cpp
    root/E03Test.cpp
    root/TestExportToROOT.cpp
    root/ImportTGeoPgon.cpp)

  build_executables("${TEST_EXECUTABLES_ROOT}")
  add_to_ctest("${TEST_UNITTESTEXECUTABLES_ROOT}")
endif()

################################################################################
# - Visualization tests
if(VECGEOM_ROOT)
  set(TEST_EXECUTABLES_VISUALIZATION
    visualization/VisualizeTrap.cpp
    visualization/VisualizePolycone.cpp
    visualization/VisualizePolyhedron.cpp
    visualization/VisualizeParboloid.cpp
    visualization/VisualizeCone.cpp
    # visualization/VisualizeTorus.cpp
    visualization/VisualizeTube.cpp
    visualization/VisualizeScaled.cpp
    visualization/DebugTube.cpp
    visualization/DebugPolyhedron.cpp
    visualization/VisualizeSphere.cpp
    visualization/VisualizeGenTrap.cpp
    visualization/VisualizeHype.cpp
    visualization/VisualizeParallelepiped.cpp
    visualization/VisualizeCutTube.cpp
    # visualization/VisualizeExtruded.cpp
    )

  if(NOT VECGEOM_ENABLE_CUDA)
    list(APPEND TEST_EXECUTABLES_VISUALIZATION visualization/VisualizeMultiUnion.cpp)
  endif()

  build_executables("${TEST_EXECUTABLES_VISUALIZATION}")
endif()

################################################################################
# - Shape Tests
set(TEST_EXECUTABLES_SHAPES
  unit_tests/TestBox.cpp
  unit_tests/TestCons.cpp
  unit_tests/TestGenTrap.cpp
  unit_tests/TestTube.cpp
  unit_tests/TestEllipticalTube.cpp
  unit_tests/TestEllipticalCone.cpp
  unit_tests/TestEllipsoid.cpp
  unit_tests/TestCoaxialCones.cpp
  unit_tests/TestGenericPolycone.cpp
  unit_tests/TestHype.cpp
  unit_tests/TestTrd.cpp
  unit_tests/TestTrap.cpp
  unit_tests/TestParallelepiped.cpp
  unit_tests/TestPolycone.cpp
  unit_tests/TestPolyhedra.cpp
  unit_tests/TestTet.cpp
  unit_tests/TestOrb.cpp
  unit_tests/TestSphere.cpp
  unit_tests/TestTorus2.cpp
  unit_tests/TestParaboloid.cpp
  unit_tests/TestReducedPolycone.cpp
  unit_tests/TestUtils3D.cpp)

if(NOT VECGEOM_ENABLE_CUDA)
  list(APPEND TEST_EXECUTABLES_SHAPES
    unit_tests/TestTessellated.cpp
    unit_tests/TestMesh.cpp)
endif()

if(CTEST OR VALIDATION)
  build_executables("${TEST_EXECUTABLES_SHAPES}")
  # add unit tests
  add_test(NAME TestBox COMMAND TestBox)
  add_test(NAME TestTrap COMMAND TestTrap)
  add_test(NAME TestTube COMMAND TestTube)
  add_test(NAME TestCone COMMAND TestCons)
  add_test(NAME TestTrd COMMAND TestTrd)
  add_test(NAME TestOrb COMMAND TestOrb)
  add_test(NAME TestPolycone COMMAND TestPolycone)
  add_test(NAME TestPolyhedra COMMAND TestPolyhedra)
  add_test(NAME TestParallelepiped COMMAND TestParallelepiped)
  add_test(NAME TestHype COMMAND TestHype)
  add_test(NAME TestGenTrap COMMAND TestGenTrap)
  if(NOT VECGEOM_ENABLE_CUDA)
    add_test(NAME TestTessellated COMMAND TestTessellated)
  endif()
endif()

################################################################################
# - CUDA tests
if(VECGEOM_ENABLE_CUDA)
  set(TEST_EXECUTABLES_CUDA
    cuda/MapTest.cpp
    cuda/MapTestClass.cpp)
  build_executables("${TEST_EXECUTABLES_CUDA}")
endif()


################################################################################
# - Shape Tester tests
set(TEST_EXECUTABLES_SHAPETESTER
  shape_tester/shape_testBox.cpp
  shape_tester/shape_testSExtru.cpp
  shape_tester/shape_testOrb.cpp
  shape_tester/shape_testSphere.cpp
  shape_tester/shape_testCone.cpp
  shape_tester/shape_testEllipticalCone.cpp
  shape_tester/shape_testTube.cpp
  shape_tester/shape_testEllipticalTube.cpp
  shape_tester/shape_testEllipsoid.cpp
  shape_tester/shape_testHype.cpp
  shape_tester/shape_testTrd.cpp
  shape_tester/shape_testTrapezoid.cpp
  shape_tester/shape_testTet.cpp
  shape_tester/shape_testParaboloid.cpp
  shape_tester/shape_testPolycone.cpp
  shape_tester/shape_testGenericPolycone.cpp
  shape_tester/shape_testGenTrap.cpp
  shape_tester/shape_testParallelepiped.cpp
  shape_tester/convention_testTube.cpp
  shape_tester/shape_testPolyhedron.cpp
  shape_tester/shape_testTorus2.cpp
  shape_tester/shape_testCutTube.cpp
  shape_tester/shape_testExtruded.cpp
  shape_tester/shape_testMultiUnion.cpp)

if(VECGEOM_ROOT)
  list(APPEND TEST_EXECUTABLES_SHAPETESTER
    shape_tester/shape_testFromROOTFile.cpp
    shape_tester/shape_debugFromROOTFile.cpp)
endif()

if(NOT VECGEOM_ENABLE_CUDA)
  list(APPEND TEST_EXECUTABLES_SHAPETESTER shape_tester/shape_testTessellated.cpp)
endif()

if(VECGEOM_BENCHMARK)
  build_executables("${TEST_EXECUTABLES_SHAPETESTER}")
endif()


if(VALIDATION AND VECGEOM_ROOT)
  macro(add_cmsshapevalidation_test TESTNAME SHAPEFILE)
    add_test(NAME ${TESTNAME} COMMAND bash -c "${PROJECT_SOURCE_DIR}/test/scripts/RunRandomValidation ${PROJECT_BINARY_DIR}/BenchmarkShapeFromROOTFile  ${PROJECT_SOURCE_DIR}/test/cmstestdata/cms2015.root  ${PROJECT_SOURCE_DIR}/test/cmstestdata/${SHAPEFILE}" )
  endmacro()

  # Adding various shapes tests for nightlies
  add_cmsshapevalidation_test( tubevalidation cmstubes.txt )
  add_cmsshapevalidation_test( trapvalidation cmstraps.txt )
  add_cmsshapevalidation_test( polyconevalidation cmspolycones.txt )
  add_cmsshapevalidation_test( polyhedravalidation cmspolyhedra.txt )
  add_cmsshapevalidation_test( conevalidation cmscones.txt )
  add_cmsshapevalidation_test( boxvalidation cmsboxes.txt )
  #taken out due to problems: add_cmsshapevalidation_test( booleanvalidation cmsbooleans.txt )
  #taken out due to Issue-133: add_cmsshapevalidation_test( torusvalidation cmstori.txt )
endif()
