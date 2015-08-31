set(FLAGS_SSE   "-msse")
set(FLAGS_SSE2  "-msse2")
set(FLAGS_SSE3  "-msse3")
set(FLAGS_SSSE3 "-mssse3")
set(FLAGS_SSE41 "-msse4.1")
set(FLAGS_SSE42 "-msse4.2")
set(FLAGS_AVX   "-mavx")
set(FLAGS_AVX2  "-mf16c -mavx2 -mfma -mlzcnt -mbmi -mbmi2")
set(FLAGS_AVX512 "-mavx512")
set(FLAGS_NATIVE "-march=native")

set(CMAKE_CXX_FLAGS_DEBUG          "-DDEBUG  -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE        "-DNDEBUG    -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DDEBUG  -g -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -fno-strict-aliasing")
set(CMAKE_EXE_LINKER_FLAGS "")
