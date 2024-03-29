##################################################################
# IQ-TREE cmake build definition
# Copyright (c) 2012-2014 Bui Quang Minh, Lam Tung Nguyen
##################################################################

# Windows example usages:
# cmake -G "Visual Studio 12" <source_dir>			(32-bit SSE3 version, compiled with MSVC)
# cmake -G "Visual Studio 12 Win64" <source_dir>	(64-bit SSE3 version, compiled with MSVC)
# cmake -G "Visual Studio 12 Win64" -T "Intel C++ Compiler XE 15.0" <source_dir>		(64-bit SSE3 version, compiled with ICC)
# cmake -G "Visual Studio 12 Win64" -T "Intel C++ Compiler XE 15.0" -DIQTREE_FLAGS="avx" <source_dir>		(64-bit AVX version, compiled with ICC)
# cmake -G "Visual Studio 12 Win64" -T "Intel C++ Compiler XE 15.0" -DIQTREE_FLAGS="omp" <source_dir>		(64-bit SSE3+OpenMP version, compiled with ICC)
# cmake -G "Visual Studio 12 Win64" -T "Intel C++ Compiler XE 15.0" -DIQTREE_FLAGS="avx omp" <source_dir>	(64-bit AVX+OpenMP version, compiled with ICC)
#
# Linux example usages:
# cmake <source_dir>						   (SSE3 version)
# cmake -DIQTREE_FLAGS="avx" <source_dir>      (AVX version)
# cmake -DIQTREE_FLAGS="omp" <source_dir>      (OpenMP version)
# cmake -DIQTREE_FLAGS="omp avx" <source_dir>  (AVX+OpenMP version)
# cmake -DIQTREE_FLAGS="m32" <source_dir>      (32-bit SSE3 version)
# cmake -DIQTREE_FLAGS="m32 omp" <source_dir>  (32-bit SSE3+OpenMP version)
#
# Mac OSX example usages:
# cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ <source_dir>							(SSE3 version)
# cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DIQTREE_FLAGS="avx" <source_dir>		(AVX version)
# To build OpenMP version one needs to download and compile Clang-OpenMP. 
# Then assuming clang2 and clang2++ are the newly built compilers, then:
# cmake -DCMAKE_C_COMPILER=clang2 -DCMAKE_CXX_COMPILER=clang2++ -DIQTREE_FLAGS="omp" <source_dir>		(SSE3+OpenMP version)
# cmake -DCMAKE_C_COMPILER=clang2 -DCMAKE_CXX_COMPILER=clang2++ -DIQTREE_FLAGS="omp avx" <source_dir>	(AVX+OpenMP version)
#


cmake_minimum_required(VERSION 2.8)

project(mpboot)
add_definitions(-DIQ_TREE)
# The version number.
set (iqtree_VERSION_MAJOR 1)
set (iqtree_VERSION_MINOR 1)
set (iqtree_VERSION_PATCH 1)

set(BUILD_SHARED_LIBS OFF)

message("IQ-TREE flags : ${IQTREE_FLAGS}")

if (NOT CMAKE_BUILD_TYPE) 
	set(CMAKE_BUILD_TYPE "Release")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release") 
	message("Builde mode   : Release")
endif()

if (CMAKE_GENERATOR MATCHES "Xcode")
    set(CMAKE_XCODE_ATTRIBUTE_DEBUG_INFORMATION_FORMAT "dwarf-with-dsym")
endif()

include_directories("${PROJECT_SOURCE_DIR}")


##################################################################
# Detect target platforms
##################################################################
if (WIN32)
	message("Target OS     : Windows")
	# build as static binary to run on most machines
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    add_definitions(-DWIN32)
elseif (APPLE) 
	message("Target OS     : Mac OS X")
	# to be compatible back to Mac OS X 10.7
	add_definitions("-mmacosx-version-min=10.7")
	if (CMAKE_SYSTEM_PROCESSOR MATCHES "(ARM)|(arm)")
		set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE INTERNAL "" FORCE)
	endif()
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mmacosx-version-min=10.7")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
elseif (UNIX) 
	message("Target OS     : Unix")
	# build as static binary to run on most machines
  if (CMAKE_BUILD_TYPE STREQUAL "Release") 	
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
  endif() 
else()
	# Note that IQ-TREE has NOT been tested on other platforms
	message("Target OS     : Unknown and untested yet")
endif()

##################################################################
# Setup compiler, currently supported GCC, CLANG, MSVC, and ICC
##################################################################

if (MSVC) 
	# MS Visual Studio environment
	message("Exporting MS Visual Studio projects...")
	add_definitions(/MP) # enable multi-processor compilation
	if (CMAKE_BUILD_TYPE STREQUAL "Release")
		add_definitions(/O2 /Ot /Oi)
	endif()
endif()

set(GCC "FALSE")   #  GNU compiler
set(CLANG "FALSE") # Clang compiler
set(ICC "FALSE")   # Intel compiler
set(VCC "FALSE")   # MS Visual C Compiler, note that it is different from MSVC variable

if (CMAKE_COMPILER_IS_GNUCXX) 	
	message("Compiler      : GNU Compiler (gcc)")
	set(GCC "TRUE")
#	set(COMBINED_FLAGS "-Wall -Wno-unused-function -Wno-sign-compare -pedantic -D_GNU_SOURCE -fms-extensions -Wno-deprecated")
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g0")
	set(CMAKE_C_FLAGS_RELEASE "-O3 -g0")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	message("Compiler      : Clang")
	set(CLANG "TRUE")
#	set(COMBINED_FLAGS "-Wall -Wno-unused-function -Wno-sign-compare -pedantic -D_GNU_SOURCE -Wno-nested-anon-types")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
	set(CMAKE_C_FLAGS_RELEASE "-O3")
	if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,-dead_strip")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
		SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    endif()	
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
	set(VCC "TRUE")
	message("Compiler      : MS Visual C++ Compiler")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
	message("Compiler      : Intel C++ Compiler (icc)")
	set(ICC "TRUE")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qstd=c99")
else()
	message("Compiler      : Unknown and untested yet")
endif()

set(EXE_SUFFIX "")

##################################################################
# Configure PLL build
##################################################################
if (IQTREE_FLAGS MATCHES "pll")
	add_definitions(-DUSING_PLL)
	set(EXE_SUFFIX "${EXE_SUFFIX}-pll")
endif()

##################################################################
# detect 32 or 64 bit binary
##################################################################
set (BINARY32 "FALSE")
if(CMAKE_SIZEOF_VOID_P EQUAL 4 OR IQTREE_FLAGS MATCHES "m32")
	set(BINARY32 "TRUE")
	message("Target binary : 32-bit")
	if (CMAKE_GENERATOR MATCHES "Win64")
		error("Both 32-bit and 64-bit mode cannot be specified")
	endif()
	SET(EXE_SUFFIX "${EXE_SUFFIX}-m32")
	if (GCC OR CLANG) 
		set(COMBINED_FLAGS "${COMBINED_FLAGS} -m32")
  	endif()
else()
	message("Target binary : 64-bit")
endif()

##################################################################
# configure OpenMP/PThreads compilation
# change the executable name if compiled for OpenMP parallel version
##################################################################
if (IQTREE_FLAGS MATCHES "omp")
	message("Parallel      : OpenMP/PThreads")
	SET(EXE_SUFFIX "${EXE_SUFFIX}-omp")	
	add_definitions(-D_USE_PTHREADS)
	if (MSVC) 
		add_definitions(/MT)
	endif()
	
	if (VCC) 
  		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
  		include_directories("${PROJECT_SOURCE_DIR}/pll") # for PThreads headers 
	elseif (ICC)
  		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qopenmp")
  		if (WIN32)
  			include_directories("${PROJECT_SOURCE_DIR}/pll") # for PThreads headers
  		endif() 
  	elseif (GCC)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
  		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  	elseif (CLANG) 
  		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  	endif()
else()
	message("Parallel      : None")
endif()

##################################################################
# configure SSE/AVX/FMA instructions
##################################################################
if (IQTREE_FLAGS MATCHES "fma") # AVX+FMA instruction set
 	message("Vectorization : AVX+FMA")
	add_definitions(-D__SSE3 -D__AVX) # define both SSE3 and AVX directive
	if (VCC)
		# Visual C++ has no /mfma flag!, FMA is only included in AVX2 
		set(COMBINED_FLAGS "${COMBINED_FLAGS} /arch:AVX2")
	elseif (CLANG)
		set(COMBINED_FLAGS "${COMBINED_FLAGS} -mavx -mfma")
	elseif (GCC)
		set(COMBINED_FLAGS "${COMBINED_FLAGS} -mavx -fabi-version=0 -mfma")
	elseif (ICC) 
		if (WIN32)
			 set(COMBINED_FLAGS "${COMBINED_FLAGS} /arch:AVX /Qfma")
		else()
			 set(COMBINED_FLAGS "${COMBINED_FLAGS} -mavx -mfma")
		endif()
	endif()
	SET(EXE_SUFFIX "${EXE_SUFFIX}-fma")

elseif (IQTREE_FLAGS MATCHES "avx2") # AVX instruction set
 	message("Vectorization : AVX2")
	add_definitions(-D__SSE3 -D__AVX) # define both SSE3 and AVX directive
	if (VCC) 
		set(COMBINED_FLAGS "${COMBINED_FLAGS} /arch:AVX2")
	elseif (CLANG)
		set(COMBINED_FLAGS "${COMBINED_FLAGS} -mavx2 -mfma")
	elseif (GCC)
		set(COMBINED_FLAGS "${COMBINED_FLAGS} -mavx2 -fabi-version=0 -mfma")
	elseif (ICC) 
		if (WIN32)
			 set(COMBINED_FLAGS "${COMBINED_FLAGS} /arch:AVX2 /Qfma")
		else()
			 set(COMBINED_FLAGS "${COMBINED_FLAGS} -mavx2 -mfma")
		endif()
	endif()
	SET(EXE_SUFFIX "${EXE_SUFFIX}-avx2")
elseif (IQTREE_FLAGS MATCHES "avx") # AVX instruction set
 	message("Vectorization : AVX")
	add_definitions(-D__SSE3 -D__AVX) # define both SSE3 and AVX directive
	if (VCC) 
		set(COMBINED_FLAGS "${COMBINED_FLAGS} /arch:AVX")
	elseif (CLANG)
		set(COMBINED_FLAGS "${COMBINED_FLAGS} -mavx")
	elseif (GCC)
		set(COMBINED_FLAGS "${COMBINED_FLAGS} -mavx -fabi-version=0")
	elseif (ICC) 
		if (WIN32)
			 set(COMBINED_FLAGS "${COMBINED_FLAGS} /arch:AVX")
		else()
			 set(COMBINED_FLAGS "${COMBINED_FLAGS} -mavx")
		endif()
	endif()
	SET(EXE_SUFFIX "${EXE_SUFFIX}-avx")
elseif (IQTREE_FLAGS MATCHES "sse4") #SSE intruction set
	message("Vectorization : SSE4.1")
	add_definitions(-D__SSE3)
	
	if (VCC)
		set(COMBINED_FLAGS "${COMBINED_FLAGS} /arch:SSE4.1 -D__SSE4.1__")
	elseif (GCC OR CLANG)
		set(COMBINED_FLAGS "${COMBINED_FLAGS} -msse4.1")
	elseif (ICC)
		if (WIN32)
			set(COMBINED_FLAGS "${COMBINED_FLAGS} /arch:SSE4.1")
		else()
			set(COMBINED_FLAGS "${COMBINED_FLAGS} -msse4.1")
		endif()
	endif()
elseif (IQTREE_FLAGS MATCHES "sse") #SSE intruction set
	message("Vectorization : SSE3")
	add_definitions(-D__SSE3)
	
	if (VCC)
		set(COMBINED_FLAGS "${COMBINED_FLAGS} /arch:SSE2 -D__SSE3__")
	elseif (GCC OR CLANG)
		set(COMBINED_FLAGS "${COMBINED_FLAGS} -msse3")
	elseif (ICC)
		if (WIN32)
			set(COMBINED_FLAGS "${COMBINED_FLAGS} /arch:SSE3")
		else()
			set(COMBINED_FLAGS "${COMBINED_FLAGS} -msse3")
		endif()
	endif()
endif()

##################################################################
# Setup compiler flags
##################################################################
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMBINED_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMBINED_FLAGS}")

if (CMAKE_BUILD_TYPE STREQUAL "Release")
	message("C flags    : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}") 
	message("CXX flags  : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}") 
endif()

set(CMAKE_CXX_FLAGS_PROFILE "-fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -fno-default-inline -fno-inline -O0 -fno-omit-frame-pointer -pg")
set(CMAKE_C_FLAGS_PROFILE "-fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -O0 -fno-omit-frame-pointer -pg")

if (GCC) 
	set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline-functions-called-once -fno-default-inline -fno-inline")
	set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline-functions-called-once -fno-default-inline -fno-inline")
	set(CMAKE_CXX_FLAGS_MEM "-g -O1")
	set(CMAKE_C_FLAGS_MEM "-g -O1")
elseif (CLANG)
	set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline")
	set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline")
	set(CMAKE_CXX_FLAGS_MEM "-g -O1")
	set(CMAKE_C_FLAGS_MEM "-g -O1")
endif()

##################################################################
# check existence of a few basic functions
##################################################################
include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
check_function_exists (gettimeofday HAVE_GETTIMEOFDAY)
check_function_exists (getrusage HAVE_GETRUSAGE)
check_function_exists (GlobalMemoryStatusEx HAVE_GLOBALMEMORYSTATUSEX)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/iqtree_config.h.in"
  "${PROJECT_BINARY_DIR}/iqtree_config.h"
  )

# add the binary tree to the search path for include files
# so that we will find iqtree_config.h
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/zlib-1.2.7")


##################################################################
# subdirectories containing necessary libraries for the build 
##################################################################
add_subdirectory(pllrepo/src)
add_subdirectory(ncl)
add_subdirectory(whtest)
add_subdirectory(sprng)
add_subdirectory(zlib-1.2.7)
add_subdirectory(vectorclass)
add_subdirectory(model)

##################################################################
# the main executable
##################################################################
add_executable(mpboot
alignment.cpp
alignmentpairwise.cpp
circularnetwork.cpp
eigendecomposition.cpp
greedy.cpp
gss.cpp
guidedbootstrap.cpp
gurobiwrapper.cpp
gzstream.cpp
hashsplitset.cpp
iqtree.cpp
maalignment.cpp
matree.cpp
mexttree.cpp
mpdablock.cpp
msetsblock.cpp
msplitsblock.cpp
mtree.cpp
mtreeset.cpp
ncbitree.cpp
ngs.cpp
node.cpp
optimization.cpp
parsmultistate.cpp
pattern.cpp
pda.cpp
pdnetwork.cpp
pdtree.cpp
pdtreeset.cpp
phyloanalysis.cpp
phylonode.cpp
phylosupertree.cpp
phylotree.cpp
phylotreesse.cpp
#phylotreeeigen.cpp
pruning.cpp
split.cpp
splitgraph.cpp
splitset.cpp
stoprule.cpp
superalignment.cpp
superalignmentpairwise.cpp
supernode.cpp
tinatree.cpp
tools.cpp
whtest_wrapper.cpp
lpwrapper.c
nnisearch.cpp
phylosupertreeplen.cpp
phylotesting.cpp
ecopd.cpp
ecopdmtreeset.cpp
graph.cpp
candidateset.cpp
checkpoint.cpp
parstree.cpp
sprparsimony.cpp
test.cpp
)

##################################################################
# setup linking flags
##################################################################

# link special lib for WIN32
if (WIN32) 
	set(PLATFORM_LIB "ws2_32")
else()
	set(PLATFORM_LIB "m")
endif()

set(THREAD_LIB "")
if (IQTREE_FLAGS MATCHES "omp") 
	link_directories(${PROJECT_SOURCE_DIR}/lib)
	if (MSVC)
		set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LIBPATH:${PROJECT_SOURCE_DIR}/lib")
		if (BINARY32) 
			set(THREAD_LIB pthreadVC2_x86)
		else()
			set(THREAD_LIB pthreadVC2)
		endif()
	elseif(CLANG AND APPLE)
		set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${PROJECT_SOURCE_DIR}/lib")
	endif()
endif()

target_link_libraries(mpboot pll ncl whtest zlibstatic sprng vectorclass model ${PLATFORM_LIB} ${THREAD_LIB})	

##################################################################
# setup the executable name 
##################################################################
set_target_properties(mpboot PROPERTIES OUTPUT_NAME "mpboot${EXE_SUFFIX}")

# strip the release build
if (CMAKE_BUILD_TYPE STREQUAL "Release" AND (GCC OR CLANG)) # strip is not necessary for MSVC
	if (WIN32)
		ADD_CUSTOM_COMMAND(TARGET mpboot POST_BUILD COMMAND strip $<TARGET_FILE:mpboot>)
	else()
		ADD_CUSTOM_COMMAND(TARGET mpboot POST_BUILD COMMAND ${CMAKE_STRIP} $<TARGET_FILE:mpboot>)
	endif()
endif()

if (MSVC)
	set (BINARY_DIR "${PROJECT_BINARY_DIR}/Release")
else()
	set (BINARY_DIR "${PROJECT_BINARY_DIR}")
endif()

if (WIN32)
	if (MSVC) 
		ADD_CUSTOM_COMMAND(TARGET mpboot POST_BUILD COMMAND copy "Release\\mpboot${EXE_SUFFIX}.exe" "Release\\mpboot${EXE_SUFFIX}-click.exe")
	else()
		ADD_CUSTOM_COMMAND(TARGET mpboot POST_BUILD COMMAND copy "mpboot${EXE_SUFFIX}.exe" "mpboot${EXE_SUFFIX}-click.exe")
	endif()	
endif()

##############################################################
# add the install targets
##############################################################
install (TARGETS mpboot DESTINATION bin)
install (FILES "${PROJECT_SOURCE_DIR}/examples/example.phy" DESTINATION .)
install (FILES "${PROJECT_SOURCE_DIR}/Documents/iqtree-manual-1.0.pdf" DESTINATION .)

if (WIN32)
	install (FILES "${BINARY_DIR}/mpboot${EXE_SUFFIX}-click.exe" DESTINATION bin)
	if (EXE_SUFFIX MATCHES "omp" AND MSVC)
		install(FILES  "${PROJECT_SOURCE_DIR}/lib/pthreadVC2.dll" DESTINATION bin)
		install(FILES  "${PROJECT_SOURCE_DIR}/lib/libiomp5md.dll" DESTINATION bin)
#		install(FILES  "${PROJECT_SOURCE_DIR}/lib/pthreadVC2_x86.dll" DESTINATION bin)
#		install(FILES  "${PROJECT_SOURCE_DIR}/lib/pthreadGC2.dll" DESTINATION bin)
#		install(FILES  "${PROJECT_SOURCE_DIR}/lib/pthreadGC2_64.dll" DESTINATION bin)
	endif()
endif()

##############################################################
# build a CPack driven installer package
##############################################################
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE  
     "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${iqtree_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${iqtree_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${iqtree_VERSION_PATCH}")
if(WIN32 OR APPLE)
  set(CPACK_GENERATOR "ZIP")
  set(CPACK_SOURCE_GENERATOR "ZIP")
else()
  set(CPACK_GENERATOR "TGZ")
  set(CPACK_SOURCE_GENERATOR "TGZ")
endif()

#set(CPACK_SOURCE_PACKAGE_FILE_NAME
#  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
set(CPACK_SOURCE_IGNORE_FILES
  "/build.*/;/debug.*/;/examples/;/manual/;/.bzr/;~$;/\\\\.svn/;/\\\\.git/;${CPACK_SOURCE_IGNORE_FILES}")

set (SYSTEM_NAME "${CMAKE_SYSTEM_NAME}")
if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin") 
	set (SYSTEM_NAME "MacOSX")
endif()

set (PROJECT_NAME_SUFFIX "${EXE_SUFFIX}")
if (NOT IQTREE_FLAGS MATCHES "omp" AND NOT IQTREE_FLAGS MATCHES "avx" AND NOT IQTREE_FLAGS MATCHES "fma")  
	set (PROJECT_NAME_SUFFIX "${PROJECT_NAME_SUFFIX}-sse") 
endif()

set(CPACK_PACKAGE_FILE_NAME 
	"${CMAKE_PROJECT_NAME}${PROJECT_NAME_SUFFIX}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-${SYSTEM_NAME}")

set(CPACK_STRIP_FILES TRUE)

include (CPack)

add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
