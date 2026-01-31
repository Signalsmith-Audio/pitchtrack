#---- boilerplate ----#

cmake_minimum_required(VERSION 3.28)

set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)
set(CMAKE_VERBOSE_MAKEFILE on)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	set(CMAKE_CXX_STANDARD 20) # MSVC is strict about designated initialisers
	set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
if (APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET 10.15)
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64")
endif ()
if (NOT GENERATOR_IS_MULTI_CONFIG AND NOT CMAKE_BUILD_TYPE)
	message(WARNING "Using single-config generator, but CMAKE_BUILD_TYPE isn't set")
endif()

#----

project(signalsmith-pitch-tracker-cmd VERSION 0.0.0)

add_executable(analyse
	${CMAKE_CURRENT_SOURCE_DIR}/source/main.cpp	
)

add_subdirectory(../ parent-build)
target_link_libraries(analyse PUBLIC signalsmith-pitch-tracker)

include(FetchContent)
FetchContent_Declare(
	signalsmith-plot
	GIT_REPOSITORY https://github.com/Signalsmith-Audio/plot.git
	GIT_TAG 1.0.2
	GIT_SHALLOW ON
)
FetchContent_Declare(
	signalsmith-util
	GIT_REPOSITORY https://github.com/geraintluff/util.git
	GIT_TAG 2025-12-13
	GIT_SHALLOW ON
)
FetchContent_Declare(
	signalsmith-dsp
	GIT_REPOSITORY https://github.com/Signalsmith-Audio/dsp.git
	GIT_TAG v1.7.0
	GIT_SHALLOW ON
)
FetchContent_Declare(
	vocadito-dataset
	URL https://zenodo.org/records/5578807/files/vocadito.zip?download=1
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data/vocadito
	DOWNLOAD_NO_PROGRESS NO
)
FetchContent_MakeAvailable(signalsmith-plot signalsmith-util signalsmith-dsp vocadito-dataset)

target_link_libraries(analyse PUBLIC signalsmith-plot signalsmith-util signalsmith-dsp)