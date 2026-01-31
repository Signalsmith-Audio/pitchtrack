cmake_minimum_required(VERSION 3.24)

add_library(signalsmith-pitch-tracker INTERFACE)
target_include_directories(signalsmith-pitch-tracker INTERFACE
	${CMAKE_CURRENT_SOURCE_DIR}/include
)

include(FetchContent)
FetchContent_Declare(
	signalsmith-linear
	GIT_REPOSITORY https://github.com/Signalsmith-Audio/linear.git
	GIT_TAG 0.3.0
	GIT_SHALLOW ON
)
FetchContent_MakeAvailable(signalsmith-linear)

target_link_libraries(signalsmith-pitch-tracker INTERFACE signalsmith-linear)