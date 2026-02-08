cmake_minimum_required(VERSION 3.24)

add_library(signalsmith-pitchtrack INTERFACE)
target_include_directories(signalsmith-pitchtrack INTERFACE
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

target_link_libraries(signalsmith-pitchtrack INTERFACE signalsmith-linear)