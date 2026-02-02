# This Makefile is just used for internal development - use CMake directly instead

all: out/analyse

clean:
	rm -rf out

out/build: CMakeLists.txt ../CMakeLists.txt
	mkdir -p out
	cmake -B out/build -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=.. -DCMAKE_BUILD_TYPE=Release

out/analyse: out/build source/main.cpp ../include/pitchtrack/pitchtrack.h
	cmake --build out/build --config Release

dev: out/analyse
	out/analyse data/vocadito/Audio/vocadito_13.wav out/vocadito_13.csv --plot=out/energy.svg --plot-ref=data/vocadito/Annotations/F0/vocadito_13_f0.csv

nodemon:
	nodemon -e h,hpp,hxx,c,cpp --watch source --watch ../include --exec "$(MAKE) dev || exit 1"
