all: out/analyse

clean:
	rm -rf out

out/build: CMakeLists.txt ../CMakeLists.txt
	mkdir -p out
	cmake -B out/build -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=.. -DCMAKE_BUILD_TYPE=Release

out/analyse: out/build source/main.cpp ../include/pitchtrack/pitchtrack.h
	cmake --build out/build --config Release

nodemon:
	nodemon -e h,hpp,hxx,c,cpp --watch source --watch ../include --exec "$(MAKE) && out/analyse data/vocadito/Audio/vocadito_13.wav out/vocadito_13.csv || exit 1"