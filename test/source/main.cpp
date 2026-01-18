#include <iostream>
#define LOG_EXPR(expr) std::cout << #expr " = " << (expr) << std::endl;

#include "pitchtrack/pitchtrack.h"

#include "signalsmith-util/simple-args.h"
#include "signalsmith-util/wav.h"
#include "signalsmith-util/stopwatch.h"

#include "signalsmith-plot/plot.h"

int main(int argc, char *argv[]) {
	SimpleArgs args{argc, argv};
	
	std::string inWavPath = args.arg("input.wav", "input WAV file");
	std::string outCsvPath = args.arg("output.csv", "output CSV", "");
	std::string svgPlotPath = args.flag("plot", "output SVG plot (optional)", "");
	args.errorExit();
	
	if (outCsvPath.empty()) outCsvPath = inWavPath + ".csv";
	std::ofstream csvOut(outCsvPath);
	if (!csvOut) args.errorExit("Couldn't open output CSV");
	
	Wav wav;
	if (!wav.read(inWavPath).warn()) args.errorExit("failed to open WAV");
	wav.makeMono();
	std::vector<float> input;
	input.resize(wav.samples.size());
	for (size_t i = 0; i < input.size(); ++i) input[i] = float(wav.samples[i]);
	
	signalsmith::pitchtrack::PitchTrack tracker;
	tracker.configure(wav.sampleRate);
	
	signalsmith::Stopwatch stopwatch;
	stopwatch.start();
	tracker.block(input.data(), input.size());
	auto seconds = stopwatch.seconds(stopwatch.lap());
	
	auto inputSeconds = wav.length()/wav.sampleRate;
	double percent = seconds/inputSeconds*100;
	std::cout << "Done in " << seconds << " seconds (" << percent << " % CPU)\n";
	
	for (auto &event : tracker.updates()) {
		csvOut << (event.offset/wav.sampleRate) << "," << event.hz << "," << event.tonality << "\n";
	}
}
