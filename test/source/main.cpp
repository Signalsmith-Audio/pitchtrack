#include <iostream>
#define LOG_EXPR(expr) std::cout << #expr " = " << (expr) << std::endl;

#include "pitchtrack/pitchtrack.h"

#include "signalsmith-util/simple-args.h"
#include "signalsmith-util/wav.h"
#include "signalsmith-util/stopwatch.h"

#include "signalsmith-plot/plot.h"
#include "signalsmith-dsp/curves.h"
#include "signalsmith-dsp/filters.h"

#include <array>

double hzToBark(double hz) {
	auto barkScale = signalsmith::curves::Reciprocal<double>::barkScale();
	return barkScale.inverse(hz);
}
double hzToLogIsh(double hz) {
	return std::log(hz + 20);
}

int main(int argc, char *argv[]) {
	SimpleArgs args{argc, argv};
	
	std::string inWavPath = args.arg("input.wav", "input WAV file");
	std::string outCsvPath = args.arg("output.csv", "output CSV", "");
	std::string svgPlotPath = args.flag("plot", "output SVG plot (optional)", "");
	double svgPlotSpeed = args.flag("plot-speed", "plot animation rate (compared to realtime)", 0.1);
	bool hasPlot = !svgPlotPath.empty();
	std::string svgPlotRef = args.flag("plot-ref", "reference (seconds,hz) CSV", "");
	double highpassHz = args.flag("highpass", "highpass (0=off)", 0);
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
	if (highpassHz > 0) {
		std::array<signalsmith::filters::BiquadStatic<float>, 2> filters;
		for (auto &f : filters) f.highpass(highpassHz/wav.sampleRate);
		for (auto &v : input) {
			for (auto &f : filters) v = f(v);
		}
	}
	
	signalsmith::pitchtrack::PitchTrack tracker;
	tracker.configure(wav.sampleRate, input.size());
	
	signalsmith::Stopwatch stopwatch;
	stopwatch.start();
	// Process everything at once
	tracker.block(input.data(), input.size());
	auto seconds = stopwatch.seconds(stopwatch.lap());

	auto inputSeconds = wav.length()/wav.sampleRate;
	double percent = seconds/inputSeconds*100;
	std::cout << "Processed in " << seconds << " seconds (" << percent << " % CPU)\n";

	if (hasPlot) { // Process again
		signalsmith::plot::Figure figure;
		auto &energyPlot = figure.plot(600, 200);
		energyPlot.x.major(0).minors(20, 100, 200, 500, 1000, 2000, 5000).range(hzToLogIsh, 20, 6000).label("Hz");
		energyPlot.y.major(0).minors(-20, -40, -60, -80).linear(-80, 0).label("dB");
		auto &scoreY = energyPlot.newY().flip().minors(-1, 1).label("score");
		auto &legend = energyPlot.legend(1, 1);

		LOG_EXPR(svgPlotRef);
		if (!svgPlotRef.empty()) {
			auto &refLine = energyPlot.line(4);
			legend.add(refLine, "reference");
			
			std::ifstream csv{svgPlotRef};
			if (!csv) args.errorExit("Failed to open CSV");
			std::string line, colStr;
			std::vector<double> row;
			float prevHz = -1;
			while (std::getline(csv, line)) {
				row.resize(0);
				std::istringstream sstr{line};
				while (std::getline(sstr, colStr, ',')) {
					row.push_back(std::stod(colStr));
				}
				
				if (row[1] == prevHz) continue;
				prevHz = row[1];
				if (row[1] > 0) {
					for (size_t h = 1; h <= 10; ++h) {
						refLine.cut();
						refLine.add(row[1]*h, 0);
						refLine.add(row[1]*h, -200);
					}
				}
				refLine.toFrame(row[0]/svgPlotSpeed);
			}
		}

		size_t overlap = 4;
		{
			size_t blockSize = 16;
			tracker.configure(wav.sampleRate, blockSize, 0, overlap);
			tracker.reset();
			auto &hz = tracker.debugHz();
			auto &energyLine = energyPlot.line();
			auto &scoreLine = energyPlot.line(energyPlot.x, scoreY);
			auto &scoreRefLine = energyPlot.line(energyPlot.x, scoreY);
			legend.add(energyLine, "energy").add(scoreLine, "score").add(scoreRefLine, "score kernel");
			
			for (size_t blockStart = 0; blockStart < input.size(); blockStart += blockSize) {
				size_t actualBlockSize = std::min<size_t>(blockSize, input.size() - blockStart);
				tracker.block(input.data() + blockStart, actualBlockSize);
				if (tracker.updates().empty()) continue;

				auto &update = tracker.updates()[0];
				int sampleIndex = int(blockStart + update.offset) - tracker.latencySamples();
				if (sampleIndex < 0) continue; // it breaks for t < 0

				auto &energy = tracker.debugEnergy();
				auto &score = tracker.debugScore();
				size_t bestB = 0;
				float bestBScore = 0;
				for (size_t b = 0; b < hz.size(); ++b) {
					auto db = 10*std::log10(energy[b] + float(1e-30));
					energyLine.add(hz[b], db);

					scoreLine.add(hz[b], score[b]);

					if (score[b] > bestBScore) {
						bestBScore = score[b];
						bestB = b;
					}
				}
				scoreLine.marker(hz[bestB], score[bestB]);

				auto frameTime = sampleIndex/svgPlotSpeed/wav.sampleRate;
				energyLine.toFrame(frameTime);
				scoreLine.toFrame(frameTime);
			}

			auto &scoreRef = tracker.debugScoreKernel();
			for (size_t b = 0; b < hz.size(); ++b) {
				scoreRefLine.add(hz[b], scoreRef[b]);
			}
		}

		// synchronise loop times
		figure.loopFrame(wav.length()/svgPlotSpeed/wav.sampleRate);

		figure.write(svgPlotPath);
	}

	for (auto &event : tracker.updates()) {
		csvOut << (event.offset/wav.sampleRate) << "," << event.hz << "," << event.tonality << "\n";
	}
}
