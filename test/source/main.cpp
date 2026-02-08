#include <iostream>
#define LOG_EXPR(expr) std::cout << #expr " = " << (expr) << std::endl;

#include "signalsmith-pitchtrack/pitchtrack.h"

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
	bool plotUnvoiced = args.hasFlag("plot-unvoiced", "whether to include unvoiced frequency results");
	double highpassHz = args.flag("highpass", "highpass (0=off)", 0);
	double wavRate = args.flag<double>("rate", "change playback rate of input audio", 1);
	args.errorExit();
	
	if (outCsvPath.empty()) outCsvPath = inWavPath + ".csv";
	std::ofstream csvOut(outCsvPath);
	if (!csvOut) args.errorExit("Couldn't open output CSV");
	
	Wav wav;
	if (!wav.read(inWavPath).warn()) args.errorExit("failed to open WAV");
	wav.sampleRate *= wavRate;
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
	tracker.analyse(input.data(), input.size());
	auto seconds = stopwatch.seconds(stopwatch.lap());

	auto inputSeconds = wav.length()/wav.sampleRate;
	double percent = seconds/inputSeconds*100;
	std::cout << "Processed in " << seconds << " seconds (" << percent << " % CPU)\n";

	if (hasPlot) { // Process again
		signalsmith::plot::Figure figure;
		figure.style.precision = 10; // the SVGs are huge, take the edge off a bit
		auto &energyPlot = figure(0, 0).plot(600, 200);
		energyPlot.x.major(0).minors(20, 100, 200, 500, 1000, 2000, 5000).range(hzToLogIsh, 20, 6000).label("Hz");
		energyPlot.y.major(0).minors(-20, -40, -60, -80).linear(-80, 0).label("dB");
		auto &scoreY = energyPlot.newY().flip().tick(0).label("score");
		auto &legend = energyPlot.legend(2, 1);

		auto &notePlot = figure(0, 1).plot(600, 150);
		notePlot.x.label("seconds");
		notePlot.y.major(60).minors(24, 36, 48, 72, 84, 96, 108).label("MIDI note");
		auto &noteLegend = notePlot.legend(2, 1);

		LOG_EXPR(svgPlotRef);
		static constexpr float note0Hz = 8.175798915643707;
		if (!svgPlotRef.empty()) {
			auto &refLine = energyPlot.line(4);
			auto &noteRefLine = notePlot.line(4);
			auto &noteMarkerLine = notePlot.line(4);
			legend.add(refLine, "reference");
			noteLegend.add(noteRefLine, "reference");
			
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
				row[0] /= wavRate;
				row[1] *= wavRate;
				
				if (row[1] == prevHz) continue;
				prevHz = row[1];
				if (row[1] > 0) {
					for (size_t h = 1; h <= 10; ++h) {
						refLine.cut();
						refLine.add(row[1]*h, 0);
						refLine.add(row[1]*h, -200);
					}
					noteRefLine.add(row[0], 12*std::log2(row[1]/note0Hz));
					noteMarkerLine.marker(row[0], 12*std::log2(row[1]/note0Hz));
				} else {
					noteRefLine.cut();
				}
				refLine.toFrame(row[0]/svgPlotSpeed);
				noteMarkerLine.toFrame(row[0]/svgPlotSpeed);
			}
		}

		size_t overlap = 4;
		{
			size_t blockSize = 16;
			tracker.configure(wav.sampleRate, blockSize, 0, overlap);
			tracker.reset();
			auto &hz = tracker.debugHz();
			auto &resultLine = energyPlot.line(energyPlot.x, energyPlot.newY().linear(0, 1).blank(), 5);
			auto &noteLine = notePlot.line();
			auto &noteMarkerLine = notePlot.line(noteLine.styleIndex);
			auto &tonalityY = notePlot.newY().flip().linear(0, 1).minors(0, 1).label("tonality");
			auto &tonalityLine = notePlot.line(notePlot.x, tonalityY);
			auto &tonalityMarkerLine = notePlot.line(notePlot.x, tonalityY, tonalityLine.styleIndex);
			auto &noteTimeLine = notePlot.line(notePlot.x, tonalityY);
			auto &energyLine = energyPlot.line();
			auto &scoreLine = energyPlot.line(energyPlot.x, scoreY);
			auto &scoreRefLine = energyPlot.line(energyPlot.x, scoreY);
			auto &scoreKernelLine = energyPlot.fill(energyPlot.x, energyPlot.newY().blank()).fillToY(0);
			legend.add(energyLine, "energy").add(scoreLine, "score (result)").add(scoreRefLine, "score (raw)").fill(scoreKernelLine, "score kernel");
			noteLegend.add(noteLine, "note estimate").add(tonalityLine, "tonality").add(noteTimeLine, "time");
			
			for (size_t blockStart = 0; blockStart < input.size(); blockStart += blockSize) {
				size_t actualBlockSize = std::min<size_t>(blockSize, input.size() - blockStart);
				tracker.analyze(input.data() + blockStart, actualBlockSize);
				if (tracker.updates().empty()) continue;

				auto &update = tracker.updates()[0];
				int sampleIndex = int(blockStart + update.offset);// - tracker.latencySamples();
				if (sampleIndex < 0) continue; // it breaks for t < 0

				auto &energy = tracker.debugEnergy();
				auto &scoreRef = tracker.debugScoreRef();
				auto &score = tracker.debugScore();
				size_t bestB = 0;
				float bestBScore = 0;
				for (size_t b = 0; b < hz.size(); ++b) {
					auto db = 10*std::log10(energy[b] + float(1e-30));
					energyLine.add(hz[b], db);

					scoreRefLine.add(hz[b], scoreRef[b]);
					scoreLine.add(hz[b], score[b]);

					if (score[b] > bestBScore) {
						bestBScore = score[b];
						bestB = b;
					}
				}
				scoreLine.marker(hz[bestB], score[bestB]);
				if (update.voiced) {
					resultLine.marker(update.hz, update.tonality);
				}
				for (size_t h = 1; h <= 10; ++h) {
					resultLine.cut();
					resultLine.add(update.hz*h, 0);
					resultLine.add(update.hz*h, update.tonality);
				}

				float midiNote = 12*std::log2(update.hz/note0Hz);
				float seconds = sampleIndex/float(wav.sampleRate);
				if (update.hz > 0 && (plotUnvoiced || update.voiced)) {
					noteLine.add(seconds, midiNote);
					tonalityLine.add(seconds, update.tonality);
					if (update.voiced) {
						noteMarkerLine.marker(seconds, midiNote);
						tonalityMarkerLine.marker(seconds, update.tonality);
					}
				} else {
					noteLine.cut();
					tonalityLine.cut();
				}
				noteTimeLine.add(seconds, 0);
				noteTimeLine.add(seconds, 1);

				auto frameTime = seconds/svgPlotSpeed;
				energyLine.toFrame(frameTime);
				scoreRefLine.toFrame(frameTime);
				scoreLine.toFrame(frameTime);
				resultLine.toFrame(frameTime);
				noteTimeLine.toFrame(frameTime);
				noteMarkerLine.toFrame(frameTime);
				tonalityMarkerLine.toFrame(frameTime);

				auto &scoreKernel = tracker.debugScoreKernel();
				for (size_t b = 0; b < scoreRef.size(); ++b) {
					size_t offsetB = (b + bestB)%scoreKernel.size();
					if (offsetB >= hz.size()) {
						scoreKernelLine.cut();
						continue;
					}
					float deltaHz = update.hz - hz[bestB];
					scoreKernelLine.add(hz[offsetB] + deltaHz, scoreKernel[b]);
				}
				scoreKernelLine.toFrame(frameTime);
			}
		}

		// synchronise loop times
		figure.loopFrame(wav.length()/svgPlotSpeed/wav.sampleRate);
		
		notePlot.x.major(0);
		for (size_t sec = 1; sec < wav.length()/wav.sampleRate; ++sec) {
			if (sec%10) {
				notePlot.x.tick(sec, "");
			} else {
				notePlot.x.tick(sec);
			}
		}

		std::cout << "writing SVG..." << std::endl;
		figure.write(svgPlotPath);
		std::cout << "\tdone\n";
	}

	for (auto &event : tracker.updates()) {
		csvOut << (event.offset/wav.sampleRate) << "," << event.hz << "," << event.tonality << "\n";
	}
}
