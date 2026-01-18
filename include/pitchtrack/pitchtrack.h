#pragma once

#include "signalsmith-linear/stft.h"
#include "signalsmith-linear/linear.h"
#include "signalsmith-dsp/curves.h"

#include <vector>

namespace signalsmith { namespace pitchtrack {

struct PitchEvent {
	float hz, tonality;
	size_t offset;
};

struct PitchTrack {
	const float minHz, maxHz, bandwidthFactor;
	PitchTrack(float minHz=20, float maxHz=8000, float bandwidthFactor=0.5f) : minHz(minHz), maxHz(maxHz), bandwidthFactor(bandwidthFactor) {
		float hz = minHz;
		while (hz <= maxHz) {
			bandHz.push_back(hz);
			hz *= 1 + 0.15f*bandwidthFactor;
		}
		
		size_t bands = bandHz.size();
		inputIndex.resize(bands);
		inputIndexFrac.resize(bands);
		inputReal.resize(bands);
		inputImag.resize(bands);
		filterGain.resize(bands);
		filterTwistReal.resize(bands);
		filterTwistImag.resize(bands);
		stateReal.resize(bands);
		stateImag.resize(bands);
		bandEnergy.resize(bands);
	}

	void configure(float sampleRate, size_t maxBlockSize=0, size_t latencySamples=0, float overlap=4) {
		if (!latencySamples) latencySamples = size_t(sampleRate*0.002);
		
		size_t intervalSamples = latencySamples/overlap;
		stft.configure(1, 0, latencySamples, 0, intervalSamples);
		
		// Input bin indices
		size_t bands = bandHz.size();
		for (size_t b = 0; b < bands; ++b) {
			auto hz = bandHz[b];
			auto bin = stft.freqToBin(hz/sampleRate);
			inputIndex[b] = size_t(bin);
			inputIndexFrac[b] = bin - inputIndex[b];
		}
		
		// Complex filter: y[n] = twist*y[n - 1] + gain*x[n];
		float subRate = sampleRate/intervalSamples;
		auto barkScale = signalsmith::curves::Reciprocal<float>::barkScale();
		for (size_t b = 0; b < bands; ++b) {
			float bark = barkScale.inverse(bandHz[b]);
			float bandwidth = barkScale.dx(bark)*bandwidthFactor;
			auto logTwist = Complex{-bandwidth, bandHz[b]}*float(2*M_PI)/subRate;
			float absTwist = std::exp(logTwist.real());
			Complex twist = std::polar(absTwist, logTwist.imag());
			filterGain[b] = 1 - absTwist; // 0dB gain at peak
			filterTwistReal[b] = twist.real();
			filterTwistImag[b] = twist.imag();
		}
		
		events.reserve(maxBlockSize/intervalSamples + 1);
		events.resize(0);
	}
	
	const std::vector<PitchEvent> & updates() const {
		return events;
	}
	
	void reset() {
		events.resize(0);
		stft.reset();
		linear(stateReal, stateImag) = 0;
		linear(bandEnergy) = 0;
	}
	
	void block(const float *buffer, size_t length) {
		events.resize(0);
		for (size_t i = 0; i < length; ++i) {
			stft.writeInput(0, 1, buffer + i);
			stft.moveInput(1);
			if (stft.samplesSinceAnalysis() >= stft.defaultInterval()) {
				analyseBlock();
			}
		}
	}
	
	const std::vector<float> & debugHz() const {
		return bandHz;
	}
	const std::vector<float> & debugEnergy() const {
		return bandEnergy;
	}
	
private:
	using Complex = std::complex<float>;

	std::vector<PitchEvent> events;
	signalsmith::linear::DynamicSTFT<float> stft;
	signalsmith::linear::Linear linear;

	std::vector<float> bandHz;
	std::vector<size_t> inputIndex;
	std::vector<float> inputIndexFrac;
	std::vector<float> inputReal, inputImag;
	std::vector<float> filterGain, filterTwistReal, filterTwistImag;
	std::vector<float> stateReal, stateImag;
	std::vector<float> bandEnergy;

	void analyseBlock() {
		for (size_t s = 0; s < analyseSteps; ++s) {
			analyseBlock(s);
		}
	}

	static constexpr size_t analyseSteps = 100;
	void analyseBlock(size_t step) {
		size_t bands = bandHz.size();

		if (!step--) {
			stft.analyse();
		} else if (!step--) {
			auto *spectrum = stft.spectrum(0);
			for (size_t b = 0; b < bands; ++b) {
				Complex bin = spectrum[inputIndex[b]];
				Complex binNext = spectrum[inputIndex[b] + 1];
				bin += (binNext - bin)*inputIndexFrac[b];
				inputReal[b] = bin.real();
				inputImag[b] = bin.imag();
			}
		} else if (!step--) {
			auto input = linear(inputReal, inputImag);
			auto gain = linear(filterGain);
			auto twist = linear(filterTwistReal, filterTwistImag);
			auto state = linear(stateReal, stateImag);
			state = state*twist + input*gain;
			linear(bandEnergy) = state.norm();
		}
	}
};

}} // namespace
