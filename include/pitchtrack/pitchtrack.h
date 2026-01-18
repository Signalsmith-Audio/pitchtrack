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
	PitchTrack(float minHz=20, float maxHz=6000, float bandwidthFactor=0.25f) : minHz(minHz), maxHz(maxHz), bandwidthFactor(bandwidthFactor) {
		float hz = minHz;
		float bandFreqRatio = 1 + 0.15f*bandwidthFactor;
		while (hz <= maxHz) {
			bandHz.push_back(hz);
			hz *= bandFreqRatio;
		}

		bands = bandHz.size();
		inputIndex.resize(bands);
		inputIndexFrac.resize(bands);
		inputReal.resize(bands);
		inputImag.resize(bands);
		filterGain.resize(bands);
		filterTwistReal.resize(bands);
		filterTwistImag.resize(bands);
		stateReal.resize(bands);
		stateImag.resize(bands);
		state2Real.resize(bands);
		state2Imag.resize(bands);

		// Generate the scoring kernel
		size_t maxScoreHarmonic = 10;
		size_t maxScoreKernelIndex = 0;
		bandScoreKernel.resize(bands);
		auto addToKernel = [&](float index, float weight){
			size_t indexLow = std::floor(index);
			float indexFrac = index - indexLow;
			bandScoreKernel[indexLow] += weight*(1 - indexFrac);
			bandScoreKernel[indexLow + 1] += weight*indexFrac;
			maxScoreKernelIndex = std::max(maxScoreKernelIndex, indexLow + 1);
		};
		for (size_t h = 1; h < maxScoreHarmonic; ++h) {
			float weight = 1 - h/float(maxScoreHarmonic);

			addToKernel(std::log(2*h - 1)/std::log(bandFreqRatio), -weight);
			addToKernel(std::log(2*h)/std::log(bandFreqRatio), weight);
		}

		size_t paddedBands = energyFft.fastSizeAbove(bands + maxScoreKernelIndex);
		energyFft.resize(paddedBands);
		bandScoreKernel.resize(paddedBands);
		bandEnergy.resize(paddedBands);
		bandScore.resize(paddedBands);

		bandScoreSpectrum.resize(paddedBands/2);
		bandScoreKernelSpectrum.resize(paddedBands/2);

		energyFft.fft(bandScoreKernel.data(), bandScoreKernelSpectrum.data());
		linear(bandScoreKernelSpectrum) /= float(paddedBands);
		linear(bandScoreKernelSpectrum) = linear(bandScoreKernelSpectrum).conj();
	}
	
	void configure(float sampleRate, size_t maxBlockSize=0, size_t latencySamples=0, float overlap=4) {
		if (!latencySamples) latencySamples = size_t(sampleRate*0.002);
		
		size_t intervalSamples = latencySamples/overlap;
		stft.configure(1, 0, latencySamples, 0, intervalSamples);
		
		// Input bin indices
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
		reset();
	}

	size_t latencySamples() const {
		return stft.blockSamples();
	}

	const std::vector<PitchEvent> & updates() const {
		return events;
	}
	
	void reset() {
		events.resize(0);
		stft.reset();
		linear(stateReal, stateImag) = 0;
		linear(state2Real, state2Imag) = 0;
		linear(bandEnergy) = 0;
		linear(bandScore) = 0;
	}
	
	void block(const float *buffer, size_t length) {
		events.resize(0);
		for (size_t i = 0; i < length; ++i) {
			stft.writeInput(0, 1, buffer + i);
			stft.moveInput(1);
			if (stft.samplesSinceAnalysis() >= stft.defaultInterval()) {
				blockOffset = i;
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
	const std::vector<float> & debugScore() const {
		return bandScore;
	}
	const std::vector<float> & debugScoreKernel() const {
		return bandScoreKernel;
	}
	
private:
	using Complex = std::complex<float>;

	std::vector<PitchEvent> events;
	signalsmith::linear::DynamicSTFT<float> stft;
	signalsmith::linear::RealFFT<float> energyFft;
	signalsmith::linear::Linear linear;

	size_t bands;
	std::vector<float> bandHz;
	std::vector<size_t> inputIndex;
	std::vector<float> inputIndexFrac;
	std::vector<float> inputReal, inputImag;
	std::vector<float> filterGain, filterTwistReal, filterTwistImag;
	std::vector<float> stateReal, stateImag;
	std::vector<float> state2Real, state2Imag;
	std::vector<float> bandEnergy;

	std::vector<float> bandScoreKernel;
	std::vector<Complex> bandScoreKernelSpectrum;
	std::vector<Complex> bandScoreSpectrum;
	std::vector<float> bandScore;

	void analyseBlock() {
		for (size_t s = 0; s < analyseSteps; ++s) {
			analyseBlock(s);
		}
	}

	size_t blockOffset = 0;
	static constexpr size_t analyseSteps = 100;
	void analyseBlock(size_t step) {
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
			auto state2 = linear(state2Real, state2Imag);
			state = state*twist + input*gain;
			state2 = state2*twist + state*gain;
			
			auto energy = linear(bandEnergy.data(), bandHz.size()); // un-padded size
			energy = state2.norm();
			energy *= 1/float(stft.blockSamples()*stft.blockSamples());
		} else if (!step--) {
			energyFft.fft(bandEnergy.data(), bandScoreSpectrum.data());
			linear(bandScoreSpectrum) *= linear(bandScoreKernelSpectrum);
			energyFft.ifft(bandScoreSpectrum.data(), bandScore.data());
		} else if (!step--) {
			for (auto &v : bandScore) v = std::max<float>(v, 0);
//			for (size_t b = 0; b < bands; ++b) {
//				bandScore[b] = std::max<float>(bandScoreKernel[b], 0);
//			}
		} else if (!step--) {
			events.push_back({0, 0, blockOffset});
		}
	}
};

}} // namespace
