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

		float maxScoreHarmonic = 20.25f, minScoreHarmonic = 0.5f;
		scoreKernelLength = std::log(maxScoreHarmonic/minScoreHarmonic)/std::log(bandFreqRatio);
		scoreKernelLeadIn = std::log(1/minScoreHarmonic)/std::log(bandFreqRatio);
		
		paddedBands = energyFft.fastSizeAbove(bands + scoreKernelLength);
		energyFft.resize(paddedBands);
		bandScoreKernel.resize(paddedBands);
		bandEnergy.resize(paddedBands);
		bandScore.resize(paddedBands);
		bandScoreRef.resize(paddedBands);
		bandScoreSpectrum.resize(paddedBands/2);
		bandScoreKernelSpectrum.resize(paddedBands/2);

		// Generate the scoring kernel
		for (size_t b = 0; b < paddedBands; ++b) {
			float ratio = std::pow(bandFreqRatio, b);
			if (b > paddedBands - scoreKernelLeadIn) {
				ratio = std::pow(bandFreqRatio, -int(paddedBands - b));
			};
			if (ratio > minScoreHarmonic && ratio < maxScoreHarmonic) {
				float w = std::cos(float(2*M_PI)*ratio);
				bandScoreKernel[b] = w;
			}
		}
		// and store in the frequency domain
		energyFft.fft(bandScoreKernel.data(), bandScoreKernelSpectrum.data());
		bandScoreKernelSpectrum[0].imag(0);
		linear(bandScoreKernelSpectrum) /= float(paddedBands);
		linear(bandScoreKernelSpectrum) = linear(bandScoreKernelSpectrum).conj();
	}
	
	void configure(float sampleRate, size_t maxBlockSize=0, size_t latencySamples=0, float overlap=4) {
		if (!latencySamples) latencySamples = size_t(sampleRate*0.002);
		
		size_t intervalSamples = latencySamples/overlap;
		stft.configure(1, 0, latencySamples, 0, intervalSamples);
		linear(stft.analysisWindow(), stft.blockSamples()) /= float(stft.blockSamples());
		
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
			filterGain[b] = (1 - absTwist); // 0dB peak
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
	const std::vector<float> & debugScoreRef() {
		return bandScoreRef;
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

	size_t bands, paddedBands;
	std::vector<float> bandHz;
	std::vector<size_t> inputIndex;
	std::vector<float> inputIndexFrac;
	std::vector<float> inputReal, inputImag;
	std::vector<float> filterGain, filterTwistReal, filterTwistImag;
	std::vector<float> stateReal, stateImag;
	std::vector<float> state2Real, state2Imag;
	std::vector<float> bandEnergy;

	size_t scoreKernelLength, scoreKernelLeadIn;
	std::vector<float> bandScoreKernel;
	std::vector<Complex> bandScoreKernelSpectrum;
	std::vector<Complex> bandScoreSpectrum;
	std::vector<float> bandScore, bandScoreRef;

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
//energy = energy.sqrt();
//		} else if (!step--) {
//			float slew = 1/(1.0f + scoreKernelLeadIn*0.05f);
//			float e = 0;
//			for (size_t b = 0; b < paddedBands; ++b) {
//				bandScore[b] = e += (bandEnergy[b] - e)*slew;
//			}
//			for (size_t b = paddedBands; b-- > 0;) {
//				bandScore[b] = e += (bandScore[b] - e)*slew;
//			}
//			linear(bandEnergy) -= linear(bandScore);
//			linear(bandEnergy) = linear.max(linear(bandEnergy), 0);
//		} else if (!step--) {
//			for (size_t b = 1; b < bands - 1; ++b) {
//				auto prev = bandEnergy[b - 1], next = bandEnergy[b + 1];
//				if (bandEnergy[b] > prev && bandEnergy[b] >= next) {
//					bandEnergy[b] *= 100;
//				}
//			}
		} else if (!step--) {
			energyFft.fft(bandEnergy.data(), bandScoreSpectrum.data());
			linear(bandScoreSpectrum) *= linear(bandScoreKernelSpectrum);
			energyFft.ifft(bandScoreSpectrum.data(), bandScore.data());
		} else if (!step--) {
			// Cumulative sum
			float sum = 1e-30f;
			for (size_t b = 0; b < bands; ++b) {
				sum += bandEnergy[b];
				bandScoreRef[b] = sum;
			}
			linear(bandScoreRef) /= sum;
			linear(bandScore) /= sum;
			for (auto &v : bandScore) v = std::max<float>(v, 0);
		} else if (!step--) {
			events.push_back({0, 0, blockOffset});
		}
	}
};

}} // namespace
