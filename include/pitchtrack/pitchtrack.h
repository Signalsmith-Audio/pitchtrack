#pragma once

#include "signalsmith-linear/stft.h"
#include "signalsmith-linear/linear.h"

#include <vector>

namespace signalsmith { namespace pitchtrack {

struct PitchEvent {
	float hz, tonality;
	size_t offset;
	bool voiced;
};

struct PitchTrack {
	const float minHz, maxHz, bandwidthFactor;
	PitchTrack(float minHz=20, float maxHz=6000, float bandwidthFactor=0.25f) : minHz(minHz), maxHz(maxHz), bandwidthFactor(bandwidthFactor) {
		float hz = minHz;
		bandFreqRatio = 1 + 0.15f*bandwidthFactor;
		octaveBands = std::log(2.0f)/std::log(bandFreqRatio);
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
		filterWeights.resize(bands);
		filterTwistReal.resize(bands);
		filterTwistImag.resize(bands);
		stateReal.resize(bands);
		stateImag.resize(bands);
		state2Real.resize(bands);
		state2Imag.resize(bands);

		scoreKernelLength = std::log(maxScoreHarmonic/minScoreHarmonic)/std::log(bandFreqRatio);
		scoreKernelLeadIn = std::log(1/minScoreHarmonic)/std::log(bandFreqRatio);
		scoreKernelLeadOut = scoreKernelLength - scoreKernelLeadIn;
		
		paddedBands = energyFft.fastSizeAbove(bands + scoreKernelLength);
		energyFft.resize(paddedBands);
		bandScoreKernel.resize(paddedBands);
		bandEnergy.resize(paddedBands);
		bandScore.resize(paddedBands);
		bandScoreResult.resize(paddedBands);
		bandScoreCumulativeSum.resize(paddedBands);
		bandScoreSpectrum.resize(paddedBands/2);
		bandScoreKernelSpectrum.resize(paddedBands/2);

		// Generate the scoring kernel
		for (size_t b = 0; b < paddedBands; ++b) {
			float ratio = std::pow(bandFreqRatio, b);
			if (b > paddedBands - scoreKernelLeadIn) {
				ratio = std::pow(bandFreqRatio, -int(paddedBands - b));
			};
			if (ratio > minScoreHarmonic && ratio < maxScoreHarmonic) {
				bandScoreKernel[b] = scoreKernel(ratio);
			}
		}

		// normalise kernel to ~1 total energy
		float kernelSum2 = 0;
		for (auto &k : bandScoreKernel) kernelSum2 += k*k;
		linear(bandScoreKernel) /= std::sqrt(kernelSum2);
		scoreWeightBelow /= std::sqrt(kernelSum2);
		scoreWeightAbove /= std::sqrt(kernelSum2);

		// and store in the frequency domain
		energyFft.fft(bandScoreKernel.data(), bandScoreKernelSpectrum.data());
		bandScoreKernelSpectrum[0].imag(0);
		linear(bandScoreKernelSpectrum) /= float(paddedBands);
		linear(bandScoreKernelSpectrum) = linear(bandScoreKernelSpectrum).conj();
		
		setTracking();
	}
	
	void configure(float sampleRate, size_t maxBlockSize=0, size_t latencySamples=0, float overlap=2.5) {
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
		for (size_t b = 0; b < bands; ++b) {
			float bandwidth = bandHz[b]*(bandFreqRatio - 1);
			auto logTwist = Complex{-bandwidth, bandHz[b]}*float(2*M_PI)/subRate;
			float absTwist = std::exp(logTwist.real());
			Complex twist = std::polar(absTwist, logTwist.imag());
			filterGain[b] = (1 - absTwist); // 0dB peak
//			filterWeights[b] = std::pow(bandwidth, 2);
			filterWeights[b] = bandHz[b]/bandwidth; // reduce overlapping weight (if we're not constant Q)
			//filterWeights[b] *= 1 - std::pow(b/float(bands), 2.0f); // taper high end
			filterTwistReal[b] = twist.real();
			filterTwistImag[b] = twist.imag();
		}
		
		events.reserve(maxBlockSize/intervalSamples + 1);
		reset();
	}

	void setTracking(float tonalityOn=0.8f, float tonalityOff=0.5f, float stickiness=0.5f, float jumpRatio=0.1f) {
		thresholdTonalityOn = tonalityOn;
		thresholdTonalityOff = tonalityOff;
		scoreStickiness = stickiness;
		thresholdJumpRatio = jumpRatio;
	}

	size_t latencySamples() const {
		return stft.blockSamples();
	}

	const std::vector<PitchEvent> & updates() const {
		return events;
	}
	
	void reset() {
		events.resize(0);
		
		// Resonator state
		stft.reset();
		linear(stateReal, stateImag) = 0;
		linear(state2Real, state2Imag) = 0;
		linear(bandEnergy) = 0;
		linear(bandScore) = 0;

		// Tracking state
		prevScoreBand = 0;
		prevScore = 0;
		prevHz = -1;
		voiced = 0;
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
		return bandScoreResult;
	}
	const std::vector<float> & debugScoreRef() const {
		return bandScore;
	}
	const std::vector<float> & debugScoreKernel() const {
		return bandScoreKernel;
	}
	
private:
	using Complex = std::complex<float>;

	float thresholdTonalityOn = 0, thresholdTonalityOff = 0, scoreStickiness = 0, thresholdJumpRatio = 0;

	float bandFreqRatio = 0, octaveBands = 0;
	float prevScoreBand = 0, prevScore = 0, prevHz = -1;
	bool voiced = false;

	std::vector<PitchEvent> events;
	signalsmith::linear::DynamicSTFT<float> stft;
	signalsmith::linear::RealFFT<float> energyFft;
	signalsmith::linear::Linear linear;

	size_t bands, paddedBands;
	std::vector<float> bandHz;
	std::vector<size_t> inputIndex;
	std::vector<float> inputIndexFrac;
	std::vector<float> inputReal, inputImag;
	std::vector<float> filterGain, filterWeights, filterTwistReal, filterTwistImag;
	std::vector<float> stateReal, stateImag;
	std::vector<float> state2Real, state2Imag;
	std::vector<float> bandEnergy;

	size_t scoreKernelLength, scoreKernelLeadIn, scoreKernelLeadOut;
	std::vector<float> bandScoreKernel;
	std::vector<Complex> bandScoreKernelSpectrum;
	std::vector<Complex> bandScoreSpectrum;
	std::vector<float> bandScore, bandScoreResult, bandScoreCumulativeSum;

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
			
			auto energy = linear(bandEnergy.data(), bands); // un-padded size
			energy = state2.norm();
		} else if (!step--) {
			linear(bandScore) = linear(bandEnergy).sqrt().sqrt().sqrt();
			linear(bandScore.data(), bands) *= linear(filterWeights);

			energyFft.fft(bandScore.data(), bandScoreSpectrum.data());
			linear(bandScoreSpectrum) *= linear(bandScoreKernelSpectrum);
			energyFft.ifft(bandScoreSpectrum.data(), bandScoreResult.data());
		} else if (!step--) {
			if (scoreWeightBelow != 0 || scoreWeightAbove != 0) {
				// Cumulative sum
				float sum = 1e-30f;
				for (size_t b = 0; b < bands; ++b) {
					sum += bandScore[b];
					bandScoreCumulativeSum[b] = sum;
				}
				if (scoreWeightBelow != 0) {
					for (size_t b = scoreKernelLeadIn; b < bands; ++b) {
						bandScoreResult[b] += scoreWeightBelow*bandScoreCumulativeSum[b - scoreKernelLeadIn];
					}
				}
				if (scoreWeightAbove != 0) {
					auto endScore = bandScoreCumulativeSum[bands - 1];
					for (size_t b = 0; b + scoreKernelLeadOut < bands; ++b) {
						auto remainingScore = endScore - bandScoreCumulativeSum[b + scoreKernelLeadOut];
						bandScoreResult[b] += scoreWeightAbove*remainingScore;
					}
				}
			}
		} else if (!step--) {
			if (voiced && prevScore > 0) {
				for (size_t b = 0; b < bands; ++b) {
					float d = (b - prevScoreBand)/octaveBands;
					bandScoreResult[b] += scoreStickiness*prevScore/(1 + d*d);
				}
			}
		} else if (!step--) {
			size_t bestB = 0, metricBLimit = 0;
			float bestFracB = 0, bestScore = -1e30, nextBestScore = bestScore;
			for (size_t b = 1; b + 1 < bands; ++b) {
				auto curr = bandScoreResult[b];
				if (curr > nextBestScore) {
					auto prev = bandScoreResult[b - 1], next = bandScoreResult[b + 1];
					if (curr >= prev && curr > next) {
						// It's a peak - is it a new best one?
						if (curr > bestScore) {
							nextBestScore = bestScore;
							bestScore = curr;
							bestB = b;
							bestFracB = 0.5f*(next - prev)/(2*curr - prev - next + 1e-30f);
						} else {
							nextBestScore = curr;
						}
					}
				}
			}
			prevScore = bestScore;
			prevScoreBand = bestB + bestFracB;
		} else if (!step--) {
			if (prevScore > 0) {
				size_t bestB = std::floor(prevScoreBand);
				float bestFracB = prevScoreBand - bestB;
				float hz = bandHz[bestB]*std::pow(bandFreqRatio, bestFracB);
				float onEnergy = 1e-30f, allEnergy = 1e-30f;
				for (size_t b = 0; b < bands; ++b) {
					float h = bandHz[b]/hz;
					if (h > 10.5) break;
					float e = bandEnergy[b]*filterWeights[b];
					allEnergy += e;
					onEnergy += e*0.5f*(1 + std::cos(h*float(2*M_PI)));
				}
				float tonality = onEnergy/allEnergy;
				if (tonality > thresholdTonalityOn) {
					voiced = true;
				} else if (tonality < thresholdTonalityOff) {
					voiced = false;
				} else {
					auto jump = std::abs(hz/prevHz - 1);
					if (jump > thresholdJumpRatio) {
						// The tonality is in the middle ground, but we jumped frequency too suddenly
						voiced = false;
					}
				}
				prevHz = hz;
				events.push_back({hz, tonality, blockOffset, voiced});
			} else {
				voiced = false;
				prevHz = -1;
				events.push_back({0, 0, blockOffset, false});
			}
		}
	}
	
	static constexpr float maxScoreHarmonic = 10.25f, minScoreHarmonic = 0.25f;
	float scoreWeightBelow = -0.01f, scoreWeightAbove = 0;
	static float scoreKernel(float harmonic) {
		auto ratio = harmonic/maxScoreHarmonic;
		auto k = std::cos(float(2*M_PI)*harmonic);
		if (harmonic < 1) k *= harmonic*harmonic*(3 - 2*harmonic);
		return k*(1 - ratio*ratio*ratio);
	}
};

}} // namespace
