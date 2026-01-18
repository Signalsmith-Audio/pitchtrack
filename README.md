# Monophonic Pitch Tracker

```cpp
#include "pitchtrack/pitchtrack.h"

signalsmith::pitchtrack::PitchTrack tracker;
tracker.config(sampleRate, maxBlockLength=0, latencySamples=<2ms>, overlap=4);
```

The `maxBlockLength` is optional, but providing it guarantees no allocations.  Latency shouldn't be greater than 10ms, but higher values use less CPU.  Reduced overlap similarly reduces CPU, because it processes FFT blocks less often.

You pass the tracker a block of audio, but it doesn't return you anything.  Instead, you can then inspect the list of timed pitch updates:

```
const float *floatSamples;
tracker.block(floatSamples);

for (auto &update : tracker.events()) {...}
```

These update events have:
* `.hz` (Hz)
* `.tonality` (0-1), where `0` means it's definitely noise.
* `.offset` (`size_t`) relative to the start of the block, not including the analysis latency.

This event array is cleared by a call to `tracker.reset()`, which also clears the internal buffers/etc.

## Architecture

For simplicity (and speed), it uses single-channel, 32-bit float.  There are four steps to the processing:

### 1. Split into sub-bands using STFT

This is done using the Linear library, so set `SIGNALSMITH_USE_...` so it uses Accelerate/PFFT/IPP/etc.

### 2. Log-spaced complex resonators

These are complex 1-pole filters run at the downsampled rate.  Although the centre frequencies are log-spaced, the widths are based on the Bark curve (critical bands).

### 3. Spectral estimator

This estimates the instantaneous pitch (and confidence), by comparing energy near harmonics against total energy.

### 4. Cleanup

The estimator produces multiple options, and this stage has some heuristics to produce smoother results in the case of noise/uncertainty.

## License and Support

The code is [MIT Licensed](LICENSE.txt).  The `test/` directory includes an audio file from the [vocadito data set](https://zenodo.org/records/5578807) which is released under [CCA-4.0](https://creativecommons.org/licenses/by/4.0/legalcode), but this is not included in any binary builds.
