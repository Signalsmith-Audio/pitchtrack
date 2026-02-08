# Monophonic Pitch Tracker

```cpp
#include "signalsmith-pitchtrack/pitchtrack.h"

signalsmith::pitchtrack::PitchTrack tracker;
tracker.config(sampleRate, maxBlockLength=0, latencySamples=<2ms>, overlap=<2.5>);
```

The `maxBlockLength` is optional, but providing it guarantees no allocations.  Latency is specified in samples, and shouldn't be greater than 10ms.  Higher latency uses less CPU, and increased overlap uses more, but I think the default configuration is pretty reasonable.

You pass the tracker a block of audio, but it doesn't give you any output.  Instead, you then check the list of timed pitch updates:

```
const float *floatSamples;
tracker.analyse(floatSamples, length);

for (auto &update : tracker.events()) {...}
```

These update events (`signalsmith::pitchtrack::PitchEvent`) have:
* `.voiced` (`bool`) - if this is `false`, then the other results should be less trusted
* `.hz` (Hz)
* `.energy` - square-root this to get an RMS (amplitude) level
* `.tonality` (0-1), where `0` means it's classified as noise, `1` means it's a note.
* `.tonalityRaw` (0-1), a raw energy ratio which never actually reaches 0-1 except in silence.
* `.offset` (`size_t`) relative to the start of the block, not including the analysis latency.

This event array is replaced by the next `.analyse()` call, and cleared by a call to `tracker.reset()` (which also clears the internal buffers/etc.)

## Architecture

For simplicity (and speed), it uses single-channel, 32-bit float.  There are four steps to the processing:

### 1. Split into sub-bands using STFT

This is done using the Linear library, so set [`SIGNALSMITH_USE_...`](https://github.com/Signalsmith-Audio/linear?tab=readme-ov-file#building) if you want Accelerate/PFFT/IPP/etc.

### 2. Log-spaced complex resonators

The resonators are complex Butterworth filters (roughly equivalent to 4th-order bandpass) run at the downsampled rate.

### 3. Spectral estimator

This estimates the instantaneous pitch (and confidence), by comparing energy near harmonics against total energy.

### 4. Cleanup

The estimator produces multiple peaks, and this stage has some heuristics to produce smoother results in the case of noise/uncertainty.

## License, testing and support

The code is [MIT Licensed](LICENSE.txt), and there may be [some quirks](SUPPORT.txt) since this is a fairly new project.

The `test/` directory includes a smoke-test which runs the tracker over an audio file from the [vocadito data set](https://zenodo.org/records/5578807).
