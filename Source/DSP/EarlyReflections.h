#pragma once

#include "DSP/VelvetNoise.h"

namespace DSP
{

/**
 * OVN early reflections (Layer 1).
 * Sparse FIR preserves string instrument transients.
 * Linear processing: no oversampling required.
 */
class EarlyReflections
{
public:
    EarlyReflections() = default;

    void prepare (double sampleRate, int maxBlockSize, uint32_t seed)
    {
        juce::ignoreUnused (maxBlockSize);
        sr = sampleRate;
        ovn.generate (sampleRate, 30.0f, 2000.0f, seed);
    }

    void process (const float* input, float* output,
                  int numSamples, float gain)
    {
        ovn.convolve (input, output, numSamples, gain);
    }

    void reset()
    {
        // Ring buffer state lives inside VelvetNoise
    }

private:
    VelvetNoise ovn;
    double sr = 44100.0;
};

}  // namespace DSP
