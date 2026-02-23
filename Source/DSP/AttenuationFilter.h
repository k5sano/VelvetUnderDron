#pragma once

#include <cmath>

namespace DSP
{

/**
 * 1 次シェルビングフィルタ（周波数依存減衰）。
 * 各 FDN ディレイライン後段に配置し、
 * 低域と高域で異なる RT60 を実現する。
 *
 * H(z) = (b0 + b1 * z^-1) / (1 + a1 * z^-1)
 */
class AttenuationFilter
{
public:
    AttenuationFilter() = default;

    /**
     * フィルタ係数を設定。
     * @param gainLow  低域ループゲイン（0–1）
     * @param gainHigh 高域ループゲイン（0–1）
     * @param crossoverFreq クロスオーバー周波数 (Hz)
     * @param sampleRate サンプルレート (Hz)
     */
    void setCoefficients (float gainLow, float gainHigh,
                          float crossoverFreq, float sampleRate)
    {
        // Jot attenuation filter (Schlecht 2018 PhD thesis)
        // H(z) = (b0 + b1*z^-1) / (1 + a1*z^-1)
        // Design: H(1) = gainLow (DC), H(-1) = gainHigh (Nyquist)

        if (std::abs (gainLow - gainHigh) < 1.0e-6f)
        {
            b0 = gainLow;
            b1 = 0.0f;
            a1Coeff = 0.0f;
            return;
        }

        float wc = 3.14159265f * crossoverFreq / sampleRate;
        float t = std::tan (wc);

        // First-order allpass coefficient for crossover
        float a1 = (t - 1.0f) / (t + 1.0f);

        // Shelving coefficients
        b0 = (gainLow * (1.0f + a1) + gainHigh * (1.0f - a1)) * 0.5f;
        b1 = (gainLow * (1.0f + a1) - gainHigh * (1.0f - a1)) * 0.5f;
        a1Coeff = a1;
    }

    float process (float input)
    {
        float output = b0 * input + b1 * z1 - a1Coeff * zOut1;
        z1 = input;
        zOut1 = output;
        return output;
    }

    void reset()
    {
        z1 = 0.0f;
        zOut1 = 0.0f;
    }

private:
    float b0 = 1.0f;
    float b1 = 0.0f;
    float a1Coeff = 0.0f;
    float z1 = 0.0f;      // 入力遅延レジスタ
    float zOut1 = 0.0f;    // 出力遅延レジスタ
};

}  // namespace DSP
