// SimplePMTResponse.hh
#pragma once

#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include <vector>
#include <algorithm>

// A tiny PMT response helper: QE probability and simple Gaussian TTS.
// You can construct with a constant QE and TTS or with wavelength-dependent arrays.
class SimplePMTResponse {
public:
    // constant QE (0..1), TTS sigma in seconds (default 2 ns)
    SimplePMTResponse(double QE = 0.21, double TTS_sigma = 2.0 * ns)
        : fQE(QE), fTTS_sigma(TTS_sigma) {}

    // wavelength-dependent constructor: wavelengths (nm) and QE (0..1)
    SimplePMTResponse(const std::vector<double>& wavelengths_nm,
                      const std::vector<double>& qe_vals,
                      double TTS_sigma = 2.0 * ns)
        : fTTS_sigma(TTS_sigma) {
        // store pairs sorted by wavelength ascending
        size_t n = std::min(wavelengths_nm.size(), qe_vals.size());
        for (size_t i = 0; i < n; ++i) {
            fWL.push_back(wavelengths_nm[i]);
            fQEvals.push_back(qe_vals[i]);
        }
        // ensure ascending
        for (size_t i = 1; i < fWL.size(); ++i) {
            if (fWL[i] < fWL[i-1]) {
                // simple fix: sort both arrays together
                std::vector<std::pair<double,double>> tmp;
                for (size_t j=0;j<fWL.size();++j) tmp.emplace_back(fWL[j], fQEvals[j]);
                std::sort(tmp.begin(), tmp.end());
                for (size_t j=0;j<fWL.size();++j) { fWL[j]=tmp[j].first; fQEvals[j]=tmp[j].second; }
                break;
            }
        }
    }

    // Decide if a photon of wavelength (nm) is detected
    bool IsDetected(double wavelength_nm = -1.0) const {
        double qe = fQE;
        if (!fWL.empty() && wavelength_nm > 0.0) qe = InterpolateQE(wavelength_nm);
        return G4UniformRand() < qe;
    }

    // smear time (returns smeared time in same units as input)
    double SmearTime(double trueTime) const {
        if (fTTS_sigma <= 0.) return trueTime;
        return G4RandGauss::shoot(trueTime, fTTS_sigma);
    }

    void SetConstantQE(double q) { fQE = q; }
    void SetTTSSigma(double tts) { fTTS_sigma = tts; }

private:
    double InterpolateQE(double wl_nm) const {
        if (fWL.empty()) return fQE;
        if (wl_nm <= fWL.front()) return fQEvals.front();
        if (wl_nm >= fWL.back()) return fQEvals.back();
        // linear search (small arrays)
        for (size_t i = 1; i < fWL.size(); ++i) {
            if (wl_nm <= fWL[i]) {
                double x0 = fWL[i-1], x1 = fWL[i];
                double y0 = fQEvals[i-1], y1 = fQEvals[i];
                double t = (wl_nm - x0) / (x1 - x0);
                return y0 + t * (y1 - y0);
            }
        }
        return fQEvals.back();
    }

    double fQE = 0.21;
    double fTTS_sigma = 2.0 * ns;
    std::vector<double> fWL;
    std::vector<double> fQEvals;
};
