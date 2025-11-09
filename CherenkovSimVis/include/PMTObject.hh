#pragma once

#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <string>

namespace CherenkovSim::PMT {

using real = double;

// Abstract base exposing what DetectorConstruction needs
class PMTBase {
public:
  virtual ~PMTBase() = default;

  // geometry
  virtual real GetExposeHeight_m() const = 0;   // expose height (m)
  virtual real GetRadius_m() const = 0;         // radius at expose height (m)
  virtual real GetGlassThickness_m() const = 0; // glass thickness (m)

  // derived geometry helper: computes sphere radius given expose & radius
  // R = (h^2 + r^2) / (2*h)
  real ComputeSphereRadius_m() const {
    const real h = GetExposeHeight_m();
    const real r = GetRadius_m();
    assert(h > 0.0 && r >= 0.0);
    return (h*h + r*r) / (2.0 * h);
  }

  // Simple QE support: wavelength (nm) & QE (0..1) arrays
  // default empty: derived classes should override if they provide QE.
  virtual std::vector<real> GetQEWavelength_nm() const { return {}; }
  virtual std::vector<real> GetQE() const { return {}; }

  // convenience: linear-interpolate QE at a given wavelength (nm)
  // returns 0 if out-of-range or no data
  real QE_at_nm(real wl_nm) const {
    auto wl = GetQEWavelength_nm();
    auto qe = GetQE();
    if (wl.empty() || qe.empty() || wl.size() != qe.size()) return 0.0;
    if (wl.size() == 1) return qe[0];
    // enforce ascending wavelengths
    if (wl.front() > wl.back()) {
      std::reverse(wl.begin(), wl.end());
      std::reverse(qe.begin(), qe.end());
    }
    if (wl_nm <= wl.front()) return qe.front();
    if (wl_nm >= wl.back()) return qe.back();
    // find segment
    auto it = std::upper_bound(wl.begin(), wl.end(), wl_nm);
    size_t idx = std::distance(wl.begin(), it);
    if (idx == 0) return qe.front();
    // interpolate between idx-1 and idx
    real x0 = wl[idx-1], x1 = wl[idx];
    real y0 = qe[idx-1], y1 = qe[idx];
    real t = (wl_nm - x0) / (x1 - x0);
    return y0 + t * (y1 - y0);
  }

  // optional: collection efficiency vs angle (deg). default flat 100%
  // angle array (deg) and efficiency array (percent)
  virtual std::vector<real> GetCollectionEfficiencyAngle_deg() const {
    return {0., 90.}; // 0..90 deg
  }
  virtual std::vector<real> GetCollectionEfficiency_percent() const {
    return {100., 100.}; // flat 100%
  }

  // convenience to get collection efficiency (0..1) at angle (deg)
  real CollectionEfficiencyAt_deg(real angle_deg) const {
    auto ang = GetCollectionEfficiencyAngle_deg();
    auto ce  = GetCollectionEfficiency_percent();
    if (ang.empty() || ce.empty() || ang.size() != ce.size()) return 1.0;
    // linear interpolate similar to QE_at_nm
    if (ang.front() > ang.back()) {
      std::reverse(ang.begin(), ang.end());
      std::reverse(ce.begin(), ce.end());
    }
    if (angle_deg <= ang.front()) return ce.front() / 100.0;
    if (angle_deg >= ang.back()) return ce.back() / 100.0;
    auto it = std::upper_bound(ang.begin(), ang.end(), angle_deg);
    size_t idx = std::distance(ang.begin(), it);
    real x0 = ang[idx-1], x1 = ang[idx];
    real y0 = ce[idx-1], y1 = ce[idx];
    real t = (angle_deg - x0) / (x1 - x0);
    return (y0 + t*(y1-y0)) / 100.0;
  }
};

///////////////////////////////////////////////////////////////////////////////
// PMT20inch: minimal implementation with realistic numbers (from WCSim)
// - all lengths in meters
// - QE provided as a small example (wavelengths in nm)
// - you can expand/replace QE later (file or higher-resolution array)
///////////////////////////////////////////////////////////////////////////////
class PMT20inch : public PMTBase {
public:
  PMT20inch() = default;
  ~PMT20inch() override = default;

  // geometry values taken from WCSim PMT20inch.cc (converted to SI)
  real GetExposeHeight_m() const override { return 0.18; }      // 0.18 m
  real GetRadius_m() const override { return 0.254; }          // 0.254 m (20" / 2)
  real GetGlassThickness_m() const override { return 0.004; }  // 0.4 cm = 0.004 m

  // small QE sample (same shape as WCSim example; numbers are fractional QE)
  std::vector<real> GetQEWavelength_nm() const override {
    return {280., 300., 320., 340., 360., 380., 400., 420., 440., 460.,
            480., 500., 520., 540., 560., 580., 600., 620., 640., 660.};
  }

  std::vector<real> GetQE() const override {
    return {0.00, 0.0139, 0.0854, 0.169, 0.203, 0.206, 0.211, 0.202, 0.188,
            0.167, 0.140, 0.116, 0.0806, 0.0432, 0.0265, 0.0146, 0.00756,
            0.00508, 0.00158, 0.00};
  }

  // collection efficiency: default WCSim-like flat 100% but you can customize
  std::vector<real> GetCollectionEfficiencyAngle_deg() const override {
    return {0., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
  }
  std::vector<real> GetCollectionEfficiency_percent() const override {
    return {100., 100., 100., 100., 100., 100., 100., 100., 100., 100.};
  }

  // convenience accessor for maximum QE
  real GetMaxQE() const { return 0.211; }
};

} // namespace CherenkovSim::PMT
