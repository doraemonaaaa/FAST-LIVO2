#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace livo {

// A shift adds to the already configured IMU tags. The effective callback
// parameter is configured_imu_time_offset_s - stamp_shift_s.
struct OffsetMatchScore {
  size_t scans = 0, valid_scans = 0, points = 0, matches = 0;
  double squared_error = 0, capped_squared_error = 0;
  double ratio() const { return points ? double(matches)/points : 0; }
  double rms() const {
    return matches ? std::sqrt(squared_error/matches)
                   : std::numeric_limits<double>::infinity();
  }
  double capped_rms() const {
    return points ? std::sqrt(capped_squared_error/points)
                  : std::numeric_limits<double>::infinity();
  }
  bool valid() const {
    return scans && valid_scans == scans && matches >= 100 &&
        matches <= points && std::isfinite(rms()) && std::isfinite(capped_rms());
  }
  OffsetMatchScore& operator+=(const OffsetMatchScore& other) {
    scans += other.scans; valid_scans += other.valid_scans;
    points += other.points; matches += other.matches;
    squared_error += other.squared_error;
    capped_squared_error += other.capped_squared_error;
    return *this;
  }
};

struct OffsetCandidateScore {
  double stamp_shift_s = 0;
  OffsetMatchScore train, validation;
  OffsetMatchScore total() const { auto out = train; out += validation; return out; }
};

struct ImuTimeOffsetResult {
  bool requested = false, estimated = false, applied = false;
  const char* reason = "disabled";
  double stamp_shift_s = 0, candidate_shift_s = 0;
  double validation_shift_s = 0, window_s = 0;
  size_t scans = 0, observations = 0, validation_observations = 0;
  OffsetMatchScore before, after, validation_before, validation_after;
  std::vector<OffsetCandidateScore> candidates;
  double window_begin_s = 0, window_end_s = 0;
};

// Candidate scores are measured AFTER IMU deskew and rigid registration, with
// frozen reference maps, query point identities, priors and distance gates.
// This is a geometric selection score, not a clock-ground-truth measurement.
class ImuTimeOffsetEstimator {
 public:
  static constexpr double kWindowSeconds = 3.;
  // Absolute imu_time_offset_s parameter range, not shifts added to IMU tags.
  static constexpr double kMinimumOffset = 0.;
  static constexpr double kMaximumOffset = .030;
  static constexpr double kGridStep = .002;
  static constexpr double kMatchDistance = .05;
  static constexpr double kResidualCap = .30;

  static bool InSearchRange(double offset) {
    return std::isfinite(offset) && offset >= kMinimumOffset-1e-12 &&
           offset <= kMaximumOffset+1e-12;
  }

  static std::vector<double> CandidateShifts(double configured_offset_s) {
    if (!std::isfinite(configured_offset_s))
      throw std::invalid_argument("nonfinite configured IMU time offset");
    std::vector<double> shifts;
    const int steps = int(std::lround((kMaximumOffset-kMinimumOffset)/kGridStep));
    for (int i = 0; i <= steps; ++i)
      shifts.push_back(configured_offset_s-(kMinimumOffset+i*kGridStep));
    // Keep the configured timing for before/after diagnostics. If outside the
    // search range it is scored, but may never be selected as the optimum.
    if (std::none_of(shifts.begin(),shifts.end(),[](double s) { return std::abs(s)<1e-12; }))
      shifts.push_back(0.);
    return shifts;
  }

  // Mean truncated quadratic reward over ALL query points. Unmatched points
  // contribute zero; a match at distance r contributes 1-r^2/tau^2.
  // This equals ratio * (1 - RMS^2/tau^2), with the existing match gate tau.
  static double BalancedScore(const OffsetMatchScore& score) {
    if (!score.valid()) return -std::numeric_limits<double>::infinity();
    return (double(score.matches)-score.squared_error/(kMatchDistance*kMatchDistance)) /
           double(score.points);
  }

  static Eigen::Matrix3d Exp(const Eigen::Vector3d& v) {
    const double angle = v.norm();
    if (angle < 1e-12) return Eigen::Matrix3d::Identity();
    return Eigen::AngleAxisd(angle, v/angle).toRotationMatrix();
  }
  static Eigen::Vector3d Log(const Eigen::Matrix3d& R) {
    const Eigen::AngleAxisd aa(R);
    return aa.angle()*aa.axis();
  }

  static ImuTimeOffsetResult Estimate(std::vector<OffsetCandidateScore> scores,
                                    double configured_offset_s = 0.) {
    ImuTimeOffsetResult out; out.requested = true;
    out.candidates = std::move(scores);
    auto reject = [&](const char* reason) { out.reason = reason; return out; };
    const OffsetCandidateScore *baseline = nullptr, *best = nullptr, *validation = nullptr;
    auto better = [](const OffsetCandidateScore& a, const OffsetCandidateScore& b,
                     bool heldout = false) {
      const auto x = heldout ? a.validation : a.total();
      const auto y = heldout ? b.validation : b.total();
      const double sx = BalancedScore(x), sy = BalancedScore(y);
      if (std::abs(sx-sy) > 1e-12) return sx > sy;
      if (std::abs(x.rms()-y.rms()) > 1e-12) return x.rms() < y.rms();
      // Exact ties retain the smallest change from the configured timing.
      return std::abs(a.stamp_shift_s) < std::abs(b.stamp_shift_s);
    };
    for (const auto& c : out.candidates) {
      if (!std::isfinite(c.stamp_shift_s)) return reject("invalid_candidate_shift");
      if (std::abs(c.stamp_shift_s) < 1e-12) baseline = &c;
      if (!InSearchRange(configured_offset_s-c.stamp_shift_s)) continue;
      if (c.total().valid() && (!best || better(c,*best))) best = &c;
      // Temporal subsets are diagnostics only; they do not veto selection.
      if (c.validation.valid() && (!validation || better(c,*validation,true))) validation = &c;
    }
    if (!baseline) return reject("missing_baseline_candidate");
    out.before = baseline->total();
    out.validation_before = baseline->validation;
    out.observations = out.before.scans;
    out.validation_observations = baseline->validation.scans;
    if (!out.before.scans || !out.before.points) return reject("no_query_scans");
    // Never reward a candidate for dropping a hard scan or changing sampling.
    for (const auto& c : out.candidates)
      if (c.train.scans != baseline->train.scans || c.train.points != baseline->train.points ||
          c.validation.scans != baseline->validation.scans ||
          c.validation.points != baseline->validation.points)
        return reject("candidate_support_mismatch");
    if (!best) return reject("no_successful_candidate_registration");
    out.candidate_shift_s = best->stamp_shift_s;
    out.after = best->total();
    out.validation_after = best->validation;
    if (validation) out.validation_shift_s = validation->stamp_shift_s;
    // Select the maximum balanced score directly. There are no peak-width,
    // improvement, boundary, temporal-agreement or residual acceptance gates.
    out.estimated = true;
    out.applied = std::abs(best->stamp_shift_s) > 1e-12;
    out.stamp_shift_s = out.applied ? best->stamp_shift_s : 0;
    out.reason = out.applied ? "maximum_balanced_score" : "configured_offset_has_maximum_balanced_score";
    return out;
  }
};
}  // namespace livo
