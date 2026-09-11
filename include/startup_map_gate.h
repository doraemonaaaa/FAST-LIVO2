#pragma once
#include <cmath>
#include <initializer_list>
#include <stdexcept>

// One-shot admission gate. It never asserts zero motion or resets the EKF.
struct StartupMapGate
{
  bool enabled = false;
  double minimum_seconds = 2, maximum_seconds = 30, maximum_gap_seconds = .5;
  int consecutive_frames = 20, minimum_matches = 100;
  double minimum_match_ratio = .2, maximum_residual = .15;
  double maximum_translation = .10, maximum_rotation = .02;
  double maximum_velocity = .15, maximum_gravity = .05;
  enum Phase { Warmup, Recheck, Open } phase = Warmup;
  enum Event { None, Rebuild, Admit };
  double start = -1, previous_time = -1;
  int stable = 0;
  bool isOpen() const { return !enabled || phase == Open; }
  void validate() const
  {
    for (double x : {minimum_seconds, maximum_seconds, maximum_gap_seconds,
         minimum_match_ratio, maximum_residual, maximum_translation,
         maximum_rotation, maximum_velocity, maximum_gravity})
      if (!std::isfinite(x) || x <= 0) throw std::runtime_error("Invalid startup_map_gate threshold");
    if (maximum_seconds <= minimum_seconds || minimum_match_ratio > 1 ||
        consecutive_frames < 1 || minimum_matches < 1)
      throw std::runtime_error("Invalid startup_map_gate configuration");
  }
  Event update(double time, int matches, int points, double residual,
               double translation, double rotation, double velocity, double gravity)
  {
    if (isOpen()) return None;
    if (!std::isfinite(time) || (previous_time >= 0 && time <= previous_time))
      throw std::runtime_error("Non-increasing startup map gate timestamp");
    if (start < 0) start = time;
    bool good = matches >= minimum_matches && points > 0 &&
                double(matches) / points >= minimum_match_ratio;
    const double values[] = {residual, translation, rotation, velocity, gravity};
    const double limits[] = {maximum_residual, maximum_translation, maximum_rotation,
                             maximum_velocity, maximum_gravity};
    for (int i = 0; i < 5; ++i)
      good = good && std::isfinite(values[i]) && values[i] >= 0 && values[i] <= limits[i];
    if (previous_time >= 0 && time - previous_time > maximum_gap_seconds) stable = 0;
    previous_time = time;
    stable = good ? stable + 1 : 0;
    if (time - start > maximum_seconds)
      throw std::runtime_error("Startup map gate timed out: no stable map admission; inspect startup_map_gate.csv");
    if (time - start < minimum_seconds || stable < consecutive_frames) return None;
    stable = 0;
    if (phase == Warmup) { phase = Recheck; return Rebuild; }
    phase = Open;
    return Admit;
  }
};
